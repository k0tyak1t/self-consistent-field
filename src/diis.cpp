// ----------------------------------------------------------- //
//             Suffix _diis_ means DIIS stage                  //
// members without _diis_ suffix realtes to the pre-DIIS stage //
// ----------------------------------------------------------- //

#include "diis.h"
#include "matrix.h"

#include <lapacke.h>

#include <assert.h>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace { // temporary implementation
Matrix expand(Buffer<Matrix> buffer, std::vector<double> coefs) {
  auto result = Matrix::zero_like(buffer[0]);
  auto bi = buffer.begin(), be = buffer.end();
  auto ci = coefs.begin(), ce = coefs.end();
  for (; bi != be; ++bi, ++ci)
    result += (*bi) * (*ci);
  return result;
}
} // namespace

DIIS::DIIS(MO &mo, StandardMatrices &std_m)
    : SCF(mo, std_m), diis_size(DEFAULT_DIIS_SIZE),
      error(Matrix::zero(diis_size)),
      extended_diis_product(Matrix{diis_size + 1}),
      fock_buffer(diis_size, Matrix::zero(diis_size)),
      error_buffer(diis_size, Matrix::zero(diis_size)),
      density_buffer(diis_size, Matrix::zero(diis_size)),
      diis_coefs(diis_size + 1) {

  error = Matrix::zero(diis_size);

  for (std::size_t i = 0; i < diis_size; ++i) {
    extended_diis_product(diis_size, i) = 1.0;
    extended_diis_product(i, diis_size) = 1.0;
  }
}

void DIIS::update_error() {
  error = fock * density * std_m.S - std_m.S * density * fock;
}

void DIIS::update_extended_error_product() {
  for (auto i = 0u; i < diis_size; ++i)
    for (auto j = 0u; j < diis_size; ++j)
      extended_diis_product(i, j) =
          Matrix::dot(error_buffer[i], error_buffer[j]);

#ifndef NDEBUG
  std::ofstream log_stream;
  log_stream.open("logs/diis_Bmatrix.log", std::ios::app);
  log_stream << "\n\n\n\n" << extended_diis_product << std::endl;
  log_stream.close();
#endif
}

// solving system: extended_diis_product * diis_coefs = [0, ..., 0, -1]
void DIIS::update_diis_coefs() {
  int N = diis_size + 1;
  int NRHS = 1;
  std::vector<int> IPIV(N);

  // B matrix is extended with [1 ... 1 0] in constructor
  // Initialize RHS: [0 ... 0 1]
  diis_coefs = std::vector<double>(diis_size + 1, 0.0);
  diis_coefs[diis_size] = 1.0;

  Matrix b_copy{extended_diis_product};
  int INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, NRHS, b_copy.data(), N,
                           IPIV.data(), diis_coefs.data(), N);

  if (INFO != 0)
    throw std::runtime_error("Failed to solve DIIS system for coefficients. "
                             "INFO = " +
                             std::to_string(INFO) +
                             ", matrix size = " + std::to_string(N));

#ifdef NORMALIZE // normalization
  double sum{};
  for (std::size_t i = 0; i < diis_size; ++i)
    sum += diis_coefs[i];

  for (std::size_t i = 0; i < diis_size; ++i)
    diis_coefs[i] /= sum;
#endif

#ifndef NDEBUG
  std::cout << "DIIS coefs: ";
  for (double c : diis_coefs)
    std::cout << " " << c;
  std::cout << std::endl;
#endif
}

void DIIS::update_diis_error() {
#ifndef NDEBUG
  std::cout << "[DIIS]: DIIS stage error matrix update...\n";
#endif
  error = Matrix::zero_like(error);
  for (std::size_t k = 0; k < diis_size; ++k)
    error += error_buffer[k] * diis_coefs[k];

#ifndef NDEBUG
  std::cout << "[DIIS]: DIIS stage error matrix updated.\n";
#endif
}

void DIIS::update_diis_fock() {
#ifndef NDEBUG
  std::cout << "[DIIS]: DIIS stage fock matrix update...\n";
#endif
  fock = Matrix::zero_like(fock);
  for (std::size_t k = 0; k < diis_size; ++k)
    fock += fock_buffer[k] * diis_coefs[k];

#ifndef NDEBUG
  std::cout << "[DIIS]: DIIS stage fock matrix updated.\n";
#endif
}

void DIIS::update_diis_density() {
  density = Matrix::zero_like(density);
  for (auto i = 0u; i < diis_size; ++i)
    density += density_buffer[i] * diis_coefs[i];
}

void DIIS::solve() {
  fock = core_guess();

#ifndef NDEBUG
  std::ofstream log_stream;
  log_stream.open("logs/hcore.log", std::ios::app);
  log_stream << "\n\n\n\n" << std_m.H << std::endl;
  log_stream.close();
#endif

  for (auto iter = 1u; iter <= max_iter; ++iter) {
    if (fabs(cur_energy - prev_energy) < etol) {
      std::cout << "DIIS-SCF converged in " << iter << " iterations.\n";
      std::cout << "Total energy is: " << cur_energy + std_m.get_total_Vnn()
                << " Eh\n";
      return;
    }

    update_lcao_coefs(); // C first because D = f(C)

    if (iter > diis_size) { // DIIS
      update_extended_error_product();
      update_diis_coefs();
      density = expand(density_buffer, diis_coefs);
      fock = expand(fock_buffer, diis_coefs);
      error = expand(error_buffer, diis_coefs);
    } else { // before DIIS
      density = update_density();
      fock = update_fock();
      update_error(); // prev -> cur; cur = new
    }

    fock_buffer.update(fock);
    error_buffer.update(error);
    density_buffer.update(density);

    update_energy();
    print_iter(iter);
  }

  std::cout << "DIIS-SCF didn't converge in " << max_iter << " iterations.\n";
}
