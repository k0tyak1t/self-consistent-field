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
#include <deque>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

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
    extended_diis_product[diis_size][i] = 1.0;
    extended_diis_product[i][diis_size] = 1.0;
  }
}

void DIIS::update_error() {
#if 1
  error = fock * density * std_m.S - std_m.S * density * fock;
#else
  Matrix cinv = lcao_coefs.inverse();
  Matrix fock_mo = cinv * fock * Matrix::transposed(cinv);
  std::size_t n_occ = std_m.get_num_el() / 2;
  std::size_t n_virt = std_m.get_nAO() - n_occ;

  error = Matrix{n_occ > n_virt ? n_occ : n_virt};

  for (auto i = 0u; i < n_occ; ++i)
    for (auto j = 0u; j < n_virt; ++j)
      error[i][j] = fock_mo[i][j + n_occ];

#ifndef NDEBUG
  std::cout << "num occ.: " << n_occ << " num virt.: " << n_virt << std::endl;
#endif // NDEBUG
#endif // ~0
}

void DIIS::update_extended_error_product() {
  for (std::size_t i = 0; i < diis_size; ++i)
    for (std::size_t j = 0; j < diis_size; ++j)
      extended_diis_product[i][j] =
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

  // Initialize RHS: [0, ..., 0, 1]
  diis_coefs = std::vector<double>(diis_size + 1);
  diis_coefs[diis_size] = 1.0;

  Matrix copy{extended_diis_product};
  int INFO = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, copy.data(), N,
                           IPIV.data(), diis_coefs.data(), N);

  if (INFO != 0)
    throw std::runtime_error("Failed to solve DIIS system for coefficients. "
                             "INFO = " +
                             std::to_string(INFO) +
                             ", matrix size = " + std::to_string(N));

#ifndef NNORMALIZE // normalization
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
  core_guess();
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

    update_lcao_coefs();
    update_density();

    if (iter > diis_size) {
      update_extended_error_product();
      update_diis_coefs();
      update_diis_fock();
      update_diis_error();
      update_diis_density();
    }

    else {
      update_fock();
      update_error();
    }

    fock_buffer.update(fock);
    error_buffer.update(error);
    density_buffer.update(density);

    update_energy();
    print_iter(iter);
  }

  std::cout << "DIIS-SCF didn't converge in " << max_iter << " iterations.\n";
}
