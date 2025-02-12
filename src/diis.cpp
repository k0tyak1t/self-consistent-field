// ----------------------------------------------------------- //
//             Suffix _diis_ means DIIS stage                  //
// members without _diis_ suffix realtes to the pre-DIIS stage //
// ----------------------------------------------------------- //

#include "diis.h"
#include "matrix.h"

#include <algorithm>
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
      extended_diis_product(Matrix{diis_size + 1}) {
  error = Matrix(diis_size);

  for (std::size_t i = 0; i < diis_size; ++i) {
    extended_diis_product[diis_size][i] = 1.0;
    extended_diis_product[i][diis_size] = 1.0;
  }

  diis_coefs = std::vector<double>(diis_size + 1);
}

void DIIS::update_error() {
  Matrix cinv = lcao_coefs.inverse();
  Matrix fock_mo = cinv * fock * cinv.transposed();
#if 1
  error = (fock * density * std_m.S) - (std_m.S * density * fock);
#else // ~0
  std::size_t n_occ = std_m.get_num_el() / 2;
  std::size_t n_virt = std_m.get_nAO() - n_occ;

  error = Matrix{std::max(n_occ, n_virt)};

  for (auto i = 0u; i < n_occ; ++i)
    for (auto j = 0u; j < n_virt; ++j)
      error[i][j] = fock_mo[i][j + n_occ];

#ifndef NDEBUG
  std::cout << "num occ.: " << n_occ << " num virt.: " << n_virt << std::endl;
#endif // NDEBUG
#endif // ~0
}

void DIIS::update_fock_buffer() {
#ifndef NDEBUG
  std::cout << "[DIIS]: fock buffer update...\n";
#endif
  if (fock_buffer.size() == diis_size) {
    fock_buffer.pop_front();
  }
  fock_buffer.push_back(fock);

#ifndef NDEBUG
  std::cout << "[DIIS]: fock buffer updated.\n";
#endif
}

void DIIS::update_error_buffer() {
#ifndef NDEBUG
  std::cout << "[DIIS]: error buffer update...\n";
#endif
  if (error_buffer.size() == diis_size) {
    error_buffer.pop_front();
  }
  error_buffer.push_back(error);
#ifndef NDEBUG
  std::cout << "[DIIS]: error buffer updated.\n";
#endif
}

void DIIS::update_extended_error_product() {
  assert(fock_buffer.size() == diis_size);
  assert(error_buffer.size() == diis_size);

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
  int INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, NRHS, copy.data(), N,
                           IPIV.data(), diis_coefs.data(), N);

  if (INFO != 0)
    throw std::runtime_error("Failed to solve DIIS system for coefficients. "
                             "INFO = " +
                             std::to_string(INFO) +
                             ", matrix size = " + std::to_string(N));

#if 1 // normalization
  double sum = 0.0;
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
#ifndef NDIIS
  std::cout << "[DIIS]: DIIS stage error matrix update...\n";
#endif
  error = Matrix::zero_like(error);
  for (std::size_t k = 0; k < diis_size; ++k)
    error += error_buffer[k] * diis_coefs[k];

#ifndef NDIIS
  std::cout << "[DIIS]: DIIS stage error matrix updated.\n";
#endif
}

void DIIS::update_diis_fock() {
#ifndef NDIIS
  std::cout << "[DIIS]: DIIS stage fock matrix update...\n";
#endif
  fock = Matrix::zero_like(fock);
  for (std::size_t k = 0; k < diis_size; ++k)
    fock += fock_buffer[k] * diis_coefs[k];

#ifndef NDIIS
  std::cout << "[DIIS]: DIIS stage fock matrix updated.\n";
#endif
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

    if (fock_buffer.size() == diis_size) {
      update_extended_error_product();
      update_diis_coefs();
      update_diis_fock();
      update_diis_error();
    }

    if (fock_buffer.size() < diis_size) {
      update_fock();
      update_error();
    }
    update_fock_buffer();
    update_error_buffer();
    update_energy();
    print_iter(iter);
  }

  std::cout << "DIIS-SCF didn't converge in " << max_iter << " iterations.\n";
}
