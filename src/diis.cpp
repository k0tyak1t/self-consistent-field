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
      extended_diis_product(Matrix{diis_size + 1}) {
  error = Matrix(diis_size);

  for (std::size_t i = 0; i < diis_size; ++i) {
    extended_diis_product[diis_size][i] = -1.0;
    extended_diis_product[i][diis_size] = -1.0;
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

  error = Matrix(std::max(n_occ, n_virt));

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
  if (fock_buffer.empty()) {
    fock_buffer = std::deque<Matrix>{fock};
    assert(fock_buffer.end()->data() != fock.data());
    return;
  }
  if (fock_buffer.size() == static_cast<std::size_t>(diis_size)) {
    fock_buffer.pop_front();
  }
  fock_buffer.push_back(fock);
  assert(fock_buffer.end()->data() != fock.data());

#ifndef NDEBUG
  std::cout << "[DIIS]: fock buffer updated.\n";
#endif
}

void DIIS::update_error_buffer() {
#ifndef NDEBUG
  std::cout << "[DIIS]: error buffer update...\n";
#endif
  if (error_buffer.empty()) {
    error_buffer = std::deque<Matrix>{error};
    assert(error_buffer.end()->data() != error.data());
    return;
  }
  if (error_buffer.size() == static_cast<std::size_t>(diis_size)) {
    error_buffer.pop_front();
  }
  error_buffer.push_back(error);
  assert(error_buffer.end()->data() != error.data());
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
          frobenius_product(error_buffer[i], error_buffer[j]);
}

// solving system: extended_diis_product * diis_coefs = [0, ..., 0, -1]
void DIIS::update_diis_coefs() {
  int N = diis_size + 1;
  int NRHS = 1;
  std::vector<int> IPIV(N);

  // Initialize RHS: [0, ..., 0, -1]
  for (std::size_t i = 0; i <= diis_size; ++i) {
    diis_coefs[i] = (diis_size == i) ? -1.0 : 0.0;
  }

  Matrix copy = extended_diis_product;
  int INFO = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, copy.data(), N,
                           IPIV.data(), diis_coefs.data(), N);

  if (INFO != 0) {
    throw std::runtime_error("Failed to solve DIIS system for coefficients. "
                             "INFO = " +
                             std::to_string(INFO) +
                             ", Matrix Size = " + std::to_string(N));
  }

#if 1
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
  error.zeroize();
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
  fock.zeroize();
  for (std::size_t k = 0; k < diis_size; ++k)
    fock += fock_buffer[k] * diis_coefs[k];

#ifndef NDIIS
  std::cout << "[DIIS]: DIIS stage fock matrix updated.\n";
#endif
}

#if 0
void DIIS::solve() {
  core_guess();

  for (int iter = 1; iter <= max_iter; ++iter) {
    update_lcao_coefs();
    update_density();
    if (iter <= diis_size) {
      update_fock();
      update_error();
    } else {
      update_extended_error_product();
      update_diis_coefs();
      update_diis_fock();
      update_diis_error();
    }
    update_fock_buffer();
    update_error_buffer();
    update_energy();
    print_iter(iter);

    if (fabs(cur_energy - prev_energy) < etol) {
      std::cout << "DIIS-SCF converged in " << iter << " iterations.\n";
      std::cout << "Total energy is: " << cur_energy + std_m.get_total_Vnn()
                << " Eh\n";
      return;
    }
  }

  std::cout << "DIIS-SCF didn't converge in " << max_iter << " iterations.\n";
}
#else
void DIIS::solve() {
  core_guess();
#ifndef NDEBUG
  std::ofstream log_stream;
  log_stream.open("logs/hcore.log", std::ios::app);
  log_stream << "\n\n" << std_m.H << std::endl;
  log_stream.close();
#endif
#endif

  for (int iter = 1; iter <= max_iter; ++iter) {
    if (fabs(cur_energy - prev_energy) < etol) {
      std::cout << "DIIS-SCF converged in " << iter << " iterations.\n";
      std::cout << "Total energy is: " << cur_energy + std_m.get_total_Vnn()
                << " Eh\n";
      return;
    }

    update_lcao_coefs();
    update_density();

    if (fock_buffer.size() == static_cast<std::size_t>(diis_size)) {
      update_extended_error_product();
      update_diis_coefs();
      update_diis_fock();
      update_diis_error();
    }

    if (fock_buffer.size() < static_cast<std::size_t>(diis_size)) {
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
