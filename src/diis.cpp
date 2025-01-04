// ----------------------------------------------------------- //
//             Suffix _diis_ means DIIS stage                  //
// members without _diis_ suffix realtes to the pre-DIIS stage //
// ----------------------------------------------------------- //

#include "diis.h"

#include <lapacke.h>

#include <assert.h>
#include <cmath>
#include <cstddef>
#include <deque>
#include <iostream>
#include <stdexcept>
#include <vector>

DIIS::DIIS(MO &mo, StandardMatrices &std_m)
    : SCF(mo, std_m), diis_size(DEFAULT_DIIS_SIZE) {
  error = matrix(diis_size);

  extended_diis_product = matrix(diis_size + 1);
  for (int i = 0; i < diis_size; ++i) {
    extended_diis_product[diis_size][i] = -1;
    extended_diis_product[i][diis_size] = -1;
  }
  extended_diis_product[diis_size][diis_size] = 0;

  diis_coefs = std::vector<double>(diis_size + 1);
}

void DIIS::update_error() {
  error = (fock * density * std_m.S) - (std_m.S * density * fock);
}

void DIIS::update_fock_buffer() {
  if (fock_buffer.empty()) {
    fock_buffer = std::deque<matrix>{fock};
    assert(fock_buffer.end()->data() != fock.data());
    return;
  }
  if (fock_buffer.size() == static_cast<std::size_t>(diis_size)) {
    fock_buffer.pop_front();
  }
  fock_buffer.push_back(fock);
  assert(fock_buffer.end()->data() != fock.data());
}

void DIIS::update_error_buffer() {
  if (error_buffer.empty()) {
    error_buffer = std::deque<matrix>{error};
    assert(error_buffer.end()->data() != error.data());
    return;
  }
  if (error_buffer.size() == static_cast<std::size_t>(diis_size)) {
    error_buffer.pop_front();
  }
  error_buffer.push_back(error);
  assert(error_buffer.end()->data() != error.data());
}

void DIIS::update_extended_error_product() {
  assert(fock_buffer.size() == static_cast<std::size_t>(diis_size));
  assert(error_buffer.size() == static_cast<std::size_t>(diis_size));

  for (int i = 0; i < diis_size; ++i) {
    for (int j = 0; j < diis_size; ++j) {
      extended_diis_product[i][j] =
          frobenius_product(error_buffer[i], error_buffer[j]);
    }
  }
}

// solving system: extended_diis_product * diis_coefs = [0, ..., 0, -1]
void DIIS::update_diis_coefs() {
  int N = diis_size + 1;
  int NRHS = 1;
  std::vector<int> IPIV(N);

  // Initialize RHS: [0, ..., 0, -1]
  for (int i = 0; i <= diis_size; ++i) {
    diis_coefs[i] = -(i == diis_size);
  }

  int INFO =
      LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, extended_diis_product.data(), N,
                    IPIV.data(), diis_coefs.data(), N);

  if (INFO != 0) {
    throw std::runtime_error("Failed to solve DIIS system for coefficients. "
                             "INFO = " +
                             std::to_string(INFO) +
                             ", Matrix Size = " + std::to_string(N));
  }
}

void DIIS::update_diis_error() {
  error.zeroize();
  for (int k = 0; k < diis_size; ++k) {
    error += error_buffer[k] * diis_coefs[k];
  }
}

void DIIS::update_diis_fock() {
  fock.zeroize();
  for (int k = 0; k < diis_size; ++k) {
    fock += fock_buffer[k] * diis_coefs[k];
  }
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
#endif

void DIIS::solve() {
  core_guess();

  for (int iter = 1; iter <= max_iter; ++iter) {
    if (fock_buffer.size() == static_cast<std::size_t>(diis_size)) {
      update_extended_error_product();
      update_diis_coefs();
      update_diis_fock();
      update_diis_error();
    }

    update_lcao_coefs();
    update_density();
    update_energy();
    print_iter(iter);

    if (fabs(cur_energy - prev_energy) < etol) {
      std::cout << "DIIS-SCF converged in " << iter << " iterations.\n";
      std::cout << "Total energy is: " << cur_energy + std_m.get_total_Vnn()
                << " Eh\n";
      return;
    }

    update_fock_buffer();
    update_error_buffer();
    update_fock();
    update_error();
  }

  std::cout << "DIIS-SCF didn't converge in " << max_iter << " iterations.\n";
}
