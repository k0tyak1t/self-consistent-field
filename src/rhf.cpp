#include "rhf.h"

#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>

#include "matrix.h"
#include "mo.h"
#include "standard_matrices.h"

RHF::RHF(StandardMatrices &std_m, MO &mo)
    : etol(1e-12), max_iter(100), iter(0), is_converged(false), std_m(std_m),
      mo(mo), prev_energy(1), cur_energy(0) {
  evec = new double[std_m.get_nAO() * std_m.get_nAO()];
  mo_energies = new double[std_m.get_nAO()];

  density(std_m.get_nAO());
  eri_matrix(std_m.get_nAO());
  fock_matrix(std_m.get_nAO());
  error_product_matrix(diis_size + 1);
  init_error_product_matrix();

  std::cout << "\n-- Running SCF procedure --\n"
            << "Electrons: " << std_m.get_num_el() << std::endl
            << "Max iterations: " << max_iter << std::endl
            << "Energy tolerance: " << etol << "\n\n";
}

RHF::~RHF() {
  delete[] mo_energies;
  delete[] evec;
}

void RHF::init_error_product_matrix() {
  for (int i = 0; i < diis_size + 1; ++i) {
    error_product_matrix[i][diis_size] = -1;
    error_product_matrix[diis_size][i] = -1;
  }

  error_product_matrix[diis_size][diis_size] = 0;
}

bool RHF::get_convergency() { return is_converged; }

matrix RHF::transform_matrix(const matrix &source_matrix) {
  return std_m.X * source_matrix * std_m.X;
}

void RHF::validate_convergency() {
  is_converged = (fabs(prev_energy - cur_energy) < etol);

  if (iter > max_iter) {
    throw std::runtime_error("SCF algorithm has not converged in " +
                             std::to_string(max_iter) + " iterations!\n");
  }
}

void RHF::calculate_density() {
  for (int i = 0; i < std_m.get_nAO(); ++i) {
    for (int j = 0; j < std_m.get_nAO(); ++j) {
      density[i][j] = 0;
      for (int k = 0; k < std_m.get_num_el() / 2; ++k) {
        density[i][j] += 2 * mo.C[i][k] * mo.C[j][k];
      }
    }
  }
}

void RHF::calculate_eri_matrix() {
  for (int m = 0; m < std_m.get_nAO(); ++m) {
    for (int v = 0; v < std_m.get_nAO(); ++v) {
      eri_matrix[m][v] = 0;
      for (int l = 0; l < std_m.get_nAO(); ++l) {
        for (int s = 0; s < std_m.get_nAO(); ++s) {
          eri_matrix[m][v] += density[l][s] * (std_m.get_eri(m, v, s, l) -
                                               std_m.get_eri(m, l, s, v) / 2);
        }
      }
    }
  }
}

void RHF::calculate_fock() { fock_matrix = std_m.H + eri_matrix; }

void RHF::update_energy() {
  prev_energy = cur_energy;
  cur_energy = 0;

  for (int i = 0; i < std_m.get_num_el() / 2; ++i) {
    for (int j = 0; j < std_m.get_nAO(); ++j) {
      for (int k = 0; k < std_m.get_nAO(); ++k) {
        cur_energy += mo.C[j][i] * mo.C[k][i] * std_m.H[j][k];
      }
    }
  }

  for (int i = 0; i < std_m.get_num_el() / 2; ++i) {
    cur_energy += mo_energies[i];
  }
}

void RHF::calculate_expansion() {
  mo.C.from_array(evec);
  mo.C = std_m.X * mo.C.T();
}

void RHF::calculate_error_product_matrix() {
  for (int i = 0; i < diis_size; ++i) {
    for (int j = 0; j < diis_size; ++j) {
      error_product_matrix[i][j] =
          frobenius_product(error_buffer[i], error_buffer[j]);
    }
  }
}

void RHF::calculate_diis_coefs() {
  std::cout << error_product_matrix;
  std::cout.flush();
  std::vector<double> b(diis_size, 0.0);
  b.push_back(-1.0);
#if 0
diis_coefs = error_product_matrix.inv() * b;
#endif
}

void RHF::calculate_diis_fock() {
  fock_matrix.zeroize();
  for (int i = 0; i < diis_size; ++i) {
    fock_matrix += (fock_buffer[i] * diis_coefs[i]);
  }
}

void RHF::calculate_diis_error() {
  error_matrix.zeroize();
  for (int i = 0; i < diis_size; ++i) {
    error_matrix += (error_buffer[i] * diis_coefs[i]);
  }
}

void RHF::verify_buffer(std::deque<matrix> &buffer) {
  if (!buffer.empty() && buffer.size() == (std::size_t)diis_size) {
    buffer.pop_front();
  }
}

void RHF::update_buffer(std::deque<matrix> &buffer, const matrix &new_matrix) {
  verify_buffer(buffer);
  buffer.push_back(new_matrix);
}

void RHF::calculate_rh_error() {
  error_matrix =
      fock_matrix * density * std_m.S - std_m.S * density * fock_matrix;
}

void RHF::print_iteration() {
  std::cout << "#" << std::setw(4) << iter << std::setw(18)
            << std::setprecision(12) << cur_energy + std_m.get_total_Vnn()
            << '\n';
}

void RHF::core_guess() {
  // mo.init(std_m.get_nAO());
  // mo(std_m.get_nAO());
  common_step();
}

void RHF::common_step() {
  transform_matrix(fock_matrix).eigen_vv(evec, mo_energies);
  mo.set_mo_energies(mo_energies);
  calculate_expansion();
  print_iteration();
  validate_convergency();
  update_energy();
  iter++;
}

void RHF::solve_rhf() {
  core_guess();
  std::cout << "--Roothan-Hall algorithm started --\n";
  while (!is_converged) {
    roothan_hall_step(); // should be deleted when diis implemented

    // if (iter == diis_size) {
    //     std::cout << "-- DIIS approximation started --\n";
    // }
    // iter < diis_size ? roothan_hall_step() : diis_step();
    update_buffer(fock_buffer, fock_matrix);
    update_buffer(error_buffer, error_matrix);
    common_step();
  }

  std::cout << "\nTotal iterations = " << iter << '\n';
  mo.set_total_energy(cur_energy + std_m.get_total_Vnn());
}

void RHF::roothan_hall_step() {
  calculate_density();
  calculate_eri_matrix();
  calculate_fock();
  calculate_rh_error();
}

void RHF::diis_step() {
  calculate_error_product_matrix();
  calculate_diis_coefs();
  calculate_diis_fock();
  calculate_diis_error();
}
