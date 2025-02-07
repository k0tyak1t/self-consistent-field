#include "scf.h"

#include <assert.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

SCF::SCF(MO &mo, StandardMatrices &std_m)
    : etol(DEFAULT_ETOL), max_iter(DEFAULT_MAX_ITER), prev_energy(0.0),
      cur_energy(1.0), mo(mo), std_m(std_m), nAO(std_m.get_nAO()) {
  mo_energies = std::vector<double>(nAO);
  density = matrix(nAO);
  fock = matrix(nAO);
  lcao_coefs = matrix(nAO);
}

double SCF::get_total_energy() { return cur_energy + std_m.get_total_Vnn(); }

// A' = XAX^{\dagger}
matrix SCF::transform_matrix(const matrix &mat) const {
  return std_m.X * mat * std_m.X;
}

void SCF::print_iter(int iter) const {
  std::cout << "#" << std::setw(4) << iter << std::setw(18)
            << std::setprecision(12) << cur_energy + std_m.get_total_Vnn()
            << '\n';
}

// F = H^{core}
void SCF::core_guess() { fock = std_m.H; }

// F'C' = C'E => C' => C
void SCF::update_lcao_coefs() {
  transform_matrix(fock).eigen_vv(lcao_coefs.data(), mo_energies.data());
  lcao_coefs = std_m.X * lcao_coefs.transposed();

#ifndef NDEBUG
  std::cout << "mo energies: ";
  for (auto i : mo_energies) {
    std::cout << " " << i;
  }
  std::cout << std::endl;
#endif
}

// D = CC^{\dagger}
// D_{\mu\nu} = 2 \sum_{a = 1}^{n_{el}/2}{C_{\mu a}C_{\nu a}}
// src: "Modern Quantum Chemistry", p.139, f-No. 3.145
void SCF::update_density() {
  for (std::size_t m = 0; m < nAO; ++m) {
    for (std::size_t n = 0; n < nAO; ++n) {
      density[m][n] = 0;
      for (std::size_t a = 0; a < std_m.get_num_el() / 2; ++a) {
        density[m][n] += 2 * lcao_coefs[m][a] * lcao_coefs[n][a];
      }
    }
  }
}

// F_{\mu\nu} = H^{core}_{\mu\nu} +
// \sum_{\lambda\sigma}{D_{\lambda\sigma} * \left[(\mu\nu|\sigma\lambda) -
// \frac{1}{2}(\mu\lambda|\sigma\nu)\right]}
// src: "Modern Quantum Chemistry", p. 141, f-No. 3.154
void SCF::update_fock() {
  for (std::size_t m = 0; m < nAO; ++m) {
    for (std::size_t n = 0; n < nAO; ++n) {
      fock[m][n] = std_m.H[m][n];
      for (std::size_t l = 0; l < nAO; ++l) {
        for (std::size_t s = 0; s < nAO; ++s) {
          fock[m][n] += density[l][s] * (std_m.get_eri(m, n, s, l) -
                                         0.5 * std_m.get_eri(m, l, s, n));
        }
      }
    }
  }
#ifndef NDEBUG
  std::cout << "Fock matrix: \n" << fock;
#endif
}

#if 0
// F_{\mu\nu} = H^{core}_{\mu\nu}
// + \sum_{a=1}^{n_{el}/2}{2(\mu\nu|aa) - (\mu a|a\nu)}
// src: "Modern Quantum Chemistry", p. 140, f-No. 3.148
// Fock matrix for restricted Hartree-Fock method
void SCF::update_fock() {
  for (int m = 0; m < nAO; ++m)
    for (int n = 0; n < nAO; ++n) {
      fock[m][n] = std_m.H[m][n];
      for (int a = 0; a < std_m.get_num_el() / 2; ++a) {
        fock[m][n] += 2 * std_m.get_eri(m, n, a, a) - std_m.get_eri(m, a, a, n);
      }
    }
}
#endif

void SCF::update_energy() {
  prev_energy = cur_energy;

  cur_energy = 0;
  for (std::size_t m = 0; m < nAO; ++m) {
    for (std::size_t n = 0; n < nAO; ++n) {
      cur_energy += density[n][m] * (std_m.H[m][n] + fock[m][n]) / 2;
    }
  }
}

void SCF::solve() {
  core_guess();
#ifndef NDEBUG
  std::cout << "H core: \n" << std_m.H << std::endl;
#endif

  for (int iter = 1; iter <= max_iter; ++iter) {
    update_lcao_coefs();
    update_density();
    update_fock();
    update_energy();
    print_iter(iter);

    if (fabs(cur_energy - prev_energy) < etol) {
      std::cout << "SCF converged in " << iter << " iterations.\n";
      std::cout << "Total energy is: " << cur_energy + std_m.get_total_Vnn()
                << " Eh\n";
      return;
    }
  }

  std::cerr << "SCF did not converge in " << max_iter << " iterations.\n";
}
