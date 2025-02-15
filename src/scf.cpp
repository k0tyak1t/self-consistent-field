#include "scf.h"
#include <assert.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

SCF::SCF(MO &mo, StandardMatrices &std_m)
    : etol(DEFAULT_ETOL), max_iter(DEFAULT_MAX_ITER), prev_energy(0.0),
      cur_energy(1.0), std_m(std_m), nAO(std_m.get_nAO()), mo(mo),
      mo_energies(nAO), density(nAO), fock(nAO), lcao_coefs(nAO) {}

double SCF::get_total_energy() { return cur_energy + std_m.get_total_Vnn(); }

// A' = XAX^{\dagger}
Matrix SCF::transform_matrix(const Matrix &mat) const {
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
#ifndef NDEBUG
  std::cout << "[SCF]: lcao coefs update...\n";
#endif
  transform_matrix(fock).eigen_vv(lcao_coefs.data(), mo_energies.data());
  lcao_coefs = std_m.X * Matrix::transposed(lcao_coefs);

#ifndef NDEBUG
  std::ofstream log_stream;
  log_stream.open("logs/mo.log", std::ios::app);
  log_stream << "mo energies: ";
  for (auto i : mo_energies) {
    log_stream << " " << i;
  }
  log_stream << std::endl;
  std::cout << "[SCF]: lcao coefs updated.\n";
#endif
}

// D = CC^{\dagger}
// D_{\mu\nu} = 2 \sum_{a = 1}^{n_{el}/2}{C_{\mu a}C_{\nu a}}
// src: "Modern Quantum Chemistry", p.139, f-No. 3.145
void SCF::update_density() {
#ifndef NDEBUG
  std::cout << "[SCF]: density matrix update...\n";
#endif
  density = Matrix::zero_like(density);
  for (std::size_t m = 0; m < nAO; ++m)
    for (std::size_t n = 0; n < nAO; ++n)
      for (std::size_t a = 0; a < std_m.get_num_el() / 2; ++a)
        density[m][n] += 2 * lcao_coefs[m][a] * lcao_coefs[n][a];
#ifndef NDEBUG
  std::cout << "[SCF]: density matrix updated.\n";
#endif
}

// F_{\mu\nu} = H^{core}_{\mu\nu} +
// \sum_{\lambda\sigma}{D_{\lambda\sigma} * \left[(\mu\nu|\sigma\lambda) -
// \frac{1}{2}(\mu\lambda|\sigma\nu)\right]}
// src: "Modern Quantum Chemistry", p. 141, f-No. 3.154
void SCF::update_fock() {
#ifndef NDEBUG
  std::cout << "[SCF]: fock matrix update...\n";
#endif
  for (std::size_t m = 0; m < nAO; ++m)
    for (std::size_t n = 0; n < nAO; ++n) {
      fock[m][n] = std_m.H[m][n];
      for (std::size_t l = 0; l < nAO; ++l)
        for (std::size_t s = 0; s < nAO; ++s)
          fock[m][n] += density[l][s] * (std_m.get_eri(m, n, s, l) -
                                         0.5 * std_m.get_eri(m, l, s, n));
    }
#ifndef NDEBUG
  std::ofstream log_stream;
  log_stream.open("logs/fock.log", std::ios::app);
  log_stream << "\n\n" << fock;
  log_stream.close();
  std::cout << "[SCF]: fock matrix updated.\n";
#endif
}

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
  std::ofstream log_stream;
  log_stream.open("logs/hcore.log", std::ios::app);
  log_stream << "\n\n" << std_m.H << std::endl;
  log_stream.close();
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

  std::cerr << "SCF didn't converge in " << max_iter << " iterations.\n";
}
