#include "mo.h"

#include <cmath>
#include <cstring>
#include <iostream>

MO::MO(const unsigned n) : n(n), mo_energies(n), irrep(n), C(n) {}

int MO::get_size() const { return n; }

double MO::get_mo_energy(const int i) const { return mo_energies[i]; }

int MO::get_irrep_characters(const int &i) const { return irrep[i]; }

int MO::set_mo_energy(const int i, const double energy_i) {
  mo_energies[i] = energy_i;
  return 0;
}

int MO::set_total_energy(const double new_total_energy) {
  total_energy = new_total_energy;
  return 0;
}

double MO::get_total_energy() const { return total_energy; }

bool MO::set_c2v(int *symmAO, const double &limit) {
  if (n == 0)
    throw std::runtime_error("Failed to set undifined MO symmetry!");

  double tmpC{};
  for (auto idMO = 0; idMO < n; idMO++) {
    irrep[idMO] = 0;
    std::cout << "MO #" << idMO + 1 << std::endl;
    for (int i = 0; i < n; i++) {
      tmpC = C(i, idMO);
      if (fabs(tmpC) > limit)
        std::cout << "    " << i + 1 << ' ' << tmpC << std::endl;
    }
  }
  return true;
}

void MO::set_mo_energies(const double *new_mo_energies) {
  std::copy(new_mo_energies, new_mo_energies + n, mo_energies.begin());
}
