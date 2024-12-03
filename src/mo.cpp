#include "mo.h"
#include <cmath>
#include <cstring>
#include <iostream>

MOs::MOs()
    : n(0)
{
}

MOs::~MOs()
{
    if (n > 0) {
        delete[] mo_energies;
        delete[] irrep;
    }
}

int MOs::init(const int n_new)
{
    if (n > 0) {
        throw std::runtime_error("Object Cl2 already have size!\n");
    }
    n = n_new;
    mo_energies = new double[n_new];
    irrep = new int[n_new];
    C.init(n_new);
    return 0;
}

int MOs::get_size()
{
    return n;
}

double MOs::get_mo_energy(const int i) const
{
    return mo_energies[i];
}

int MOs::get_irrep(const int& i) const
{
    return irrep[i];
}

int MOs::set_mo_energy(const int i, const double energy_i)
{
    mo_energies[i] = energy_i;
    return 0;
}

int MOs::set_total_energy(const double new_total_energy)
{
    total_energy = new_total_energy;
    return 0;
}

double MOs::get_total_energy() const
{
    return total_energy;
}

bool MOs::set_c2v(int* symmAO, const double& limit)
{
    if (n == 0) {
        throw std::runtime_error("Failed to set undifined MO symmetry!\n");
    }

    double tmpC;
    for (int idMO = 0; idMO < n; idMO++) {
        irrep[idMO] = 0;
        std::cout << "MO #" << idMO + 1 << '\n';
        for (int i = 0; i < n; i++) {
            tmpC = C.get_element(i, idMO);
            if (fabs(tmpC) > limit)
                std::cout << "    " << i + 1 << ' ' << tmpC << '\n';
        }
    }
    return true;
}

void MOs::set_mo_energies(const double* new_mo_energies)
{
    memcpy(mo_energies, new_mo_energies, n * sizeof(double));
}
