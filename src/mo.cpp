#include "mo.h"
#include <cmath>
#include <iostream>
#include <cstring>

MOs::MOs(): n(0) {}

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

double MOs::check_ij(const int i, const int j, const char* name) const
{
    if (n == 0) {
        std::cerr << "** ERROR *** array not defined in " << name << '\n';
        return 1;
    };
    if ((i < 0) || (i > n - 1)) {
        std::cerr << "** ERROR *** invalid number of string in " << name << '\n';
        return 1;
    };
    if ((j < 0) || (j > n - 1)) {
        std::cerr << "** ERROR *** invalid number of element in " << name << '\n';
        return 1;
    };
    return 0;
}

double MOs::get_mo_energy(const int i) const
{
    if (check_ij(i, 0, "get_mo_energy") > 0)
        return 1;
    return mo_energies[i];
}

int MOs::get_irrep(const int& i) const
{
    if (check_ij(i, 0, "get_irrep") > 0)
        return 1;
    return irrep[i];
}

int MOs::set_mo_energy(const int i, const double energy_i)
{
    if (check_ij(i, 0, "set_mo_energy") > 0)
        return 1;
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
        std::cerr << "Попытка определить симметрию МО, когда они еще не заданы!\n";
        return false;
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