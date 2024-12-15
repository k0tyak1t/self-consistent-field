#include "standard_matrices.h"
#include "matrix.h"
#include <iostream>

standard_matrices::~standard_matrices()
{
    if (nAO > 0) {
        delete[] Vee;
    }
}

standard_matrices::standard_matrices(const unsigned int nAO): nAO(nAO)
{
    if (!nAO) {return; }

    S = matrix(nAO);
    T = matrix(nAO);
    H = matrix(nAO);
    Ven = matrix(nAO);
    Vee = new double[nAO * nAO * nAO * nAO];
}

standard_matrices::standard_matrices(): standard_matrices(0) {};

// standard_matrices::standard_matrices(int new_nAO)
// {
//     if (new_nAO <= 0) {
//         std::cerr << "size of nAO <=0 \n";
//     }
//     nAO = new_nAO;
//     S.init(new_nAO);
//     T.init(new_nAO);
//     H.init(new_nAO);
//     Ven.init(new_nAO);
//     Vee = new double[new_nAO * new_nAO * new_nAO * new_nAO];
// }

void standard_matrices::init(const unsigned int n_new)
{
    if (nAO) {
        throw std::runtime_error("AO already exists!");
    }

    nAO = n_new;
    S.init(n_new);
    T.init(n_new);
    H.init(n_new);
    Ven.init(n_new);
    Vee = new double[n_new * n_new * n_new * n_new];
}

void standard_matrices::set_total_Vnn(const double new_total_Vnn) {total_Vnn = new_total_Vnn;}

double standard_matrices::get_total_Vnn() const {return total_Vnn;}

const int standard_matrices::get_nAO() const {return nAO;}

const double standard_matrices::get_Vee(int i, int j, int k, int l) const
{ // (ij|kl)
    if (((i < 0) || (i >= nAO)) || ((j < 0) || (j >= nAO)) || ((k < 0) || (k >= nAO)) || ((l < 0) || (l >= nAO))) {
        throw std::runtime_error("Failed to access element!");
    }

    return Vee[i + nAO * j + nAO * nAO * k + nAO * nAO * nAO * l];
}

void standard_matrices::set_Vee(int i, int j, int k, int l, double Vijkl)
{ // (ij|kl)
    if (((i < 0) || (i >= nAO)) || ((j < 0) || (j >= nAO)) || ((k < 0) || (k >= nAO)) || ((l < 0) || (l >= nAO))) {
        throw std::runtime_error("Failed to access element!");
    }
    Vee[i + nAO * j + nAO * nAO * k + nAO * nAO * nAO * l] = Vijkl;
}

const int standard_matrices::get_num_el() const {return num_el;}

void standard_matrices::set_num_el(const unsigned int new_num_el)
{
    if (new_num_el % 2) {
        throw std::runtime_error("Failed to solve RHF with odd number of electrons!");
    }

    num_el = new_num_el;
}
