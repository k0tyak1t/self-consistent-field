#ifndef RHF_H
#define RHF_H
#include "matrix.h"
#include "mo.h"
#include "standard_matrices.h"
#include <deque>
#include <vector>

class RHF {
public:
    RHF(standard_matrices&, MOs&);
    ~RHF();
    void init_error_product_matrix();
    bool get_convergency();
    matrix transform_matrix(const matrix&);
    void calculate_error_product_matrix();
    void validate_convergency();
    void print_iteration();
    void calculate_fock();
    void calculate_rh_error();
    void calculate_diis_coefs();
    void calculate_diis_fock();
    void calculate_diis_error();
    void update_buffer(std::deque<matrix>&, matrix);
    void calculate_density();
    void calculate_eri_matrix();
    void update_energy();
    void calculate_expansion();
    void core_guess();
    void common_step();
    void roothan_hall_step();
    void diis_step();
    void solve_rhf();

private:
    double etol;
    int max_iter, diis_size, iter;
    bool is_converged;
    standard_matrices& std_m;
    MOs& mo;
    double prev_energy, cur_energy;
    double* evec;
    double* mo_energies;
    matrix density, eri_matrix, fock_matrix, error_matrix, error_product_matrix;
    std::deque<matrix> fock_buffer, error_buffer;
    std::vector<double> diis_coefs;
};

#endif
