#include "rhf.h"
#include "matrix.h"
#include "mo.h"
#include "standard_matrices.h"
#include <cmath>
#include <iomanip>
#include <iostream>

RHF::RHF(standard_matrices& std_m, MOs& mo)
    : etol(1e-12)
    , max_iter(100)
    , iter(0)
    , is_converged(false)
    , std_m(std_m)
    , mo(mo)
    , prev_energy(1)
    , cur_energy(0)
{
    evec = new double[std_m.get_nAO() * std_m.get_nAO()];
    mo_energies = new double[std_m.get_nAO()];

    density(std_m.get_nAO());
    eri_matrix(std_m.get_nAO());
    fock_matrix(std_m.get_nAO());

    std::cout << "\n-- Running SCF procedure --\n"
              << "Electrons: " << std_m.get_num_el() << std::endl
              << "Max iterations: " << max_iter << std::endl
              << "Energy tolerance: " << etol << std::endl;
}

RHF::~RHF()
{
    delete[] mo_energies;
    delete[] evec;
}

bool RHF::get_convergency()
{
    return is_converged;
}

matrix RHF::transform_matrix(const matrix& source_matrix)
{
    return std_m.X * source_matrix * std_m.X;
}

void RHF::validate_convergency()
{
    is_converged = (fabs(prev_energy - cur_energy) < etol);

    if (iter > max_iter) {
        throw std::runtime_error("SCF algorithm has not converged in " + std::to_string(max_iter) + " iterations!\n");
    }
}

void RHF::calculate_density()
{
    for (int i = 0; i < std_m.get_nAO(); ++i) {
        for (int j = 0; j < std_m.get_nAO(); ++j) {
            density[i][j] = 0;
            for (int k = 0; k < std_m.get_num_el() / 2; ++k) {
                density[i][j] += 2 * mo.C[i][k] * mo.C[j][k];
            }
        }
    }
}

void RHF::calculate_eri_matrix()
{
    for (int m = 0; m < std_m.get_nAO(); ++m) {
        for (int v = 0; v < std_m.get_nAO(); ++v) {
            eri_matrix[m][v] = 0;
            for (int l = 0; l < std_m.get_nAO(); ++l) {
                for (int s = 0; s < std_m.get_nAO(); ++s) {
                    eri_matrix[m][v] += density[l][s] * (std_m.get_Vee(m, v, s, l) - std_m.get_Vee(m, l, s, v) / 2);
                }
            }
        }
    }
}

void RHF::calculate_fock()
{
    fock_matrix = std_m.H + eri_matrix;
}

void RHF::update_energy()
{
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

void RHF::calculate_expansion()
{
    mo.C.from_array(evec);
    mo.C = std_m.X * mo.C.T();
}

void RHF::calculate_diis_coefs()
{
    // TODO: to be implemented
}

void RHF::calculate_diis_fock()
{
    fock_matrix.zeroize();
    for (int i = 0; i < diis_size; ++i) {
        fock_matrix += (fock_buffer[i] * diis_coefs[i]);
    }
}

void RHF::calculate_diis_error()
{
    error_matrix.zeroize();
    for (int i = 0; i < diis_size; ++i) {
        error_matrix += (error_buffer[i] * diis_coefs[i]);
    }
}

void RHF::update_error_buffer()
{
    if ((int)error_buffer.size() == diis_size) {
        error_buffer.pop_front();
    }
    error_buffer.push_back(error_matrix);
}

void RHF::update_fock_buffer()
{
    if ((int)fock_buffer.size() == diis_size) {
        fock_buffer.pop_front();
    }
    fock_buffer.push_back(fock_matrix);
}

void RHF::calculate_pre_diis_error()
{
    error_matrix = fock_matrix * density * std_m.S - std_m.S * density * fock_matrix;
}

void RHF::print_iteration()
{
    std::cout << "# "
              << std::setw(5) << iter << std::setw(20)
              << std::setprecision(12) << cur_energy + std_m.get_total_Vnn() << '\n';
}

void RHF::core_guess()
{
    mo.init(std_m.get_nAO());
    direct_iteration();
}

void RHF::direct_iteration()
{
    transform_matrix(fock_matrix).eigen_vv(evec, mo_energies);
    mo.set_mo_energies(mo_energies);
    calculate_expansion();
}

void RHF::solve_rhf()
{
    core_guess();
    std::cout << "--Roothan-Hall algorithm started --\n";
    while (!is_converged) {
        if (iter == diis_size) {
            std::cout << "-- DIIS approximation started --\n";
        }
        iter < diis_size ? roothan_hall_step() : diis_step();
        update_fock_buffer();
        update_error_buffer();
        direct_iteration();
        iter++;
        print_iteration();
        validate_convergency();
        update_energy();
    }
}

void RHF::roothan_hall_step()
{
    calculate_density();
    calculate_eri_matrix();
    calculate_fock();
    calculate_pre_diis_error();
}
void RHF::diis_step()
{
    calculate_diis_coefs();
    calculate_diis_fock();
    calculate_diis_error();
}

void RHF::roothan_hall()
{
    std::cout << "\n-- Roothan-Hall algorithm started --\n\n";

    core_guess();

    while (!is_converged) {

        calculate_density();

        calculate_eri_matrix();

        calculate_fock();

        direct_iteration();

        iter++;

        print_iteration();

        validate_convergency();

        update_energy();
    }

    std::cout << "\nTotal iterations = " << iter << '\n';
    mo.set_total_energy(cur_energy + std_m.get_total_Vnn());
}
