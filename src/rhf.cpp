#include <iomanip>
#include <iostream>
#include <cmath>
#include "rhf.h"
#include "matrix.h"
#include "standard_matrices.h"
#include "mo.h"

RHF::RHF(standard_matrices& std_m, MOs& mo)
  : etol(1e-12), max_iter(100), diis_size(5), is_converged(false), diis_mode(false),
    std_m(std_m), mo(mo), E_old(0), E_new(0)
{
  evec = new double[std_m.get_nAO() * std_m.get_nAO()];
  eval = new double[std_m.get_nAO()];

  density(std_m.get_nAO());
  eri_matrix(std_m.get_nAO());
  fock_matrix(std_m.get_nAO());

	std::cout << "\n-- Running SCF procedure --\n"
	          << "Electrons: " << std_m.get_num_el() << std::endl
	          << "Max iterations: " << max_iter << std::endl
	          << "Energy tolerance: " << etol << std::endl;
}

RHF::~RHF(){
  delete[] eval;
  delete[] evec;
}

bool RHF::get_convergency() {
	return is_converged;
}

void RHF::validate_diis_condition(const int rh_iter)
{
  diis_mode = ((int)error_buffer.size() >= diis_size);
}

void RHF::validate_convergency(const int iter) {
	is_converged = (fabs(E_old - E_new) < etol);

  if (iter > max_iter) {
    throw std::runtime_error("SCF algorithm has not converged in " + std::to_string(max_iter) + "iterations!\n");
  }
}

void RHF::build_density()
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
  					eri_matrix[m][v] += density[l][s] * (std_m.get_Vee(m,v,s,l) - std_m.get_Vee(m,l,s,v) / 2);
  				}
  			}
		}
	}
}

void RHF::calculate_fock_transformed()
{
  fock_matrix= std_m.H + eri_matrix;
  fock_matrix= std_m.X * fock_matrix * std_m.X;
}

double RHF::recalculate_energy()
{
  double E_new = 0;

	for (int i = 0; i < std_m.get_num_el() / 2; ++i) {
		for (int j = 0; j < std_m.get_nAO(); ++j) {
			for (int k = 0; k < std_m.get_nAO(); ++k) {
				E_new += mo.C[j][i] * mo.C[k][i] * std_m.H[j][k];
			}
		}
	}

	for (int i = 0; i < std_m.get_num_el() / 2; ++i) {
		E_new += eval[i];
	}

  return E_new;
}

void RHF::calculate_coef_matrix()
{
  mo.C.from_array(evec);
  mo.C = std_m.X * mo.C.T();
}

void RHF::verbose_iteration(const int iter)
{
  std::cout << "#"
            << std::setw(5) << iter << std::setw(20)
            << std::setprecision(12) << E_new + std_m.get_total_Vnn()<<'\n';
}

void RHF::core_guess()
{
	mo.init(std_m.get_nAO());
  fock_matrix= std_m.X * fock_matrix * std_m.X;
	fock_matrix.eigen_vv(evec, eval);
	mo.C.from_array(evec);
	mo.C = std_m.X * mo.C.T();
	
	for (int i = 0; i < std_m.get_nAO(); ++i) {
		mo.set_mo_energy(i, eval[i]);
	}
}

void RHF::roothan_hall(){
  std::cout << "\n** Roothan-Hall algorithm started **\n\n";
	int iter = 0;

  core_guess();

	while(!is_converged && !diis_mode) {

    build_density();

    calculate_eri_matrix();

    calculate_fock_transformed();

		fock_matrix.eigen_vv(evec, eval);

    calculate_coef_matrix();

    calculate_error();

		for (int i = 0; i < std_m.get_nAO(); ++i) {
			mo.set_mo_energy(i, eval[i]);
    }

		E_old = E_new;
    E_new = recalculate_energy();

		iter++;
    verbose_iteration(iter);
		validate_convergency(iter);
	}

	std::cout << "\nTotal iterations = " << iter << '\n';
	mo.set_total_energy(E_new + std_m.get_total_Vnn());
}

void RHF::calculate_error()
{
  matrix error_matrix(std_m.get_nAO());
  error_matrix = fock_matrix * density * std_m.S - std_m.S * density * fock_matrix;
  error_buffer.push_back(error_matrix);
}

void RHF::diis()
{
  std::cout << "** DIIS approximation started **\n";
  int iter = 0;
}
