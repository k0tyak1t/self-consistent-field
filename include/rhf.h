#ifndef RHF_H
#define RHF_H
#include <vector>
#include <deque>
#include "matrix.h"
#include "standard_matrices.h"
#include "mo.h"

class RHF {
public:
  RHF(standard_matrices&, MOs&);
	~RHF();
	bool get_convergency();
  matrix transform_matrix(const matrix&);
	void validate_convergency();
  void print_iteration();
  void calculate_fock();
  void calculate_pre_diis_error();
  void calculate_diis_coefs();
  void calculate_diis_fock();
  void calculate_diis_error();
  void update_buffer(std::deque<matrix>&, matrix);
  void calculate_density();
  void calculate_eri_matrix();
  void update_energy();
  void calculate_expansion();
  void core_guess();
  void direct_iteration();
	void roothan_hall_step();
	void roothan_hall();
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
  matrix density, eri_matrix, fock_matrix, error_matrix;
  std::deque<matrix> fock_buffer, error_buffer;
  std::vector<double> diis_coefs;
};

#endif
