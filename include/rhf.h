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
	void validate_convergency();
  void calculate_density();
  void calculate_eri_matrix();
  void calculate_fock_transformed();
  void recalculate_energy();
  void calculate_coef_matrix();
  void verbose_iteration();
  void core_guess();
  void direct_iteration();
	void roothan_hall();

private:
	double etol;
	int max_iter, diis_size, iter;
	bool is_converged;
	standard_matrices& std_m;
	MOs& mo;
  double prev_energy, cur_energy;
  double* evec;
  double* eval;
  matrix density, eri_matrix, fock_matrix;
  //std::deque<matrix> error_buffer;
  //std::vector<double> diis_coefs;
};

#endif
