#ifndef RHF_H
#define RHF_H
#include <vector>
#include "matrix.h"
#include "standard_matrices.h"
#include "mo.h"

class RHF {
public:
  RHF(standard_matrices&, MOs&);
	~RHF();
	bool get_convergency();
  void validate_diis_condition(const int);
	void validate_convergency(const int);
  void build_density();
  void calculate_eri_matrix();
  void calculate_fock_transformed();
  double recalculate_energy();
  void calculate_coef_matrix();
  void verbose_iteration(const int);
  void core_guess();
	void roothan_hall();
  void calculate_error();
  void diis();

private:
	double etol;
	int max_iter, diis_size;
	bool is_converged, diis_mode;
	standard_matrices& std_m;
	MOs& mo;
  double E_old, E_new;
  double* evec;
  double* eval;
  matrix density, eri_matrix, fock_matrix;
  std::vector<matrix> error_buffer;
};

#endif
