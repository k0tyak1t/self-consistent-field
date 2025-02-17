#pragma once

#include "matrix.h"
#include "mo.h"
#include "standard_matrices.h"
#include <vector>

#ifndef DEFAUlT_ETOL
#define DEFAULT_ETOL 1e-12
#endif

#ifndef DEFAULT_MAX_ITER
#define DEFAULT_MAX_ITER 100
#endif

class SCF {
public:
  SCF(MO &, StandardMatrices &);

public: // selectors
  double get_total_energy();

public: // utils
  Matrix transform_matrix(const Matrix &) const;
  void print_iter(int) const;

public: // solver subfunctions
  Matrix core_guess();
  void update_lcao_coefs();
  Matrix update_density() const;
  Matrix update_fock() const;
  void update_energy();

public: // solver
  virtual void solve();

protected:
  double etol;
  unsigned max_iter;
  double prev_energy, cur_energy;
  StandardMatrices &std_m;
  std::size_t nAO;
  MO &mo;
  std::vector<double> mo_energies;
  Matrix density, fock, lcao_coefs;
};
