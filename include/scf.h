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

  // getters
  double get_total_energy();

  // utilities
  Matrix transform_matrix(const Matrix &) const;
  void print_iter(int) const;

  // helpers
  void core_guess();
  void update_lcao_coefs();
  void update_density();
  void update_fock();
  void update_energy();

  // solver
  virtual void solve();

protected:
  double etol;
  int max_iter;
  double prev_energy, cur_energy;
  std::vector<double> mo_energies;
  MO &mo;
  StandardMatrices &std_m;
  std::size_t nAO{};
  Matrix density{}, fock{}, lcao_coefs{};
};
