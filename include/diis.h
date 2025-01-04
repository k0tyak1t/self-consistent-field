#pragma once

#include "scf.h"
#include "standard_matrices.h"
#include <deque>
#include <vector>

#ifndef DEFAULT_DIIS_SIZE
#define DEFAULT_DIIS_SIZE 5
#endif

class DIIS : public SCF {
public:
  DIIS(MO &, StandardMatrices &);

  // pre-diis
  void update_error();

  // diis
  void update_fock_buffer();
  void update_error_buffer();
  void update_extended_error_product();
  void update_diis_coefs();
  void update_diis_error();
  void update_diis_fock();

  // solver
  void solve() override;

private:
  int diis_size;
  matrix error, extended_diis_product;
  std::deque<matrix> fock_buffer, error_buffer;
  std::vector<double> diis_coefs;
};
