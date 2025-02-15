#pragma once

#include "buffer.h"
#include "scf.h"
#include "standard_matrices.h"
#include <deque>
#include <vector>

#ifndef DEFAULT_DIIS_SIZE
#define DEFAULT_DIIS_SIZE 5
#endif

class DIIS : public SCF {
public: // constructor
  DIIS(MO &, StandardMatrices &);

public: // pre-diis stage
  void update_error();

public: // diis stage
  void update_extended_error_product();
  void update_diis_coefs();
  void update_diis_error();
  void update_diis_fock();
  void update_diis_density();

  // solver
  void solve() override;

private:
  std::size_t diis_size;
  Matrix error, extended_diis_product;
  Buffer<Matrix> fock_buffer, error_buffer, density_buffer;
  std::vector<double> diis_coefs;
};
