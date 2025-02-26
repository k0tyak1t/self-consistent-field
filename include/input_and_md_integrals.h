#pragma once

#include <map>
#include <vector>

#include "atom.h"
#include "single_basis_function.h"
#include "standard_matrices.h"

class InputAndMDIntegrals {
public:
  bool readGeom(char *filename);
  bool readBasisLib(char *filename);
  bool calc(StandardMatrices &);
  bool makeBasis();
  void printBasis();
  bool set_c2v_z(int *);

private:
  double calc_Vnn();
  int get_nuclear_charge(const std::string &) const;
  int get_orbital_momentum(const std::string &) const;
  std::vector<Atom> atoms;
  std::map<int,
           std::vector<std::pair<int, std::vector<std::pair<double, double>>>>>
      basisLib;
  std::vector<std::vector<std::pair<int, SingleBasisFunction>>> basisFunctions;
};
