#pragma once

#include <map>
#include <vector>

#include "atom.h"
#include "single_basis_function.h"
#include "standard_matrices.h"

using std::make_pair;
using std::map;
using std::pair;
using std::vector;

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
  vector<Atom> atoms;
  map<int, vector<pair<int, vector<pair<double, double>>>>> basisLib;
  vector<vector<pair<int, SingleBasisFunction>>> basisFunctions;
};
