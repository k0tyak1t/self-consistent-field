#ifndef TEST
#include <iostream>
#include <stdexcept>

#include "diis.h"
#include "input_and_md_integrals.h"
#include "mo.h"
#include "standard_matrices.h"
#include <assert.h>
#include <cmath>

int main(int argc, char **argv) {
  char *mol_name = argv[1];
  char *basis = argv[2];

  if (argc < 2) {
    throw std::runtime_error("Filename should be specified!\n");
  }

  StandardMatrices sm;
  InputAndMDIntegrals loader;

  std::cout << "Preparing...\n";
  if (!loader.readGeom(mol_name) || !loader.readBasisLib(basis)) {
    throw std::runtime_error("Failed to load geometry or basis!\n");
  }

  loader.makeBasis();

  if (!loader.calc(sm)) {
    throw std::runtime_error("Failed to calculate standard matrices!\n");
  }

  MO mo(sm.get_nAO());

#ifdef NDIIS
  SCF rhf(mo, sm);
#else
  DIIS rhf(mo, sm);
#endif
  rhf.solve();

  return 0;
}

#endif // !TEST
