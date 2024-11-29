#include <iostream>
#include <iomanip>
#include <string>
#include "standard_matrices.h"
#include "input_and_md_integrals.h"
#include "mo.h"
#include "output.h"
#include "rhf.h"

int main(int argc, char** argv){
  char* mol_name = argv[1];
  char* basis = argv[2];

  if (argc < 2) {
    throw std::runtime_error("Filename should be specified!\n");
  }

	standard_matrices sm;
	InputAndMDIntegrals a;
	MOs mo;

	if (!a.readGeom(mol_name) || !a.readBasisLib(basis)) {
    throw std::runtime_error("Failed to read load or basis!\n");
  }

	a.makeBasis();
	a.printBasis();

	if(!a.calc(sm)) {
    throw std::runtime_error("Failed to calculate standard matrices!\n");
  }

  RHF rhf(sm, mo);
  rhf.roothan_hall();

	//output().run(mo);
	output().print_energies_only(mo);
	return 0;
}
