#include "input_and_md_integrals.h"
#include "mo.h"
#include "rhf.h"
#include "standard_matrices.h"
#include "utils.h"
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    char* mol_name = argv[1];
    char* basis = argv[2];

    if (argc < 2) {
        throw std::runtime_error("Filename should be specified!\n");
    }

    standard_matrices sm;
    InputAndMDIntegrals loader;
    MOs mo;

    if (!loader.readGeom(mol_name) || !loader.readBasisLib(basis)) {
        throw std::runtime_error("Failed to load geometry or basis!\n");
    }

    loader.makeBasis();

    if (!loader.calc(sm)) {
        throw std::runtime_error("Failed to calculate standard matrices!\n");
    }

    RHF rhf(sm, mo);
    rhf.roothan_hall();

    print_energies(mo);
    return 0;
}
