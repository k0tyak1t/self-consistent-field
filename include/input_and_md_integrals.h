#ifndef INPUT_AND_MD_INTEGRALS_H
#define INPUT_AND_MD_INTEGRALS_H

#include "single_basis_function.h"
#include "standard_matrices.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

using std::cerr;
using std::cout;
using std::getline;
using std::ifstream;
using std::make_pair;
using std::map;
using std::pair;
using std::string;
using std::transform;
using std::vector;

class InputAndMDIntegrals {
public:
    bool readGeom(char* filename);
    bool readBasisLib(char* filename);
    bool calc(standard_matrices&);
    bool makeBasis();
    void printBasis();
    bool set_c2v_z(int*);

private:
    double calcVnn();
    int get_nuclear_charge(const std::string&) const;
    int get_orbital_momentum(const std::string&) const;
    struct Atom {
        int q;
        double x, y, z;
    };
    vector<Atom> atoms;
    map<int, vector<pair<int, vector<pair<double, double>>>>> basisLib;
    vector<vector<pair<int, SingleBasisFunction>>> basisFunctions;
};

#endif // INPUT_AND_MD_INTEGRALS_H
