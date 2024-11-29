#ifndef INPUT_AND_MD_INTEGRALS_H
#define INPUT_AND_MD_INTEGRALS_H

#include<vector>
#include<map>
#include<algorithm>
#include<iostream>
#include<fstream>
#include"standard_matrices.h"
#include"single_basis_function.h"

using std::map;
using std::vector;
using std::pair;
using std::make_pair;
using std::transform;
using std::string;
using std::cout;
using std::cerr;
using std::ifstream;
using std::getline;

class InputAndMDIntegrals{
public:
	bool readGeom(char* filename);
	bool readBasisLib(char* filename);
	bool calc (standard_matrices& );
	bool makeBasis();
	void printBasis();
	void test() const;
	bool set_c2v_z(int* );
private:
	double calcVnn();
	int get_Z(const std::string&) const;
	int get_L(const std::string&) const;
	struct Atom{
		double x,y,z;
		int q;
	};
	std::vector<Atom> atoms;
	map<int, vector<pair<int, vector<pair<double, double>>>>> basisLib;
	vector<vector<pair<int, SingleBasisFunction>>> basisFunctions;
};

#endif // INPUT_AND_MD_INTEGRALS_H
