#ifndef OUTPUT
#define OUTPUT_H
#include "mo.h"
#include "standard_matrices.h"
#include<iostream>
#include<iomanip>

class output {
public:	
	~output();
	output();
	int run(const MOs&) const;
	void print_energies_only(const MOs&) const;
};
#endif
