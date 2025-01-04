#pragma once
#include <string>
#include <vector>

#include "mo.h"

void print_orbitals(const MO &);
void print_energies(const MO &);
void display_progress(int, const std::string &);
double dot_product(const std::vector<double> &, const std::vector<double> &);
