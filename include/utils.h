#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <string>
#include "mo.h"
#include <vector>
#include <deque>

void print_orbitals(const MOs&);
void print_energies(const MOs&);
void display_progress(int, const std::string&);
double scalar_product(const std::vector<double>&, const std::vector<double>&);

#endif // UTILS_H
