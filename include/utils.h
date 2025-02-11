#pragma once
#include <string>

#include "mo.h"

void print_orbitals(const MO &);
void print_energies(const MO &);
void display_progress(int, const std::string &);
template <typename ItA, typename ItB>
double dot_product(const ItA &, const ItA &, const ItB &, const ItB &);
