#pragma once
#include "matrix.h"

class MO {
public: // constructor
  MO(const unsigned = 0);

public:     // fields
  Matrix C; // MO to AO expansion. MO in columns

public: // getters
  double get_mo_energy(const int) const;
  int get_size() const;
  double get_total_energy() const;
  int get_irrep_characters(const int &) const;

public: // setters
  int set_mo_energy(const int, const double);
  int set_total_energy(const double);
  bool set_c2v(int *symmAO, const double &limit);
  void set_mo_energies(const double *);

private:         // fields
  std::size_t n; // basis dimentionality
  std::vector<double> mo_energies;
  std::vector<int> irrep;
  double total_energy;
};
