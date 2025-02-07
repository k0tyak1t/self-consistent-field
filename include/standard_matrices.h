#pragma once
#include "matrix.h"

class StandardMatrices {
public: // members
  StandardMatrices(const unsigned int = 0);
  ~StandardMatrices();

  void init(const unsigned int);

  // getters
  double get_eri(int, int, int, int) const;
  double get_total_Vnn() const;
  int get_nAO() const;
  int get_num_el() const;

  // setters
  void set_eri(int, int, int, int, double);
  void set_total_Vnn(const double);
  void set_num_el(const unsigned int);

private: // fields
  int nAO;
  int num_el;
  double total_Vnn;
  double *eri;

public: // fields
  matrix S, T, H, Ven, X;
};
