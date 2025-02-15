#pragma once
#include "matrix.h"

class StandardMatrices {
public: // members
  StandardMatrices(std::size_t = 0);
  ~StandardMatrices();

  void init(std::size_t);

  // getters
  double get_eri(std::size_t, std::size_t, std::size_t, std::size_t) const;
  double get_total_Vnn() const;
  std::size_t get_nAO() const;
  std::size_t get_num_el() const;

  // setters
  void set_eri(std::size_t, std::size_t, std::size_t, std::size_t, double);
  void set_total_Vnn(double);
  void set_num_el(std::size_t);

private: // fields
  std::size_t nAO;
  std::size_t num_el;
  double total_Vnn;
  double *eri;

public: // fields
  Matrix S, T, H, Ven, X;
};
