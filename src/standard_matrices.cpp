#include "standard_matrices.h"
#include "matrix.h"

StandardMatrices::~StandardMatrices() {
  if (nAO > 0) {
    delete[] eri;
  }
}

StandardMatrices::StandardMatrices(std::size_t nAO) : nAO(nAO) {
  if (!nAO) {
    return;
  }

  S = matrix(nAO);
  T = matrix(nAO);
  H = matrix(nAO);
  Ven = matrix(nAO);
  eri = new double[nAO * nAO * nAO * nAO];
}

void StandardMatrices::init(std::size_t n_new) {
  if (nAO) {
    throw std::runtime_error("AO already exists!");
  }

  nAO = n_new;
  S = matrix{nAO};
  T = matrix{nAO};
  H = matrix{nAO};
  S = matrix{nAO};
  Ven = matrix{nAO};
  eri = new double[n_new * n_new * n_new * n_new]{};
}

void StandardMatrices::set_total_Vnn(const double new_total_Vnn) {
  total_Vnn = new_total_Vnn;
}

double StandardMatrices::get_total_Vnn() const { return total_Vnn; }

std::size_t StandardMatrices::get_nAO() const { return nAO; }

double StandardMatrices::get_eri(std::size_t i, std::size_t j, std::size_t k,
                                 std::size_t l) const {
  if (((i < 0) || (i >= nAO)) || ((j < 0) || (j >= nAO)) ||
      ((k < 0) || (k >= nAO)) || ((l < 0) || (l >= nAO))) {
    throw std::runtime_error("Failed to access element!");
  }

  return eri[i + nAO * j + nAO * nAO * k + nAO * nAO * nAO * l];
}

void StandardMatrices::set_eri(std::size_t i, std::size_t j, std::size_t k,
                               std::size_t l, double Vijkl) {
  if (((i < 0) || (i >= nAO)) || ((j < 0) || (j >= nAO)) ||
      ((k < 0) || (k >= nAO)) || ((l < 0) || (l >= nAO))) {
    throw std::runtime_error("Failed to access element!");
  }
  eri[i + nAO * j + nAO * nAO * k + nAO * nAO * nAO * l] = Vijkl;
}

std::size_t StandardMatrices::get_num_el() const { return num_el; }

void StandardMatrices::set_num_el(std::size_t new_num_el) {
  if (new_num_el % 2) {
    throw std::runtime_error(
        "Failed to solve RHF with odd number of electrons!");
  }

  num_el = new_num_el;
}
