#include "matrix.h"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <vector>

extern "C" {
// LU decomoposition of a general matrix
void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
// generate inverse of a matrix given its LU decomposition
void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork,
             int *INFO);
}

// * * * * ------------ * * * * //
// * * * * constructors * * * * //
// * * * * ------------ * * * * //
matrix::matrix() : n(0), data_(nullptr) {}

matrix::matrix(std::size_t n_new) : n(n_new) { data_ = new double[n * n]; }

matrix::matrix(const matrix &other) {
  n = other.n;
  data_ = new double[n * n];
  memcpy(data_, other.data_, n * n * sizeof(double));
#ifndef NDEBUG
  std::cout << "matrix created!\n";
#endif
}

matrix::matrix(std::size_t init_size, double const *init_arr) : n(init_size) {
  data_ = new double[n * n];
  memcpy(data_, init_arr, n * n * sizeof(double));
}

matrix::~matrix() {
  delete[] data_;
#ifndef NDEBUG
  std::cout << "matrix destroyed!\n";
#endif // NDEBUG
}

// * * * * ********* * * * * //
// * * * * operators * * * * //
// * * * * ********* * * * * //

matrix &matrix::operator+=(const matrix &other) {
  for (auto i = 0; i < n * n; ++i) {
    data_[i] += other.data_[i];
  }

  return *this;
}

matrix matrix::operator+(const matrix &other) {
  matrix result = *this;
  result += other;
  return result;
}

matrix &matrix::operator-=(const matrix &other) {
  for (int i = 0; i < n; ++i) {
    data_[i] -= other.data_[i];
  }
  return *this;
}

matrix matrix::operator-(const matrix &other) {
  matrix result = *this;
  result -= other;
  return result;
}

matrix matrix::operator*(const matrix &other) {
  matrix result(n);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      result[i][j] = 0;
      for (int k = 0; k < n; ++k)
        result[i][j] += (*this)[i][k] * other[k][j];
    }

  return result;
}

matrix matrix::operator*(const double num) {
  matrix result(n);
  for (int i = 0; i < n * n; ++i)
    result.data_[i] *= num;

  return result;
}

matrix matrix::operator/(const double num) {
  matrix result(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result[i][j] /= num;
    }
  }
  return result;
}

matrix &matrix::operator=(const matrix &other) {
  if (this != &other) {
    delete[] data_;
    n = other.n;
    data_ = new double[n * n];
    memcpy(data_, other.data_, n * n * sizeof(double));
  }
  return *this;
}

void matrix::operator()(int new_size) {
  if (n) {
    delete[] data_;
  }
  n = new_size;
  data_ = new double[n * n];
}

std::ostream &operator<<(std::ostream &os, const matrix &mat) {
  for (int i = 0; i < mat.n; ++i) {
    for (int j = 0; j < mat.n; ++j) {
      os << std::setprecision(4) << std::setw(8) << mat[i][j] << " ";
    }
    os << '\n';
  }
  return os;
}

std::vector<double> matrix::operator*(const std::vector<double> &vec) const {
  std::vector<double> result;

  double ri = 0;
  for (int i = 0; i < n; ++i) {
    ri = 0;
    for (int j = 0; j < n; ++j) {
      ri += (*this)[i][j] * vec[j];
    }
    result.push_back(ri);
  }
  return result;
}

bool matrix::operator==(const matrix &other) {
  for (int i = 0; i < n * n; ++i) {
    if (data_[i] != other.data_[i])
      return false;
  }

  return true;
}

// * * * * ---------- * * * * //
// * * * * operations * * * * //
// * * * * ---------- * * * * //
double trace(const matrix &mat) {
  double result = 0;
  for (int i = 0; i < mat.n; ++i) {
    result += mat[i][i];
  }
  return result;
}

double matrix::trace() {
  double result = 0;
  for (int i = 0; i < n; ++i) {
    result += data_[i * (1 + n)];
  }
  return result;
}

double frobenius_product(matrix &mat1, matrix &mat2) {
  double result = 0;
  for (int i = 0; i < mat1.n; ++i)
    for (int k = 0; k < mat1.n; ++k)
      result += mat1[k][i] * mat2[k][i];

  return result;
}

matrix matrix::minor(int i, int j) const {
  matrix result(n - 1);
  int linear_idx = 0;
  for (int k = 0; k < n; ++k) {
    for (int l = 0; l < n; ++l) {
      if (k != i && l != j) {
        result.data_[linear_idx] = (*this)[k][l];
        ++linear_idx;
      }
    }
  }

  return result;
}

double det(const matrix &mat) {
  if (mat.n == 0) {
    return 0;
  }

  if (mat.n == 1) {
    return mat[0][0];
  }

  if (mat.n == 2) {
    return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
  }

  double result = 0.0;

  for (int i = 0; i < mat.n; ++i) {
    if (mat[0][i]) {
      result += mat[0][i] * (i % 2 == 0 ? 1 : -1) * det(mat.minor(0, i));
    }
  }

  return result;
}

matrix matrix::T() {
  matrix result(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      result[i][j] = (*this)[j][i];
    }
  }
  return result;
}

extern "C" {
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
            double *work, int *lwork, int *info);
}

void matrix::eigen_vv(double *evec, double *eval) {
  char jobz, uplo;
  int lwork;
  int info;
  if (n == 0) {
    throw std::runtime_error("Failed to get eigensystem of undefined matrix!");
  }
  lwork = 3 * n;
  double *work = new double[3 * n];
  memcpy(evec, this->data(), sizeof(double) * n * n);
  jobz = 'V';
  uplo = 'U';
  dsyev_(&jobz, &uplo, &n, evec, &n, eval, work, &lwork, &info);
  if (info != 0) {
    throw std::runtime_error("Failed to diagonalize matrix!\n");
  }
  delete[] work;
}

matrix matrix::inverse() const {
  matrix inv = *this;
  int N = n;
  int *IPIV = new int[n];
  int LWORK = n * n;
  double *WORK = new double[LWORK];
  int INFO;
  dgetrf_(&N, &N, inv.data(), &N, IPIV, &INFO);
  dgetri_(&N, inv.data(), &N, IPIV, WORK, &LWORK, &INFO);
  delete[] IPIV;
  delete[] WORK;
  return inv;
}

// * * * * ------------ * * * * //
// * * * * initializers * * * * //
// * * * * ------------ * * * * //
void matrix::zeroize() {
  for (int i = 0; i < n * n; ++i) {
    data_[i] = 0;
  }
}

int matrix::init(const int n_new) {
  if (n > 0) {
    throw std::runtime_error("Matrix already initialized.\n");
  }
  n = n_new;
  data_ = new double[n * n];
  return 0;
}

void matrix::from_array(const double *src) {
  memcpy(data_, src, n * n * sizeof(double));
}

matrix zero_like(const matrix &other) {
  matrix result = other;
  int n = other.n;
  for (int i = 0; i < n * n; ++i) {
    result.data_[i] = 0;
  }

  return result;
}

// * * * * ------------ * * * * //
// * * * *   getters    * * * * //
// * * * * ------------ * * * * //
double *matrix::data() { return data_; }

int matrix::get_size() const { return n; }
