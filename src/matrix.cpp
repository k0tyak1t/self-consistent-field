#include "matrix.h"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <vector>
#define CHECK_RANGES 1

extern "C" {
// LU decomoposition of a general matrix
void dgetrf_(int *M, int *N, double *A, int *lda, int *IPIV, int *INFO);
// generate inverse of a matrix given its LU decomposition
void dgetri_(int *N, double *A, int *lda, int *IPIV, double *WORK, int *lwork,
             int *INFO);
}

const double &Matrix::RowProxy::operator[](std::size_t i) const {
#ifndef NDEBUG
#if CHECK_RANGES
  if (i >= n || i < 0)
    throw std::out_of_range("Column index is out of range! " +
                            std::to_string(i) + " > " + std::to_string(n) +
                            '\n');

#endif // CHECK_RANGES
#endif // NDEBUG
  return data_[i];
}

double &Matrix::RowProxy::operator[](std::size_t i) {
#ifndef NDEBUG
#if CHECK_RANGES
  if (i >= n || i < 0)
    throw std::out_of_range("Column index is out of range! " +
                            std::to_string(i) + " > " + std::to_string(n) +
                            '\n');

#endif // CHECK_RANGES
#endif // NDEBUG
  return data_[i];
}

// * * * * ------------ * * * * //
// * * * * constructors * * * * //
// * * * * ------------ * * * * //
Matrix::Matrix() : n(0), data_(nullptr) {
#ifndef NDEBUG
  std::cout << "[matrix]: default initialization.\n";
#endif
}

Matrix::Matrix(const std::size_t n_new) : n(n_new), data_(new double[n * n]{}) {
#ifndef NDEBUG
  std::cout << "[matrix]: direct initialization.\n";
#endif
}

Matrix::Matrix(const Matrix &other) {
  n = other.n;
  data_ = new double[n * n]{};
  std::copy(other.begin(), other.end(), data_);
#ifndef NDEBUG
  std::cout << "[matrix]: copy initialization.\n";
#endif
}

Matrix::Matrix(Matrix &&other) {
  std::swap(n, other.n);
  std::swap(data_, other.data_);
#ifndef NDEBUG
  std::cout << "[matrix]: move initialization\n";
#endif
}

Matrix::~Matrix() {
#ifndef NDEBUG
  std::cout << "[matrix]: destructor.\n";
#endif // NDEBUG
  delete[] data_;
}

// selectors
const Matrix::RowProxy Matrix::operator[](std::size_t i) const {
  if (n <= i || i < 0)
    throw std::out_of_range("Row index is out of range!");

  return RowProxy{n, data_ + i * n};
}

// modifiers
Matrix::RowProxy Matrix::operator[](std::size_t i) {
  if (n <= i || i < 0)
    throw std::out_of_range("Row index is out of range!");

  return RowProxy{n, data_ + i * n};
}

Matrix &Matrix::operator=(Matrix &&other) {
#ifndef NDEBUG
  std::cout << "[matrix]: move assignment.\n";
#endif

  if (this == &other)
    return *this;

  std::swap(n, other.n);
  std::swap(data_, other.data_);
  return *this;
}

Matrix &Matrix::operator+=(const Matrix &other) {
  for (std::size_t i = 0; i < n * n; ++i) {
    data_[i] += other.data_[i];
  }

  return *this;
}

Matrix Matrix::operator+(const Matrix &other) const {
  Matrix result = *this;
  result += other;
  return result;
}

Matrix &Matrix::operator-=(const Matrix &other) {
  for (std::size_t i = 0; i < n; ++i) {
    data_[i] -= other.data_[i];
  }
  return *this;
}

Matrix Matrix::operator-(const Matrix &other) const {
  Matrix result = *this;
  result -= other;
  return result;
}

Matrix Matrix::operator*(const Matrix &other) const {
  Matrix result{n};

  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      result[i][j] = 0;
      for (std::size_t k = 0; k < n; ++k)
        result[i][j] += (*this)[i][k] * other[k][j];
    }

  return result;
}

Matrix Matrix::operator*(double num) const {
  Matrix result = *this;
  for (auto &i : result)
    i *= num;

  return result;
}

Matrix Matrix::transpose() {
  for (auto i = 0u; i < n; ++i)
    for (auto j = i; j < n; ++j)
      std::swap((*this)[i][j], (*this)[j][i]);

  return *this;
}

Matrix Matrix::operator/(const double num) const {
  Matrix result(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      result[i][j] /= num;
    }
  }
  return result;
}

Matrix &Matrix::operator=(const Matrix &other) {
  if (this == &other)
    return *this;

  delete[] data_;
  n = other.n;
  data_ = new double[n * n];
  std::copy(other.begin(), other.end(), data_);
  return *this;
}

std::ostream &operator<<(std::ostream &os, const Matrix &mat) {
  for (std::size_t i = 0; i < mat.size(); ++i) {
    for (std::size_t j = 0; j < mat.size(); ++j) {
      os << std::setprecision(4) << std::setw(8) << mat[i][j] << " ";
    }
    os << '\n';
  }
  return os;
}

std::vector<double> Matrix::operator*(const std::vector<double> &vec) const {
  std::vector<double> result;

  double ri = 0;
  for (std::size_t i = 0; i < n; ++i) {
    ri = 0;
    for (std::size_t j = 0; j < n; ++j) {
      ri += (*this)[i][j] * vec[j];
    }
    result.push_back(ri);
  }
  return result;
}

bool Matrix::operator==(const Matrix &other) const {
  for (std::size_t i = 0; i < n * n; ++i)
    if (data_[i] != other.data_[i])
      return false;

  return true;
}

double Matrix::trace() const {
  double result{};
  for (std::size_t i = 0; i < n; ++i)
    result += data_[i * (1 + n)];
  return result;
}

double Matrix::dot(const Matrix &mat1, const Matrix &mat2) {
  double result{};
  if (mat1.size() != mat2.size())
    throw std::invalid_argument(
        "Matrices must have the same number of elements.");

  for (auto i = 0u; i < mat1.n; ++i)
    result += mat1.data()[i] * mat2.data()[i];

  return result;
}

Matrix Matrix::minor(std::size_t i, std::size_t j) const {
  Matrix result(n - 1);
  int linear_idx = 0;
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t l = 0; l < n; ++l)
      if (k != i && l != j) {
        result.data_[linear_idx] = (*this)[k][l];
        ++linear_idx;
      }

  return result;
}

double Matrix::det(const Matrix &mat) { // good luck with O(n!) complexity :)
  if (mat.n == 0)
    return 0;

  if (mat.n == 1)
    return mat[0][0];

  if (mat.n == 2)
    return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];

  double result = 0.0;

  for (std::size_t i = 0; i < mat.n; ++i)
    if (mat[0][i])
      result += mat[0][i] * (i % 2 == 0 ? 1 : -1) * det(mat.minor(0, i));

  return result;
}

Matrix Matrix::transposed() {
  Matrix result = *this;
  result.transpose();
  return result;
}

extern "C" {
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
            double *work, int *lwork, int *info);
}

void Matrix::eigen_vv(double *evec, double *eval) const {
  if (n == 0)
    throw std::runtime_error("Failed to get eigensystem of undefined matrix!");

  char jobz = 'V', uplo = 'U';
  int N = n;
  int lwork = 3 * N;
  int info;
  double *work = new double[3 * N];
  std::copy(this->begin(), this->end(), evec);
  dsyev_(&jobz, &uplo, &N, evec, &N, eval, work, &lwork, &info);
  if (info != 0) {
    throw std::runtime_error("Failed to diagonalize matrix! lapack info: " +
                             std::to_string(info) + '\n');
  }
  delete[] work;
}

Matrix Matrix::inverse() const {
  Matrix inv{*this};
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

// methods to be dead...

void Matrix::zeroize() {
  for (std::size_t i = 0; i < n * n; ++i)
    data_[i] = 0;
}

int Matrix::init(const int n_new) {
  if (n > 0)
    throw std::runtime_error("Matrix already initialized.\n");
  n = n_new;
  data_ = new double[n * n]{};
  return 0;
}

void Matrix::from_array(const double *src) {
  std::copy(src, src + n * n, data_);
}

// static factory

Matrix Matrix::zero_like(const Matrix &other) { return Matrix{other.n}; }

Matrix Matrix::zero(const std::size_t n) { return Matrix{n}; }

Matrix Matrix::identity(const std::size_t n) {
  Matrix identity{n};
  for (auto i = 0u; i < n; ++i)
    identity[i][i] = 1.0;
  return identity;
}

Matrix Matrix::identity_like(const Matrix &reference) {
  return identity(reference.n);
}

// static methods

Matrix Matrix::transposed(const Matrix &matrix) {
  Matrix transposed = matrix;
  transposed.transpose();
  return transposed;
}
