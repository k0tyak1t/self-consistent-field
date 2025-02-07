#pragma once

#include "utilities.hpp"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <type_traits>

// Declaration
namespace linalg {
template <typename T> class Matrix {
  typedef T value_type;
  struct Row; // proxy for accessors

public: // raw constructors and destructor
  ~Matrix();
  Matrix(std::size_t, std::size_t); // create matrix with fixed shape
  Matrix(std::size_t n)
      : Matrix(n, n) {};     // create squared matrix with fixed shape
  Matrix(const Matrix<T> &); // copy constructor

public: // static factory methods
  static Matrix<T> zero(std::size_t, std::size_t);
  static Matrix<T> zeros_like(const Matrix<T> &);
  static Matrix<T> identity(std::size_t);
  static Matrix<T> identity_like(const Matrix<T> &);
  static Matrix<T> transposed(const Matrix<T> &);
  template <typename It> static Matrix<T> diagonal(const It &);

public: // selectors
  const T *begin() const { return data; }
  const T *end() const { return data + ncols * nrows; }
  const Row operator[](std::size_t) const;

public: // modifiers
  T *begin() { return data; }
  T *end() { return data + ncols * nrows; }
  Row operator[](std::size_t);

public: // properties
  bool is_squared() const { return ncols == nrows; }
  std::size_t size() const { return ncols * nrows; }
  std::tuple<std::size_t, std::size_t> shape() const { return {nrows, ncols}; }
  std::size_t get_nrows() const { return nrows; }
  std::size_t get_ncols() const { return ncols; }

public: // static methods
  static T det(const Matrix<T> &);
  static T trace(const Matrix<T> &);

public: // non-const operations
  Matrix<T> &operator+=(const Matrix<T> &);
  Matrix<T> &operator-=(const Matrix<T> &);
  Matrix<T> &operator*=(const Matrix<T> &);

public: // const operations
  Matrix<T> operator+(const Matrix<T> &) const;
  Matrix<T> operator-(const Matrix<T> &) const;
  Matrix<T> operator*(const Matrix<T> &) const;
  T trace() const;
  T det() const; // TODO: implement LU

private: // implementation details
  void resize(std::size_t, std::size_t);

private: // fields and proxy
  std::size_t nrows = 0, ncols = 0;
  T *data = nullptr;

  struct Row {
    typedef T value_type;
    ~Row() = default;
    Row() = delete;
    Row(T *, std::size_t);
    T &operator[](std::size_t);
    const T &operator[](std::size_t) const;
    const T *begin() const { return row_data; }
    const T *end() const { return row_data + ncols; }
    T *begin() { return row_data; }
    T *end() { return row_data + ncols; }

  private:
    T *row_data;
    std::size_t ncols;
  };
};
} // namespace linalg

// Implementation
namespace linalg {
// internal constructors & destructor
template <typename T> Matrix<T>::~Matrix() { delete[] data; }

template <typename T>
Matrix<T>::Matrix(std::size_t nrows, std::size_t ncols)
    : nrows(nrows), ncols(ncols), data(new T[ncols * nrows]{}) {}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &other) : Matrix(other.nrows, other.ncols) {
  std::copy(other.begin(), other.end(), data);
}

// static factory methods
template <typename T>
Matrix<T> Matrix<T>::zero(std::size_t nrows, std::size_t ncols) {
  static_assert(std::is_arithmetic<T>::value,
                "Zero() can be used only with arithmetical types");
  return Matrix<T>{nrows, ncols};
}

template <typename T> Matrix<T> Matrix<T>::zeros_like(const Matrix<T> &source) {
  return Matrix<T>::zero(source.nrows, source.ncols);
}

template <typename T> Matrix<T> Matrix<T>::identity(std::size_t n) {
  static_assert(std::is_arithmetic<T>::value,
                "Identity() can be used only with arithmetical types");

  Matrix<T> identity = Matrix<T>::zero(n, n);
  for (auto i = 0; i < n; ++i) {
    identity.data[i * n + i] = 1;
  }
  return identity;
}

template <typename T>
Matrix<T> Matrix<T>::identity_like(const Matrix<T> &prototype) {
  if (!prototype.is_squared())
    throw std::invalid_argument("prototype matrix must be squared!\n");

  return Matrix<T>::identity(prototype.nrows);
}

template <typename T>
Matrix<T> Matrix<T>::transposed(const Matrix<T> &original) {
  Matrix<T> transposed = Matrix<T>{original.ncols, original.nrows};
  // TODO: implement with in-place transposition of copy

  for (auto i = 0; i < original.nrows; ++i)
    for (auto j = 0; j < original.ncols; ++j)
      transposed[j][i] = original[i][j];

  return transposed;
}

template <typename T>
template <typename It>
Matrix<T> Matrix<T>::diagonal(const It &src) {
  std::size_t n = src.size();
  Matrix<T> result{n, n};
  for (auto i = 0; i < n; ++i)
    result[i][i] = src[i];

  return result;
}

// selectors and modifiers

template <typename T>
typename Matrix<T>::Row Matrix<T>::operator[](std::size_t row) {
  if (row >= this->nrows)
    throw std::out_of_range(
        "Row index is out of range: " + std::to_string(row) + "\n");

  return Row(data + row * ncols, ncols);
}

template <typename T>
const typename Matrix<T>::Row Matrix<T>::operator[](std::size_t row) const {
  if (row >= this->nrows)
    throw std::out_of_range(
        "Row index is out of range: " + std::to_string(row) + "\n");

  return Row(data + row * ncols, ncols);
}

template <typename T>
Matrix<T>::Row::Row(T *row_data, std::size_t ncols)
    : row_data(row_data), ncols(ncols) {}

template <typename T> T &Matrix<T>::Row::operator[](std::size_t col) {
  if (col >= ncols)
    throw std::out_of_range(
        "Column index is out of range: " + std::to_string(col) + "\n");

  return row_data[col];
}

template <typename T>
const T &Matrix<T>::Row::operator[](std::size_t col) const {
  if (col >= ncols)
    throw std::out_of_range(
        "Column index is out of range: " + std::to_string(col) + "\n");

  return row_data[col];
}

// internal implementation
template <typename T>
void Matrix<T>::resize(std::size_t new_nrows, std::size_t new_ncols) {
  if (new_nrows * new_ncols != nrows * ncols)
    throw std::invalid_argument(
        "On resizing a number of elements shouldn't change!");

  nrows = new_nrows;
  ncols = new_ncols;
}

// non-const operations
template <typename T> Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &other) {
  for (auto &this_other_element : zip(*this, other)) {
    std::get<0>(this_other_element) += std::get<0>(this_other_element);
  }

  return *this;
}

template <typename T> Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &other) {
  for (auto &this_other_element : zip(*this, other)) {
    std::get<0>(this_other_element) -= std::get<0>(this_other_element);
  }

  return *this;
}

template <typename T> Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &other) {
  if (this->ncols != other.nrows)
    throw std::invalid_argument("Matrix must have matching dimentions");

  Matrix<T> result = Matrix(this->nrows, other.ncols);
  Matrix<T> other_T = transposed(other);

  for (auto i = 0; i < nrows; ++i)
    for (auto j = 0; j < other.ncols; ++j)
      result[i][j] = dot((*this)[i], other_T[j]);
  *this = result;

  return *this;
}

// const operations
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const {
  Matrix<T> result = *this;
  result += other;
  return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const {
  Matrix<T> result = *this;
  result -= other;
  return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const {
  Matrix<T> result{nrows, other.ncols};
  Matrix<T> other_T = transposed(other);

  for (auto i = 0; i < result.nrows; ++i)
    for (auto j = 0; j < result.ncols; ++j)
      result[i][j] = dot((*this)[i], other_T[j]);

  return result;
}

template <typename T> T Matrix<T>::det() const {
  throw std::runtime_error("det method not yet implemented!");
}

template <typename T> T Matrix<T>::trace() const {
  if (!this->is_squared())
    throw std::invalid_argument("Matrix must be squared to get trace of it");

  T trace{};
  for (int i = 0; i < nrows; ++i)
    trace += (*this)[i][i];
  return trace;
}

} // namespace linalg

#ifdef USE_NEW_MARTICES
using matrix = linalg::Matrix<double>;
#endif
