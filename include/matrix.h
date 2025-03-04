#pragma once

#include <cstddef>
#include <ostream>
#include <vector>

class Matrix {
  class RowProxy;

public:
  // constructors & destructor
  Matrix();                  // default constructor
  Matrix(const std::size_t); // direct constructor
  Matrix(const Matrix &);    // copy constructor
  Matrix(Matrix &&);         // move constructor
  ~Matrix();                 // destructor

public: // static factory methods
  static Matrix zero(const std::size_t);
  static Matrix identity(const std::size_t);
  static Matrix zero_like(const Matrix &);
  static Matrix identity_like(const Matrix &);

public: // selectors
  const double operator()(const unsigned i, const unsigned j) const {
    return data_[i * n + j];
  }
  const double *begin() const { return data_; }
  const double *end() const { return data_ + n * n; }

public: // modifiers
  double &operator()(const unsigned i, const unsigned j) {
    return data_[i * n + j];
  }
  double *begin() { return data_; }
  double *end() { return data_ + n * n; }

public:
  Matrix &operator=(const Matrix &); // copy-assignment
  Matrix &operator=(Matrix &&);      // move-assignment

public: // non-const operations
  Matrix &operator+=(const Matrix &);
  Matrix &operator-=(const Matrix &);
  Matrix &operator/=(double);
  Matrix transpose(); // in-place transposition

public: // const operations
  Matrix operator+(const Matrix &) const;
  Matrix operator-(const Matrix &) const;
  Matrix operator*(const Matrix &) const;
  Matrix operator*(double) const;
  Matrix operator/(double) const;
  bool operator==(const Matrix &) const;
  void eigen_vv(double *, double *) const;

public: // static functions:
  static Matrix transposed(const Matrix &);
  static Matrix inversed(const Matrix &);
  static double dot(const Matrix &, const Matrix &);
  static double trace(const Matrix &);

public:
  Matrix minor(std::size_t, std::size_t) const;

  // initializers
  int init(const int n);
  void from_array(const double *);

  // getters
  double *data() { return data_; }
  const double *data() const { return data_; }
  std::size_t size() const { return n; }

private: // fields
  std::size_t n;
  double *data_;
};

std::ostream &operator<<(std::ostream &os, const Matrix &mat);
