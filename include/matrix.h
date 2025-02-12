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
  const RowProxy operator[](std::size_t) const;
  const double *begin() const { return data_; }
  const double *end() const { return data_ + n * n; }

public: // modifiers
  RowProxy operator[](std::size_t);
  double *begin() { return data_; }
  double *end() { return data_ + n * n; }

public:
  Matrix &operator=(const Matrix &); // copy-assignment
  Matrix &operator=(Matrix &&);      // move-assignment

public: // non-const operations
  Matrix &operator+=(const Matrix &);
  Matrix &operator-=(const Matrix &);
  Matrix &operator*=(double);
  Matrix &operator/=(double);
  Matrix transpose(); // in-place transposition

public: // const operations
  Matrix operator+(const Matrix &) const;
  Matrix operator-(const Matrix &) const;
  Matrix operator*(const Matrix &) const;
  Matrix operator*(double) const;
  Matrix operator/(double) const;
  bool operator==(const Matrix &) const;
  double trace() const;
  Matrix inverse() const;
  void eigen_vv(double *, double *) const;
  std::vector<double> operator*(const std::vector<double> &) const;

public: // static functions:
  static Matrix transposed(const Matrix &);
  static Matrix inversed(const Matrix &);
  static double det(const Matrix &);
  static double dot(const Matrix &, const Matrix &);

public:
  Matrix minor(std::size_t, std::size_t) const;
  Matrix transposed();

  // initializers
  void zeroize();
  int init(const int n);
  void from_array(const double *);

  // getters
  double *data() { return data_; }
  const double *data() const { return data_; }
  std::size_t size() const { return n; }

private: // fields
  std::size_t n;
  double *data_;

private: // proxy for getters
  class RowProxy {
  public: // constructors, default destructor
    RowProxy() : n(0), data_(nullptr) {}
    RowProxy(std::size_t n, double *data_) : n(n), data_(data_) {}

  public: // selectors
    const double &operator[](std::size_t) const;
    const double *begin() const { return data_; }
    const double *end() const { return data_ + n * n; }

  public: // modifiers
    double &operator[](std::size_t);
    double *begin() { return data_; }
    double *end() { return data_ + n * n; }

  private: // fields
    std::size_t n;
    double *data_;
  };
};

std::ostream &operator<<(std::ostream &os, const Matrix &mat);
