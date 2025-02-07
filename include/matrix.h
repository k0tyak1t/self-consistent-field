#ifndef USE_NEW_MATRICES
#pragma once

#include <cstddef>
#include <ostream>
#include <vector>

class matrix {
  class RowProxy;

public:
  // constructors & destructor
  matrix();
  matrix(std::size_t);
  matrix(const matrix &);
  matrix(std::size_t, double const *);
  ~matrix();

public: // static factory methods
  static matrix zeros(std::size_t);
  static matrix identity(std::size_t);
  static matrix zero_like(const matrix &);
  static matrix identity_like(const matrix &);
  template <typename It> static matrix from_array(const It &, const It &);

public: // selectors
  const RowProxy operator[](std::size_t) const;
  const double *begin() const { return data_; }
  const double *end() const { return data_ + n * n; }

public: // modifiers
  RowProxy operator[](std::size_t);
  double *begin() { return data_; }
  double *end() { return data_ + n * n; }

public: // operators
  matrix &operator+=(const matrix &);
  matrix operator+(const matrix &);
  matrix &operator-=(const matrix &);
  matrix operator-(const matrix &);
  matrix operator*(const matrix &);
  matrix operator*(const double);
  matrix operator/(const double);
  matrix &operator=(const matrix &);
  bool operator==(const matrix &);
  void operator()(int);
  friend std::ostream &operator<<(std::ostream &, const matrix &);
  std::vector<double> operator*(const std::vector<double> &) const;

  // operations
  friend double trace(const matrix &);
  double trace();
  friend double frobenius_product(matrix &, matrix &);
  matrix minor(std::size_t, std::size_t) const;
  friend double det(const matrix &);
  matrix T();
  void eigen_vv(double *, double *);
  matrix inverse() const;

  // initializers
  void zeroize();
  int init(const int n);
  void from_array(const double *);
  friend matrix zero_like(const matrix &);

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
#endif
