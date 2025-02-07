#ifndef USE_NEW_MATRICES
#pragma once
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

  // operators
  double *operator[](int);
  const double *operator[](int) const;
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
  matrix minor(int, int) const;
  friend double det(const matrix &);
  matrix T();
  void eigen_vv(double *, double *);
  matrix inverse() const;

  // initializers
  void zeroize();
  int init(const int n);
  void from_array(const double *);
  friend matrix zero_like(const matrix &);

  // ugly getters :)
  double *data();
  int get_size() const;

private:
  std::size_t n;
  double *data_;

private:
  class RowProxy {
  public: // selectors
    double operator[](std::size_t) const;
    const double *begin() const;
    const double *end() const;

  public: // modifiers
    double operator[](std::size_t);
    double *begin();
    double *end();

  private:
    std::size_t n;
    double *data;
  };
};
#endif
