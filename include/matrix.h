#pragma once
#include <initializer_list>
#include <ostream>
#include <vector>

class matrix {
public:
  // constructors & destructor
  matrix();
  matrix(int);
  matrix(const matrix &);
  matrix(int, double const *);
  matrix(std::initializer_list<double>);
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
  int n;
  double *_matrix_elements;
};
