#ifndef MATRIX_H
#define MATRIX_H
#include <initializer_list>
#include <ostream>

class matrix {
public:
    matrix();
    matrix(int);
    matrix(const matrix&);
    matrix();
    matrix(int, double const*);
    matrix(std::initializer_list<double>);
    ~matrix();
    double* operator[](int);
    const double* operator[](int) const;
    matrix T();
    void eigen_vv(double*, double*);
    matrix inv(const matrix&);
    int init(const int n); // TODO: deprecate.
    int get_size() const;
    void get_matrix_elements(double*);
    void from_array(const double*);
    friend std::ostream& operator<<(std::ostream&, const matrix&);
    void zeroize();
    matrix zero_like(const matrix&);
    matrix operator+(const matrix&);
    matrix& operator+=(const matrix&);
    matrix operator-(const matrix&);
    matrix operator*(const matrix&);
    matrix operator*(const double);
    matrix operator/(const double);
    matrix& operator=(const matrix&);
    void operator()(int);
    double trace();
    friend double trace(const matrix&);
    friend double frobenius_product(matrix&, matrix&);
    double minor(int, int) const;
    friend double det(const matrix&);

private:
    int n;
    double* _matrix_elements;
};
#endif
