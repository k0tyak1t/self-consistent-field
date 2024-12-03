#ifndef MATRIX_H
#define MATRIX_H
class matrix {
public:
    matrix();
    matrix(int);
    matrix(const matrix&);
    matrix(int, double const*);
    ~matrix();
    double* operator[](int);
    const double* operator[](int) const;
    matrix T();
    int eigen_vv(double*, double*);
    int init(const int n); // TODO: deprecate.
    int get_size() const;
    int check_ij(const int, const int, const char*) const; // legacy
    int get_matrix_elements(double*) const; // legacy
    int get_row(const int, double*) const; // legacy
    double get_element(const int, const int) const; // legacy
    void from_array(const double*);
    int set_row(const int, const double*); // legacy
    int set_element(const int, const int, const double); // legacy
    int print() const; // TODO: change return type.
    void zeroize();
    // matrix zero_like(const matrix&); // TODO: implement.
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
    friend double frobenius_product(const matrix&, const matrix&);

private:
    int n;
    double* _matrix_elements;
};
#endif
