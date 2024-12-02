#include "matrix.h"
#include <cstring>
#include <iomanip>
#include <iostream>

matrix::matrix()
    : n(0)
    , _matrix_elements(nullptr)
{
}

matrix::matrix(int n_new)
    : n(n_new)
{
    _matrix_elements = new double[n * n];
}

matrix::matrix(const matrix& other)
{
    n = other.n;
    _matrix_elements = new double[n * n];
    memcpy(_matrix_elements, other._matrix_elements, n * n * sizeof(double));
}

matrix::matrix(int init_size, double const* init_arr)
    : n(init_size)
{
    _matrix_elements = new double[n * n];
    memcpy(_matrix_elements, init_arr, n * n * sizeof(double));
}

matrix::~matrix()
{
    if (n > 0) {
        delete[] _matrix_elements;
        _matrix_elements = nullptr;
    };
}

double* matrix::operator[](int row)
{
    return &_matrix_elements[row * n];
}

const double* matrix::operator[](int row) const
{
    return &_matrix_elements[row * n];
}

void matrix::zeroize()
{
    for (int i = 0; i < n * n; ++i) {
        _matrix_elements[i] = 0;
    }
}

int matrix::init(const int n_new)
{
    if (n > 0) {
        throw std::runtime_error("** ERROR *** object Cl2 yet have n\n");
    };
    n = n_new;
    _matrix_elements = new double[n * n];
    return 0;
}

int matrix::print() const
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            std::cout << std::setw(10) << std::setprecision(6) << _matrix_elements[i * n + j] << ' ';
        std::cout << '\n';
    }
    return 0;
}

int matrix::get_size() const { return n; }

int matrix::check_ij(const int i, const int j, const char* name) const
{
    if (n == 0) {
        std::cerr << "** ERROR *** array not defined in " << name << '\n';
        return 1;
    }

    if ((i < 0) || (i > n - 1)) {
        std::cerr << "** ERROR *** invalid number of string in " << name << '\n';
        return 1;
    }

    if ((j < 0) || (j > n - 1)) {
        std::cerr << "** ERROR *** invalid number of element in " << name << '\n';
        return 1;
    }

    return 0;
}

double matrix::get_element(const int i, const int j) const
{
    if (check_ij(i, j, "get_element") > 0)
        return 1;
    return _matrix_elements[i * n + j];
}

int matrix::get_row(const int i, double* a) const
{
    if (check_ij(i, 0, "get_row") > 0)
        return 1;
    for (int j = 0; j < n; j++)
        a[j] = _matrix_elements[i * n + j];
    return 0;
}

int matrix::get_matrix_elements(double* a) const
{
    if (n == 0) {
        std::cerr << "** ERROR *** array not defined in "
                  << "get_matrix_elements" << '\n';
        return 1;
    };
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[i * n + j] = _matrix_elements[i * n + j];
    return 0;
}

int matrix::set_element(const int i, const int j, const double aij)
{
    if (check_ij(i, j, "set_element") > 0)
        return 1;
    _matrix_elements[i * n + j] = aij;
    return 0;
}

int matrix::set_row(const int i, const double* a)
{
    if (check_ij(i, 0, "set_row") > 0)
        return 1;
    for (int j = 0; j < n; j++)
        _matrix_elements[i * n + j] = a[j];
    return 0;
}

void matrix::from_array(const double* source_array)
{
    memcpy(_matrix_elements, source_array, n * n * sizeof(double));
}

matrix matrix::T()
{
    matrix result(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = (*this)[j][i];
        }
    }
    return result;
}

extern "C" {
void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
    double* work, int* lwork, int* info);
}

int matrix::eigen_vv(double* evec, double* eval)
{
    char jobz, uplo;
    int lwork;
    int info;
    if (n == 0) {
        std::cerr << "** ERROR *** matrix is not defined" << '\n';
        return 1;
    };
    lwork = 3 * n;
    double* work;
    work = new double[3 * n];
    get_matrix_elements(evec);
    jobz = 'V';
    uplo = 'U';
    dsyev_(&jobz, &uplo, &n, evec, &n, eval, work, &lwork, &info);
    if (info != 0) {
        std::cerr << "** ERROR *** matrix has not been diagonalized " << '\n';
        return 1;
    };
    delete[] work;
    return 0;
}

matrix matrix::operator+(const matrix& other)
{
    matrix result(n);
    for (int i = 0; i < n * n; ++i) {
        result._matrix_elements[i] = _matrix_elements[i] + other._matrix_elements[i];
    }
    return result;
}

matrix& matrix::operator+=(const matrix& other)
{
    for (int i = 0; i < n * n; ++i) {
        _matrix_elements[i] += other._matrix_elements[i];
    }

    return *this;
}

matrix matrix::operator-(const matrix& other)
{
    matrix result(n);
    for (int i = 0; i < n * n; ++i) {
        result._matrix_elements[i] = _matrix_elements[i] - other._matrix_elements[i];
    }
    return result;
}

matrix matrix::operator*(const matrix& other)
{
    matrix result(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = 0;
            for (int k = 0; k < n; ++k) {
                result[i][j] += (*this)[i][k] * other[k][j];
            }
        }
    }

    return result;
}

matrix matrix::operator*(const double num)
{
    for (int i = 0; i < n * n; ++i) {
        _matrix_elements[i] *= n;
    }
    return *this;
}

matrix& matrix::operator=(const matrix& other)
{
    if (this != &other) {
        delete[] _matrix_elements;
        n = other.n;
        _matrix_elements = new double[n * n];
        memcpy(_matrix_elements, other._matrix_elements, n * n * sizeof(double));
    }
    return *this;
}

matrix matrix::operator/(const double num)
{
    matrix result(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] /= n;
        }
    }
    return result;
}

void matrix::operator()(int new_size)
{
    if (n) {
        delete[] _matrix_elements;
    }
    n = new_size;
    _matrix_elements = new double[n * n];
}
