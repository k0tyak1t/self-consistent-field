#include "matrix.h"
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <cmath>
#include <vector>

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

matrix::matrix(std::initializer_list<double> init_list)
{
    n = static_cast<int>(std::sqrt(init_list.size()));

    if (n * n != static_cast<int>(init_list.size())) {
        throw std::invalid_argument("Failed to create non-square matrix");
    }

    _matrix_elements = new double[n * n];
    std::copy(init_list.begin(), init_list.end(), _matrix_elements);
}

matrix::matrix(int init_size, double const* init_arr)
    : n(init_size)
{
    _matrix_elements = new double[n * n];
    memcpy(_matrix_elements, init_arr, n * n * sizeof(double));
}

matrix::~matrix() { delete[] _matrix_elements; }

double* matrix::operator[](int row) { return &_matrix_elements[row * n]; }

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
        throw std::runtime_error("Matrix already initialized.\n");
    }
    n = n_new;
    _matrix_elements = new double[n * n];
    return 0;
}

std::ostream& operator<<(std::ostream& os, const matrix& mat)
{
    for (int i = 0; i < mat.n; ++i) {
        for (int j = 0; j < mat.n; ++j) {
            os << mat[i][j] << " ";
        }
        os << '\n';
    }
    return os;
}

matrix matrix::zero_like(const matrix& other)
{
    matrix result = other;
    for (int i = 0; i < n * n; ++i) {
        _matrix_elements[i] = 0;
    }

    return result;
}

const int matrix::get_size() const { return n; }

void matrix::get_matrix_elements(double* dest_arr)
{
    if (n == 0) {
        throw std::runtime_error("Failed to get elements of not defined matrix!\n");
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest_arr[i * n + j] = _matrix_elements[i * n + j];
        }
    }
}

void matrix::from_array(const double* src)
{
    memcpy(_matrix_elements, src, n * n * sizeof(double));
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

void matrix::eigen_vv(double* evec, double* eval)
{
    char jobz, uplo;
    int lwork;
    int info;
    if (n == 0) {
        throw std::runtime_error(
            "Failed to get eigensystem of undefined matrix!\n");
    }
    lwork = 3 * n;
    double* work = new double[3 * n];
    get_matrix_elements(evec);
    jobz = 'V';
    uplo = 'U';
    dsyev_(&jobz, &uplo, &n, evec, &n, eval, work, &lwork, &info);
    if (info != 0) {
        throw std::runtime_error("Failed to diagonalize matrix!\n");
    }
    delete[] work;
}

matrix matrix::inv() const {
    if (n <= 0) {
        throw std::invalid_argument("Matrix size must be greater than zero");
    }

    double determinant = det(*this);
    if (determinant == 0) {
        throw std::runtime_error("Matrix is singular and cannot be inverted");
    }

    // Create identity matrix
    matrix identity(n);
    for (int i = 0; i < n; ++i) {
        identity[i][i] = 1.0;
    }

    // Step 1: LU Decomposition with partial pivoting
    std::vector<int> pivot(n);
    matrix LU = *this;
    for (int i = 0; i < n; ++i) pivot[i] = i;

    for (int k = 0; k < n; ++k) {
        // Find pivot row
        double max_val = 0.0;
        int max_row = k;
        for (int i = k; i < n; ++i) {
            if (std::abs(LU[pivot[i]][k]) > max_val) {
                max_val = std::abs(LU[pivot[i]][k]);
                max_row = i;
            }
        }
        if (max_val == 0.0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted");
        }

        // Swap rows in pivot
        std::swap(pivot[k], pivot[max_row]);

        // Perform LU factorization
        for (int i = k + 1; i < n; ++i) {
            LU[pivot[i]][k] /= LU[pivot[k]][k];
            for (int j = k + 1; j < n; ++j) {
                LU[pivot[i]][j] -= LU[pivot[i]][k] * LU[pivot[k]][j];
            }
        }
    }

    // Step 2: Invert the matrix using forward and backward substitution
    matrix inverse(n);
    for (int col = 0; col < n; ++col) {
        // Forward substitution: solve LY = I[:, col]
        std::vector<double> Y(n, 0.0);
        for (int i = 0; i < n; ++i) {
            Y[i] = identity[pivot[i]][col];
            for (int j = 0; j < i; ++j) {
                Y[i] -= LU[pivot[i]][j] * Y[j];
            }
        }

        // Backward substitution: solve UX = Y
        std::vector<double> X(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            X[i] = Y[i];
            for (int j = i + 1; j < n; ++j) {
                X[i] -= LU[pivot[i]][j] * X[j];
            }
            X[i] /= LU[pivot[i]][i];
        }

        // Set the column of the inverse matrix
        for (int i = 0; i < n; ++i) {
            inverse[i][col] = X[i];
        }
    }

    return inverse;
}


matrix& matrix::operator+=(const matrix& other)
{
    for (int i = 0; i < n * n; ++i) {
        _matrix_elements[i] += other._matrix_elements[i];
    }

    return *this;
}

matrix matrix::operator+(const matrix& other)
{
    matrix result = *this;
    result += other;
    return result;
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
        _matrix_elements[i] *= num;
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
            result[i][j] /= num;
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

double matrix::trace()
{
    double result = 0;
    for (int i = 0; i < n; ++i) {
        result += _matrix_elements[i * (1 + n)];
    }
    return result;
}

double trace(const matrix& mat)
{
    double result = 0;
    for (int i = 0; i < mat.n; ++i) {
        result += mat[i][i];
    }
    return result;
}

double frobenius_product(matrix& mat1, matrix& mat2)
{
    return trace(mat1.T() * mat2);
}

matrix matrix::minor(int i, int j) const
{
    matrix result(n - 1);
    int linear_idx = 0;
    for (int k = 0; k < n; ++k) {
        for (int l = 0; l < n; ++l) {
            if (k != i && l != j) {
                result._matrix_elements[linear_idx] = (*this)[k][l];
                ++linear_idx;
            }
        }
    }

    return result;
}

double det(const matrix& mat)
{
    if (mat.n == 0) {
        return 0;
    }

    if (mat.n == 1) {
        return mat[0][0];
    }

    if (mat.n == 2) {
        return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    }

    double result = 0.0;

    for (int i = 0; i < mat.n; ++i) {
        if (mat[0][i]) {
            result += mat[0][i] * (i % 2 == 0 ? 1 : -1) * det(mat.minor(0, i));
        }
    }

    return result;
}

std::vector<double> matrix::operator*(const std::vector<double>& vec) const
{
    std::vector<double> result;

    double ri = 0;
    for (int i = 0; i < n; ++i) {
        ri = 0;
        for (int j = 0; j < n; ++j) {
            ri += (*this)[i][j] * vec[j];
        }
        result.push_back(ri);
    }
    return result;
}