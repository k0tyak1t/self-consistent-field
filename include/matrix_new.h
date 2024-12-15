#ifndef MATRIXv2_H
#define MATRIXv2_H
#include <vector>
#include <initializer_list>
#include <tuple>
#include <stdexcept>

template <typename T>
class Matrix{
public:
    Matrix(const int, const int, const std::initializer_list<T>&);
    Matrix(const int, const int);
    Matrix(const int);
    Matrix();
    Matrix<T>& operator=(const std::initializer_list<T>&);
    std::tuple<int, int> get_size() const noexcept;
    T& operator[](const std::tuple<int, int>&);
    const T operator[](const std::tuple<int, int>&) const;
    Matrix<T>& operator+=(const Matrix<T>&);
    Matrix<T> operator+(const Matrix<T>&) const;
    Matrix<T> operator-() const noexcept;
    Matrix<T>& operator-=(const Matrix<T>&);
    Matrix<T> operator-(const Matrix<T>&) const;
    Matrix<T>& operator*=(const Matrix<T>&);
    Matrix<T> operator*(const Matrix<T>&) const;
    void set_row(const int, const std::vector<T>&);
    void set_col(const int, const std::vector<T>&);
    std::tuple<Matrix<T>, Matrix<T>> get_eigensystem();
private:
    int nrows, ncols;
    std::vector<T> data;
};

#endif // MATRIXv2_H

template <typename T>
Matrix<T>::Matrix(const int nrows, const int ncols, const std::initializer_list<T>& items) : Matrix(nrows, ncols)
{
    auto item = items.begin();
    for (int i = 0; i < ncols * nrows; ++i, ++item) {
        data[i] = *item;
    }
}

template <typename T>
Matrix<T>::Matrix(const int nrows, const int ncols)
    : nrows(nrows), ncols(ncols)
{
    data = std::vector<T>(ncols * nrows, 0);
}

template <typename T>
Matrix<T>::Matrix(const int n): Matrix(n, n) {}

template <typename T>
Matrix<T>::Matrix(): Matrix(0, 0) {}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const std::initializer_list<T>& items)
{
    if (ncols * nrows != items.size()) {
        throw std::runtime_error("Failed to assign matrix elements with non-matching size of elements");
    }
    auto item = items.begin();
    for (int i = 0; i < data.size(); ++i, ++item) {
        data[i] = *item;
    }

    return *this;
}

template<typename T>
std::tuple<int, int> Matrix<T>::get_size() const noexcept
{
    return std::tuple<int, int>{nrows, ncols};
}

template <typename T>
T& Matrix<T>::operator[](const std::tuple<int, int>& ij)
{
    int i = std::get<0>(ij), j = std::get<1>(ij);
    return data.at(nrows * i + j);
}

template <typename T>
const T Matrix<T>::operator[](const std::tuple<int, int>& ij) const
{
    int i = std::get<0>(ij), j = std::get<1>(ij);
    return data.at(nrows * i + j);
}

template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other)
{
    if (nrows != other.nrows || ncols != other.ncols) {
        throw std::runtime_error("Failed to add matrices with different sizes!");
    }

    for (auto this_it = data.begin(), other_it = other.data.begin();
              this_it != data.end(); ++this_it, ++other_it) {
        
        *this_it += *other_it;
    }

    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const
{
    Matrix<T> result = *this;
    result += other;
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-() const noexcept
{
    Matrix<T> result = *this;
    for (auto& matrix_element : result.data) {
        matrix_element = -matrix_element;
    }

    return result;
} 

template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other)
{
    if (nrows != other.nrows || ncols != other.ncols) {
        throw std::runtime_error("Failed to add matrices with different sizes!");
    }

    for (auto this_it = data.begin(), other_it = other.data.begin();
              this_it != data.end(); ++this_it, ++other_it) {
        
        *this_it -= *other_it;
    }

    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const
{
    Matrix<T> result = *this;
    result -= other;
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const
{
    if (ncols != other.nrows) {
        throw std::runtime_error("Failed to multiply matrices with non-matching dimentions!");
    }

    Matrix<T> result(nrows, other.ncols);

    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < other.ncols; ++j) {
            for (int k = 0; k < ncols; ++k) {
                result[{i, j}] += (*this)[{i, k}] * other[{k, j}];
            }
        }
    }

    return result;
}

template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& other)
{
    *this = (*this) * other;
    return *this;
}

template <typename T>
void Matrix<T>::set_row(const int i, const std::vector<T>& new_row)
{
    if (new_row.size() != ncols) {
        throw std::runtime_error("Failed to set row with non-matching size!");
    }

    if (i < 0 || i >= nrows) {
        throw std::runtime_error("Incorrect index: " + std::to_string(i));
    }
}

template <typename T>
void Matrix<T>::set_col(const int j, const std::vector<T>& new_col)
{
    if (new_col.size() != nrows) {
        throw std::runtime_error("Failed to set column with non-matching size!");
    }

    if (j < 0 || j >= ncols) {
        throw std::runtime_error("Incorrect index: " + std::to_string(j));
    }
}

template <typename T>
std::tuple<Matrix<T>, Matrix<T>> Matrix<T>::get_eigensystem()
{
    if (ncols != nrows) {
        throw std::runtime_error("Failed to solve eigensystem for non-squared matrix!");
    }
}