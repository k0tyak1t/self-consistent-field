#include "matrix_new.h"
#include <string>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <iomanip>

void wrapper(const int callback, const std::string& msg)
{
    std::cout << std::left << std::setw(35) << msg + " test";
    if (callback) {
        std::cout << "failed!\n";
        return;
    }

    std::cout <<  "passed!\n";
}

int test_int_int_constructor()
{
    try
    {
        Matrix<double> mat(2, 2);
    }
    catch(const std::exception& e)
    {
        return 1;
    }
    return 0;
} 

int test_int_constructor()
{
    try
    {
        Matrix<double> mat(2);
    }
    catch(const std::exception& e)
    {
        return 1;
    }

    return 0;
}

int test_default_constructor()
{
    try
    {
        Matrix<double> mat;
    }
    catch(const std::exception& e)
    {
        return 1;
    }

    return 0;
}

int test_getter()
{
    try
    {
        Matrix<double> mat(2, 2);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                if (mat[{i, j}] != 0) {
                    return 1;
                } 
            }
        }
    }
    catch(const std::exception& e)
    {
        return 1;
    }

    return 0;
}

int test_setter()
{
    try
    {
        Matrix<double> mat(2, 2);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                mat[{i, j}] = i + j;
            }
        }
        
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                if (mat[{i, j}] != i + j) {
                    return 1;
                }
            }
        }

        return 0;
    }
    catch(const std::exception& e)
    {
        return 1;
    }

    return 0;
}

int test_default_assignement()
{
    Matrix<double> mat1(2, 2);
    try
    {
        Matrix<double> mat2 = mat1;

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                if (mat2[{i, j}] != mat1[{i, j}]) {return 1;}
            }
        }
    }
    catch(const std::exception& e)
    {
        return 1;
    }

    return 0;
}

int test_init_list_assignement()
{

    Matrix<double> mat(2);
    try
    {
        mat = {1, 2, 3, 4};
    }

    catch(const std::exception& e)
    {
        return 1;
    }

    return 0;
    
}

int test_sum()
{
    try
    {
        Matrix<double> mat1(2, 2);
        mat1 = {1, 2, 3, 4};
        Matrix<double> mat2(2, 2);
        mat2 = {4, 3, 2, 1};
        Matrix<double> result = mat1 + mat2;

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                if (result[{i,j}] != mat1[{i, j}] + mat2[{i, j}]) {
                    return 1;
                }
            }
        }
    }
    catch(const std::exception& e)
    {
        return 1;
    }

    return 0;
}

template<typedef return_type>
void new_wrapper(return_type (*test)(), const std::string& name)
{
    std::cout << std::left << std::setw(35) << name + " test";
    try
    {
        test();
        std::cout << "passed!";
    }
    catch(std::exception& e)
    {
        std::cout << "failed!";
    }
}

int test_init_list_constructor()
{
    try {
        Matrix<double> mat(3, 2) = {1, 2, 3, 4, 5, 6};
    }
    catch (std::exception& e) {
        return 1;
    }

    return 0;
}

int main()
{
    wrapper(test_int_int_constructor(), "Matrix(int, int)");
    wrapper(test_int_constructor(), "Matrix(int)");
    wrapper(test_default_constructor(), "Matrix()");
    wrapper(test_getter(), "Matrix getter");
    wrapper(test_setter(), "Matrix setter");
    wrapper(test_default_assignement(), "Matrix assignement");
    wrapper(test_init_list_assignement(), "init_list assignement");
    wrapper(test_sum(), "sum");
    return 0;
}