#include "matrix.h"
#include <iostream>
#include <string>

void wrap_test(int callback, std::string test_name)
{
    std::cout << test_name + "test: ";
    if (callback) {
        std::cout << "failed!\n";
    } else {
        std::cout << "passed!\n";
    }
}

int test_matrix_init_list()
{
    try {
        matrix mat = { 1, 2, 3, 4 };
        std::cout << mat;
    } catch (std::exception& e) {
        return 1;
    }
    return 0;
}

int test_det()
{
    try {
        matrix mat = {
            1, 2, 3,
            2, 1, 3,
            3, 1, 2
        };

        std::cout << mat << std::endl
                  << det(mat) << std::endl;
    } catch (std::exception& e) {
        return 1;
    }
    return 0;
}

int test_inv()
{
    try {
        matrix mat = {
            1, 0, 0,
            0, 2, 0,
            0, 0, 3
        };

        std::cout << mat << "\n\n"
                  << mat.inv();
    } catch (const std::exception& e) {
        return 1;
    }

    return 0;
}

int main()
{
    wrap_test(test_matrix_init_list(), "matrix(init_list)");
    wrap_test(test_det(), "det");
    wrap_test(test_inv(), "inv");
    return 0;
}