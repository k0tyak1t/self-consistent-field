#include "utils.h"
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

void print_orbitals(const MO& mo) { std::cout << "\n-- Orbitals --\n"
                                              << mo.C; }

void print_energies(const MO& mo)
{
    std::cout << "\n-- Final energies --\n";
    int nAO = mo.C.get_size();
    for (int i = 0; i < nAO; ++i) {
        std::cout << std::setw(8) << std::setprecision(8) << mo.get_mo_energy(i)
                  << ' ';
    }
    std::cout << '\n';
    std::cout << "\nTotal_energy = " << mo.get_total_energy() << '\n';
}

void display_progress(int progress, const std::string& leading_str)
{
    int progress_bar_width = 30;
    int cursor_position = progress_bar_width * progress / 100;
    std::cout << std::setw(20) << leading_str;
    for (int i = 0; i < progress_bar_width; ++i) {
        std::cout << (i <= cursor_position ? "█" : "▒");
        std::cout.flush();
    }
    std::cout << " " << std::setw(5) << int(progress) << " %\r";
    std::cout.flush();

    if (progress == 100) {
        std::cout << std::endl;
    }
}

double dot_product(const std::vector<double>& vec_1,
    const std::vector<double>& vec_2)
{
    if (vec_1.size() != vec_2.size()) {
        throw std::runtime_error(
            "\nFailed to multiply vectors with different length!\n");
    }
    double result = 0;
    for (int i = 0; i < (int)vec_1.size(); ++i) {
        result += vec_1[i] * vec_2[i];
    }
    return result;
}
