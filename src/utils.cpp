#include "utils.h"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

void print_orbitals(const MO &mo) { std::cout << "\n-- Orbitals --\n" << mo.C; }

void print_energies(const MO &mo) {
  std::cout << "\n-- Final energies --\n";
  int nAO = mo.C.size();
  for (int i = 0; i < nAO; ++i) {
    std::cout << std::setw(8) << std::setprecision(8) << mo.get_mo_energy(i)
              << ' ';
  }
  std::cout << '\n';
  std::cout << "\nTotal_energy = " << mo.get_total_energy() << '\n';
}

void display_progress(int progress, const std::string &leading_str) {
  int progress_bar_width = 30;
  int cursor_position = progress_bar_width * progress / 100;
  std::cout << std::setw(20) << std::left << leading_str;
  for (int i = 0; i < progress_bar_width; ++i) {
    std::cout << (i <= cursor_position ? "█" : "▒");
    std::cout.flush();
  }
  std::cout << " " << int(progress) << "%\r";
  std::cout.flush();

  if (progress == 100) {
    std::cout << std::endl;
  }
}

template <typename ItA, typename ItB>
double dot_product(const ItA &first_start, const ItA &first_end,
                   const ItB &second_start, const ItB &second_end) {
  if (std::distance(first_start, first_end) !=
      std::distance(second_start, second_end))
    throw std::invalid_argument(
        "Sizes of containers must be the same to find dot product!");
  double dot{};
  for (auto first = first_start, second = second_start; first != first_end;
       ++first, ++second)
    dot += first * second;

  return dot;
}
