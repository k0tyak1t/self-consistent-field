#include "utils.h"

#include <iomanip>
#include <iostream>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <string>

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

/**
 * @brief print progress bar to the screen
 *
 * @param progress progress state for progress bar
 * @param leading_str string before prgoress bar
 */
void display_progress(const int progress, const std::string &leading_str) {
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
