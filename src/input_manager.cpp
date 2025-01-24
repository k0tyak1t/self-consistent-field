#include "input_manager.h"
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

InputManager::~InputManager() {
  if (input.is_open()) {
    input.close();
  }
}

void InputManager::read_geometry(const std::string &filename) {
  unsigned n_atoms = 0;
  std::string cur_line;
  input.open(filename);
  input >> n_atoms;
  std::getline(input, cur_line); // comment line skip

  for (unsigned i = 0; i < n_atoms; ++i) {
    Atom cur_atom;
    input >> cur_line >> cur_atom.x >> cur_atom.y >> cur_atom.z;
    cur_atom.q = get_atomic_charge(cur_line);
    atoms.push_back(cur_atom);
  }

  input.close();
}

unsigned InputManager::get_atomic_charge(std::string &element) {

  std::transform(element.begin(), element.end(), element.begin(), ::toupper);

  if (element == "H" || element == "HYDROGEN" || element == "1")
    return 1;
  if (element == "HE" || element == "HELIUM" || element == "2")
    return 2;
  if (element == "LI" || element == "LITHIUM" || element == "3")
    return 3;
  if (element == "BE" || element == "BERILIUM" || element == "4")
    return 4;
  if (element == "B" || element == "BORON" || element == "5")
    return 5;
  if (element == "C" || element == "CARBON" || element == "6")
    return 6;
  if (element == "N" || element == "NITROGEN" || element == "7")
    return 7;
  if (element == "O" || element == "OXYGEN" || element == "8")
    return 8;
  if (element == "F" || element == "FLUORINE" || element == "9")
    return 9;

  throw std::runtime_error("Failed to read element with sign: " + element);
}
