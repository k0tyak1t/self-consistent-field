#pragma once
#include "atom.h"
#include <fstream>
#include <string>
#include <vector>

class InputManager {
public:
  ~InputManager();
  void read_geometry(const std::string &);
  std::vector<Atom> get_geometry();
  unsigned get_atomic_charge(std::string &);

private:
  std::ifstream input;
  std::vector<Atom> atoms;
};
