#pragma once
#include "matrix.h"
#include <deque>

struct Buffer : std::deque<Matrix> {
  Buffer(std::size_t = 0);

  void update(const Matrix &);
};
