#include "matrix_buffer.h"

Buffer::Buffer(std::size_t size)
    : std::deque<Matrix>(size, Matrix::zero(size)) {}

void Buffer::update(const Matrix &rhs) {
  this->pop_front();
  this->push_back(rhs);
}
