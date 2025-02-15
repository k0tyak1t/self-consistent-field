#pragma once
#include <deque>

template <typename T> struct Buffer : public std::deque<T> {
  Buffer(std::size_t size, T const &value = T{})
      : std::deque<T>(size, value) {};
  void update(const T &rhs) {
    this->pop_front();
    this->push_back(rhs);
  }
};
