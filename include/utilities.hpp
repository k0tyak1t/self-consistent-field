#include <iostream>
#include <iterator>
#include <stdexcept>

namespace linalg {
template <typename ItA, typename ItB, typename T = typename ItA::value_type>
T dot(ItA first, ItB second) {
  if (std::distance(first.begin(), first.end()) !=
      std::distance(second.begin(), second.end()))
    throw std::invalid_argument(
        "Arguments or dot product must have the same length!");
  T result{};
  for (auto it1 = first.begin(), it2 = second.begin(); it2 != second.end();
       it1++, it2++)
    result += (*it1) * (*it2);
  return result;
}

template <typename It1d> void print1d(const It1d &cont1d) {
  for (auto &i : cont1d)
    std::cout << i << ' ';
}

template <typename It2d> void print2d(const It2d &matrix) {
  auto nrows = matrix.get_nrows();
  for (auto i = 0; i < nrows; ++i) {
    print1d(matrix[i]);
    std::cout << '\n';
  }
  std::cout << '\n';
}
} // namespace linalg
