#include "matrix.hpp"

int main() {
  linalg::Matrix<double> a{2, 2};
  auto b = linalg::Matrix<double>::identity(2);
  auto c = a * b;
  linalg::print2d(a);
  linalg::print2d(b);
  linalg::print2d(c);
  auto d = linalg::Matrix<double>::diagonal(std::vector<double>{1, 2, 3});
  linalg::print2d(d);
  std::cout << "trace: " << d.trace() << std::endl;
}
