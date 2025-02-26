#pragma once
#include <vector>

struct SingleBasisFunction {
  SingleBasisFunction(const int nx = 0, const int ny = 0, const int nz = 0,
                      const int x = 0, const int y = 0, const int z = 0)
      : nx(nx), ny(ny), nz(nz), x(x), y(y), z(z) {};

  void setLC(const int nx_, const int ny_, const int nz_, const double x_,
             const double y_, const double z_);

  void renorm_ci();

  int nx, ny, nz;
  double x, y, z;
  std::vector<std::pair<double, double>> ai_ci;
};
