#ifndef SINGLEBASISFUNCTION_H
#define SINGLEBASISFUNCTION_H

#include <vector>

using std::make_pair;
using std::pair;
using std::vector;

class SingleBasisFunction {
public:
    void setLC(const int& nx_, const int& ny_, const int& nz_, const double& x_, const double& y_, const double& z_);
    double x, y, z;
    int nx, ny, nz;
    vector<pair<double, double>> ai_ci;
    void renorm_ci();
};

#endif
