#include "single_basis_function.h"
#include <cmath>
#include <vector>

using std::vector;

void SingleBasisFunction::setLC(const int& nx_, const int& ny_, const int& nz_, const double& x_, const double& y_, const double& z_)
{
    x = x_;
    y = y_;
    z = z_;
    nx = nx_;
    ny = ny_;
    nz = nz_;
}

void SingleBasisFunction::renorm_ci()
{
    for (vector<pair<double, double>>::iterator i = ai_ci.begin(); i != ai_ci.end(); i++) {
        double a = (*i).first;
        double p = 1;
        double q = 1;
        for (int j = 1; j <= nx; j++) {
            p *= (2 * j - 1);
            q *= 4 * a;
        };
        for (int j = 1; j <= ny; j++) {
            p *= (2 * j - 1);
            q *= 4 * a;
        };
        for (int j = 1; j <= nz; j++) {
            p *= (2 * j - 1);
            q *= 4 * a;
        };
        (*i).second /= sqrt(p / q * pow(M_PI * 0.5 / a, 1.5));
    }
}
