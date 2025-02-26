#include "mur_dav_prim_integrals.h"
#include <boost/math/special_functions/erf.hpp>
#include <vector>

int MurDavPrimIntegrals::calcEijt(double *E, int imax, int jmax, double p,
                                  double mu, double xPA, double xPB) {
  int tmax = imax + jmax + 1;
  std::vector<double> Etmp(tmax * tmax);

  for (int i = 0; i < tmax * tmax; i++) {
    Etmp[i] = 0;
  }

  Etmp[0] = exp(-mu * (xPA - xPB) * (xPA - xPB));

  for (int i = 0; i < imax; i++)
    for (int t = 0; t <= i; t++) {
      Etmp[(i + 1) * tmax + t] += Etmp[i * tmax + t] * xPA;
      Etmp[(i + 1) * tmax + t + 1] += Etmp[i * tmax + t] * 0.5 / p;
      Etmp[(i + 1) * tmax + t - 1] += Etmp[i * tmax + t] * t;
    }

  for (int j = 0; j < jmax; j++)
    for (int t = 0; t <= (j + imax); t++) {
      Etmp[(imax + j + 1) * tmax + t] += Etmp[(imax + j) * tmax + t] * xPB;
      Etmp[(imax + j + 1) * tmax + t + 1] +=
          Etmp[(imax + j) * tmax + t] * 0.5 / p;
      Etmp[(imax + j + 1) * tmax + t - 1] += Etmp[(imax + j) * tmax + t] * t;
    }

  for (int t = 0; t < tmax; t++)
    E[t] = Etmp[tmax * tmax - tmax + t];

  return 0;
}

int MurDavPrimIntegrals::calcEij3(double &Eij, double &Eij2p, double &Eij2m,
                                  int imax, int jmax, double p, double mu,
                                  double xPA, double xPB) {
  int tmax = imax + jmax + 3;
  std::vector<double> Etmp(tmax * tmax);

  for (int i = 0; i < tmax * tmax; i++)
    Etmp[i] = 0;
  Etmp[0] = exp(-mu * (xPA - xPB) * (xPA - xPB));
  for (int i = 0; i < imax; i++) {
    for (int t = 0; t <= i; t++) {
      Etmp[(i + 1) * tmax + t] += Etmp[i * tmax + t] * xPA;
      Etmp[(i + 1) * tmax + t + 1] += Etmp[i * tmax + t] * 0.5 / p;
      Etmp[(i + 1) * tmax + t - 1] += Etmp[i * tmax + t] * t;
    }
  }
  for (int j = 0; j < jmax + 2; j++) {
    for (int t = 0; t <= (j + imax); t++) {
      Etmp[(imax + j + 1) * tmax + t] += Etmp[(imax + j) * tmax + t] * xPB;
      Etmp[(imax + j + 1) * tmax + t + 1] +=
          Etmp[(imax + j) * tmax + t] * 0.5 / p;
      Etmp[(imax + j + 1) * tmax + t - 1] += Etmp[(imax + j) * tmax + t] * t;
    }
  }

  Eij = Etmp[(tmax - 3) * tmax];
  Eij2p = Etmp[(tmax - 1) * tmax];
  Eij2m = 0;

  if (jmax > 1) {
    Eij2m = Etmp[(tmax - 5) * tmax];
  }

  return 0;
}

int MurDavPrimIntegrals::calcRntuv(double *R, int tmax, int umax, int vmax,
                                   double p, double xPA, double yPA,
                                   double zPA) {
  const int nmax = tmax + umax + vmax,
            ndim = (tmax + 1) * (umax + 1) * (vmax + 1),
            tdim = (umax + 1) * (vmax + 1), udim = (vmax + 1);
  std::vector<double> Fn(nmax + 1);
  double x2 = p * (xPA * xPA + yPA * yPA + zPA * zPA);
  if (x2 < 1.e-4) {
    for (int n = 0; n < (nmax + 1); n++)
      Fn[n] = 1.0 / (2.0 * n + 1.0);
  } else {
    Fn[0] = erf(sqrt(x2)) * sqrt(acos(0.0) * 0.5 / x2);
    for (int n = 1; n < (nmax + 1); n++)
      Fn[n] = ((2 * n - 1) * Fn[n - 1] - exp(-x2)) * 0.5 / x2;
  }
  for (int n = 0; n <= nmax; n++) {
    R[n * ndim] = pow(-2.0 * p, n) * Fn[n];
  }

  int xR;
  if (tmax > 0) {
    for (int n = 0; n < nmax; n++)
      R[n * ndim + tdim] = xPA * R[(n + 1) * ndim];
    for (int t = 1; t < tmax; t++)
      for (int n = 0; n < nmax - t; n++) {
        xR = (n + 1) * ndim + t * tdim;
        R[xR - ndim + tdim] = xPA * R[xR] + t * R[xR - tdim];
      }
  }

  if (umax > 0) {
    for (int t = 0; t <= tmax; t++)
      for (int n = 0; n <= (nmax - t); n++)
        R[n * ndim + t * tdim + udim] = yPA * R[(n + 1) * ndim + t * tdim];

    for (int u = 1; u < umax; u++)
      for (int t = 0; t <= tmax; t++)
        for (int n = 0; n < nmax - t - u; n++) {
          xR = (n + 1) * ndim + t * tdim + udim * u;
          R[xR - ndim + udim] = yPA * R[xR] + u * R[xR - udim];
        }
  }

  if (vmax > 0) {
    for (int t = 0; t <= tmax; t++)
      for (int u = 0; u <= umax; u++)
        for (int n = 0; n <= (nmax - t - u); n++) {
          R[n * ndim + t * tdim + u * udim + 1] =
              zPA * R[(n + 1) * ndim + t * tdim + u * udim];
        }

    for (int v = 1; v < vmax; v++)
      for (int t = 0; t <= tmax; t++)
        for (int u = 0; u <= umax; u++)
          for (int n = 0; n < nmax - t - u - v; n++) {
            xR = (n + 1) * ndim + t * tdim + udim * u + v;
            R[xR - ndim + 1] = zPA * R[xR] + v * R[xR - 1];
          }
  }

  return 0;
}

double MurDavPrimIntegrals::Sij(int nx1, int ny1, int nz1, double x1, double y1,
                                double z1, double alpha1, int nx2, int ny2,
                                int nz2, double x2, double y2, double z2,
                                double alpha2) {

  const double p = alpha1 + alpha2, mu = alpha1 * alpha2 / p,
               xp = (alpha1 * x1 + alpha2 * x2) / p,
               yp = (alpha1 * y1 + alpha2 * y2) / p,
               zp = (alpha1 * z1 + alpha2 * z2) / p;

  std::vector<double> Eij(nx1 + nx2 + 1);
  std::vector<double> Ekl(ny1 + ny2 + 1);
  std::vector<double> Emn(nz1 + nz2 + 1);
  calcEijt(Eij.data(), nx1, nx2, p, mu, xp - x1, xp - x2);
  calcEijt(Ekl.data(), ny1, ny2, p, mu, yp - y1, yp - y2);
  calcEijt(Emn.data(), nz1, nz2, p, mu, zp - z1, zp - z2);
  const double S = pow(acos(0.0) * 2.0 / p, 1.5) * Eij[0] * Ekl[0] * Emn[0];
  return S;
}

double MurDavPrimIntegrals::Tij(int nx1, int ny1, int nz1, double x1, double y1,
                                double z1, double alpha1, int nx2, int ny2,
                                int nz2, double x2, double y2, double z2,
                                double alpha2) {
  const double p = alpha1 + alpha2, mu = alpha1 * alpha2 / p,
               xp = (alpha1 * x1 + alpha2 * x2) / p,
               yp = (alpha1 * y1 + alpha2 * y2) / p,
               zp = (alpha1 * z1 + alpha2 * z2) / p;
  double S_ij, S_ij2p, S_ij2m, S_kl, S_kl2p, S_kl2m, S_mn, S_mn2p, S_mn2m;
  calcEij3(S_ij, S_ij2p, S_ij2m, nx1, nx2, p, mu, xp - x1, xp - x2);
  calcEij3(S_kl, S_kl2p, S_kl2m, ny1, ny2, p, mu, yp - y1, yp - y2);
  calcEij3(S_mn, S_mn2p, S_mn2m, nz1, nz2, p, mu, zp - z1, zp - z2);
  double T_ij = -2 * alpha2 * alpha2 * S_ij2p + alpha2 * (2 * nx2 + 1) * S_ij -
                nx2 * (nx2 - 1) * S_ij2m / 2;
  double T_kl = -2 * alpha2 * alpha2 * S_kl2p + alpha2 * (2 * ny2 + 1) * S_kl -
                ny2 * (ny2 - 1) * S_kl2m / 2;
  double T_mn = -2 * alpha2 * alpha2 * S_mn2p + alpha2 * (2 * nz2 + 1) * S_mn -
                nz2 * (nz2 - 1) * S_mn2m / 2;
  return (T_ij * S_kl * S_mn + S_ij * T_kl * S_mn + S_ij * S_kl * T_mn) *
         pow(acos(0.0) * 2.0 / p, 1.5);
}
double MurDavPrimIntegrals::Vij(int nx1, int ny1, int nz1, double x1, double y1,
                                double z1, double alpha1, int nx2, int ny2,
                                int nz2, double x2, double y2, double z2,
                                double alpha2, double xq, double yq,
                                double zq) {
  const double p = alpha1 + alpha2, mu = alpha1 * alpha2 / p,
               xp = (alpha1 * x1 + alpha2 * x2) / p,
               yp = (alpha1 * y1 + alpha2 * y2) / p,
               zp = (alpha1 * z1 + alpha2 * z2) / p;

  const int tmax = nx1 + nx2 + 1;
  const int umax = ny1 + ny2 + 1;
  const int vmax = nz1 + nz2 + 1;

  std::vector<double> Rntuv(tmax * umax * vmax * (tmax + umax + vmax - 2));
  std::vector<double> Eij(tmax), Ekl(umax), Emn(vmax);

  calcRntuv(Rntuv.data(), tmax - 1, umax - 1, vmax - 1, p, xp - xq, yp - yq,
            zp - zq);
  calcEijt(Eij.data(), nx1, nx2, p, mu, xp - x1, xp - x2);
  calcEijt(Ekl.data(), ny1, ny2, p, mu, yp - y1, yp - y2);
  calcEijt(Emn.data(), nz1, nz2, p, mu, zp - z1, zp - z2);

  double vt = 0;
  for (int t = 0; t < tmax; t++)
    for (int u = 0; u < umax; u++)
      for (int v = 0; v < vmax; v++)
        vt += Eij[t] * Ekl[u] * Emn[v] * Rntuv[t * umax * vmax + u * vmax + v];
  return vt * acos(0.0) * 4.0 / p;
}

double MurDavPrimIntegrals::Vijkl(
    int nx1, int ny1, int nz1, double x1, double y1, double z1, double alpha1,
    int nx2, int ny2, int nz2, double x2, double y2, double z2, double alpha2,
    int nx3, int ny3, int nz3, double x3, double y3, double z3, double alpha3,
    int nx4, int ny4, int nz4, double x4, double y4, double z4, double alpha4) {
  const double p = alpha1 + alpha2, mup = alpha1 * alpha2 / p,
               xp = (alpha1 * x1 + alpha2 * x2) / p,
               yp = (alpha1 * y1 + alpha2 * y2) / p,
               zp = (alpha1 * z1 + alpha2 * z2) / p, q = alpha3 + alpha4,
               muq = alpha3 * alpha4 / q, xq = (alpha3 * x3 + alpha4 * x4) / q,
               yq = (alpha3 * y3 + alpha4 * y4) / q,
               zq = (alpha3 * z3 + alpha4 * z4) / q, al = p * q / (p + q);

  const int t1max = nx1 + nx2 + 1, u1max = ny1 + ny2 + 1, v1max = nz1 + nz2 + 1;

  std::vector<double> E1ij(t1max);
  std::vector<double> E1kl(u1max);
  std::vector<double> E1mn(v1max);
  calcEijt(E1ij.data(), nx1, nx2, p, mup, xp - x1, xp - x2);
  calcEijt(E1kl.data(), ny1, ny2, p, mup, yp - y1, yp - y2);
  calcEijt(E1mn.data(), nz1, nz2, p, mup, zp - z1, zp - z2);

  const int t2max = nx3 + nx4 + 1;
  const int u2max = ny3 + ny4 + 1;
  const int v2max = nz3 + nz4 + 1;
  std::vector<double> E2ij(t2max);
  std::vector<double> E2kl(u2max);
  std::vector<double> E2mn(v2max);
  calcEijt(E2ij.data(), nx3, nx4, q, muq, xq - x3, xq - x4);
  calcEijt(E2kl.data(), ny3, ny4, q, muq, yq - y3, yq - y4);
  calcEijt(E2mn.data(), nz3, nz4, q, muq, zq - z3, zq - z4);

  const int tmax = t1max + t2max - 1;
  const int umax = u1max + u2max - 1;
  const int vmax = v1max + v2max - 1;

  std::vector<double> Rntuv(tmax * umax * vmax * (tmax + umax + vmax - 2));

  calcRntuv(Rntuv.data(), tmax - 1, umax - 1, vmax - 1, al, xp - xq, yp - yq,
            zp - zq);

  double vt = 0;
  for (int t1 = 0; t1 < t1max; t1++)
    for (int u1 = 0; u1 < u1max; u1++)
      for (int v1 = 0; v1 < v1max; v1++)
        for (int t2 = 0; t2 < t2max; t2++)
          for (int u2 = 0; u2 < u2max; u2++)
            for (int v2 = 0; v2 < v2max; v2++) {
              vt += E1ij[t1] * E1kl[u1] * E1mn[v1] * E2ij[t2] * E2kl[u2] *
                    E2mn[v2] *
                    Rntuv[(t1 + t2) * umax * vmax + (u2 + u1) * vmax +
                          (v1 + v2)] *
                    pow(-1, v2 + u2 + t2);
            }

  return vt * 8 * acos(0.0) * acos(0.0) / p / q * sqrt(acos(0.0) * 2 / (p + q));
}

#if 0
double prim_integrals::Sij(const SingleBasisFunction &f1, const double alpha1,
                           const SingleBasisFunction &f2, const double alpha2) {

  return MurDavPrimIntegrals().Sij(f1.nx, f1.ny, f1.nz, f1.x, f1.y, f1.z,
                                   alpha1, f2.nx, f2.ny, f2.nz, f2.x, f2.y,
                                   f2.z, alpha2);
}

double prim_integrals::Tij(const SingleBasisFunction &f1, const double alpha1,
                           const SingleBasisFunction &f2, const double alpha2) {

  return MurDavPrimIntegrals().Tij(f1.nx, f1.ny, f1.nz, f1.x, f1.y, f1.z,
                                   alpha1, f2.nx, f2.ny, f2.nz, f2.x, f2.y,
                                   f2.z, alpha2);
}

double prim_integrals::Vijkl(const SingleBasisFunction &f1, const double alpha1,
                             const SingleBasisFunction &f2, const double alpha2,
                             const SingleBasisFunction &f3, const double alpha3,
                             const SingleBasisFunction &f4,
                             const double alpha4) {

  return MurDavPrimIntegrals().Vijkl(
      f1.nx, f1.ny, f1.nz, f1.x, f1.y, f1.z, alpha1, f2.nx, f2.ny, f2.nz, f2.x,
      f2.y, f2.z, alpha2, f3.nx, f3.ny, f3.nz, f3.x, f3.y, f3.z, alpha3, f4.nx,
      f4.ny, f4.nz, f4.x, f4.y, f4.z, alpha4);
}
#endif
