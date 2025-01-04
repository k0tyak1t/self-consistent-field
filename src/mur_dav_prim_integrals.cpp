#include "mur_dav_prim_integrals.h"

#include <boost/math/special_functions/erf.hpp>

int MurDavPrimIntegrals::calcEijt(double *E, const int &imax, const int &jmax,
                                  const double &p, const double &mu,
                                  const double &xPA, const double &xPB) {
  int tmax = imax + jmax + 1;
  double *Etmp = new double[tmax * tmax];

  for (int i = 0; i < tmax * tmax; i++) {
    Etmp[i] = 0;
  }

  Etmp[0] = exp(-mu * (xPA - xPB) * (xPA - xPB));

  for (int i = 0; i < imax; i++) {
    for (int t = 0; t <= i; t++) {
      Etmp[(i + 1) * tmax + t] += Etmp[i * tmax + t] * xPA;
      Etmp[(i + 1) * tmax + t + 1] += Etmp[i * tmax + t] * 0.5 / p;
      Etmp[(i + 1) * tmax + t - 1] += Etmp[i * tmax + t] * t;
    }
  }

  for (int j = 0; j < jmax; j++) {
    for (int t = 0; t <= (j + imax); t++) {
      Etmp[(imax + j + 1) * tmax + t] += Etmp[(imax + j) * tmax + t] * xPB;
      Etmp[(imax + j + 1) * tmax + t + 1] +=
          Etmp[(imax + j) * tmax + t] * 0.5 / p;
      Etmp[(imax + j + 1) * tmax + t - 1] += Etmp[(imax + j) * tmax + t] * t;
    }
  }

  for (int t = 0; t < tmax; t++) {
    E[t] = Etmp[tmax * tmax - tmax + t];
  }

  delete[] Etmp;
  return 0;
}

int MurDavPrimIntegrals::calcEij3(double &Eij, double &Eij2p, double &Eij2m,
                                  const int &imax, const int &jmax,
                                  const double &p, const double &mu,
                                  const double &xPA, const double &xPB) {
  int tmax = imax + jmax + 3;
  double *Etmp = new double[tmax * tmax];

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

  delete[] Etmp;
  return 0;
}

int MurDavPrimIntegrals::calcRntuv(double *R, const int &tmax, const int &umax,
                                   const int &vmax, const double &p,
                                   const double &xPA, const double &yPA,
                                   const double &zPA) {
  const int nmax = tmax + umax + vmax;
  const int ndim = (tmax + 1) * (umax + 1) * (vmax + 1);
  const int tdim = (umax + 1) * (vmax + 1);
  const int udim = (vmax + 1);
  double *Fn = new double[nmax + 1];
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

  delete[] Fn;

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
    for (int t = 0; t <= tmax; t++) {
      for (int u = 0; u <= umax; u++) {
        for (int n = 0; n <= (nmax - t - u); n++) {
          R[n * ndim + t * tdim + u * udim + 1] =
              zPA * R[(n + 1) * ndim + t * tdim + u * udim];
        }
      }
    }

    for (int v = 1; v < vmax; v++) {
      for (int t = 0; t <= tmax; t++) {
        for (int u = 0; u <= umax; u++) {
          for (int n = 0; n < nmax - t - u - v; n++) {
            xR = (n + 1) * ndim + t * tdim + udim * u + v;
            R[xR - ndim + 1] = zPA * R[xR] + v * R[xR - 1];
          }
        }
      }
    }
  }

  return 0;
}

double MurDavPrimIntegrals::Sij(const int &nx1, const int &ny1, const int &nz1,
                                const double &x1, const double &y1,
                                const double &z1, const double &alpha1,
                                const int &nx2, const int &ny2, const int &nz2,
                                const double &x2, const double &y2,
                                const double &z2, const double &alpha2) {
  const double p = alpha1 + alpha2;
  const double mu = alpha1 * alpha2 / p;
  const double xp = (alpha1 * x1 + alpha2 * x2) / p;
  const double yp = (alpha1 * y1 + alpha2 * y2) / p;
  const double zp = (alpha1 * z1 + alpha2 * z2) / p;
  double *Eij = new double[nx1 + nx2 + 1];
  double *Ekl = new double[ny1 + ny2 + 1];
  double *Emn = new double[nz1 + nz2 + 1];
  calcEijt(Eij, nx1, nx2, p, mu, xp - x1, xp - x2);
  calcEijt(Ekl, ny1, ny2, p, mu, yp - y1, yp - y2);
  calcEijt(Emn, nz1, nz2, p, mu, zp - z1, zp - z2);
  const double S = pow(acos(0.0) * 2.0 / p, 1.5) * Eij[0] * Ekl[0] * Emn[0];
  delete[] Eij;
  delete[] Ekl;
  delete[] Emn;
  return S;
}

double MurDavPrimIntegrals::Tij(const int &nx1, const int &ny1, const int &nz1,
                                const double &x1, const double &y1,
                                const double &z1, const double &alpha1,
                                const int &nx2, const int &ny2, const int &nz2,
                                const double &x2, const double &y2,
                                const double &z2, const double &alpha2) {
  const double p = alpha1 + alpha2;
  const double mu = alpha1 * alpha2 / p;
  const double xp = (alpha1 * x1 + alpha2 * x2) / p;
  const double yp = (alpha1 * y1 + alpha2 * y2) / p;
  const double zp = (alpha1 * z1 + alpha2 * z2) / p;
  double S_ij, S_ij2p, S_ij2m, S_kl, S_kl2p, S_kl2m, S_mn, S_mn2p, S_mn2m;
  calcEij3(S_ij, S_ij2p, S_ij2m, nx1, nx2, p, mu, xp - x1, xp - x2);
  calcEij3(S_kl, S_kl2p, S_kl2m, ny1, ny2, p, mu, yp - y1, yp - y2);
  calcEij3(S_mn, S_mn2p, S_mn2m, nz1, nz2, p, mu, zp - z1, zp - z2);
  double T_ij = -2 * alpha2 * alpha2 * S_ij2p + alpha2 * (2 * nx2 + 1) * S_ij -
                0.5 * nx2 * (nx2 - 1) * S_ij2m;
  double T_kl = -2 * alpha2 * alpha2 * S_kl2p + alpha2 * (2 * ny2 + 1) * S_kl -
                0.5 * ny2 * (ny2 - 1) * S_kl2m;
  double T_mn = -2 * alpha2 * alpha2 * S_mn2p + alpha2 * (2 * nz2 + 1) * S_mn -
                0.5 * nz2 * (nz2 - 1) * S_mn2m;
  return (T_ij * S_kl * S_mn + S_ij * T_kl * S_mn + S_ij * S_kl * T_mn) *
         pow(acos(0.0) * 2.0 / p, 1.5);
}
double MurDavPrimIntegrals::Vij(const int &nx1, const int &ny1, const int &nz1,
                                const double &x1, const double &y1,
                                const double &z1, const double &alpha1,
                                const int &nx2, const int &ny2, const int &nz2,
                                const double &x2, const double &y2,
                                const double &z2, const double &alpha2,
                                const double &xq, const double &yq,
                                const double &zq) {
  const double p = alpha1 + alpha2;
  const double mu = alpha1 * alpha2 / p;
  const double xp = (alpha1 * x1 + alpha2 * x2) / p;
  const double yp = (alpha1 * y1 + alpha2 * y2) / p;
  const double zp = (alpha1 * z1 + alpha2 * z2) / p;

  const int tmax = nx1 + nx2 + 1;
  const int umax = ny1 + ny2 + 1;
  const int vmax = nz1 + nz2 + 1;

  double *Rntuv = new double[tmax * umax * vmax * (tmax + umax + vmax - 2)];
  double *Eij = new double[tmax];
  double *Ekl = new double[umax];
  double *Emn = new double[vmax];

  calcRntuv(Rntuv, tmax - 1, umax - 1, vmax - 1, p, xp - xq, yp - yq, zp - zq);
  calcEijt(Eij, nx1, nx2, p, mu, xp - x1, xp - x2);
  calcEijt(Ekl, ny1, ny2, p, mu, yp - y1, yp - y2);
  calcEijt(Emn, nz1, nz2, p, mu, zp - z1, zp - z2);

  double vt = 0;
  for (int t = 0; t < tmax; t++)
    for (int u = 0; u < umax; u++)
      for (int v = 0; v < vmax; v++)
        vt += Eij[t] * Ekl[u] * Emn[v] * Rntuv[t * umax * vmax + u * vmax + v];

  delete[] Eij;
  delete[] Ekl;
  delete[] Emn;
  delete[] Rntuv;
  return vt * acos(0.0) * 4.0 / p;
}

double MurDavPrimIntegrals::Vijkl(
    const int &nx1, const int &ny1, const int &nz1, const double &x1,
    const double &y1, const double &z1, const double &alpha1, const int &nx2,
    const int &ny2, const int &nz2, const double &x2, const double &y2,
    const double &z2, const double &alpha2, const int &nx3, const int &ny3,
    const int &nz3, const double &x3, const double &y3, const double &z3,
    const double &alpha3, const int &nx4, const int &ny4, const int &nz4,
    const double &x4, const double &y4, const double &z4,
    const double &alpha4) {
  const double p = alpha1 + alpha2;
  const double mup = alpha1 * alpha2 / p;
  const double xp = (alpha1 * x1 + alpha2 * x2) / p;
  const double yp = (alpha1 * y1 + alpha2 * y2) / p;
  const double zp = (alpha1 * z1 + alpha2 * z2) / p;
  const double q = alpha3 + alpha4;
  const double muq = alpha3 * alpha4 / q;
  const double xq = (alpha3 * x3 + alpha4 * x4) / q;
  const double yq = (alpha3 * y3 + alpha4 * y4) / q;
  const double zq = (alpha3 * z3 + alpha4 * z4) / q;
  const double al = p * q / (p + q);
  const int t1max = nx1 + nx2 + 1;
  const int u1max = ny1 + ny2 + 1;
  const int v1max = nz1 + nz2 + 1;
  double *E1ij = new double[t1max];
  double *E1kl = new double[u1max];
  double *E1mn = new double[v1max];
  calcEijt(E1ij, nx1, nx2, p, mup, xp - x1, xp - x2);
  calcEijt(E1kl, ny1, ny2, p, mup, yp - y1, yp - y2);
  calcEijt(E1mn, nz1, nz2, p, mup, zp - z1, zp - z2);

  const int t2max = nx3 + nx4 + 1;
  const int u2max = ny3 + ny4 + 1;
  const int v2max = nz3 + nz4 + 1;
  double *E2ij = new double[t2max];
  double *E2kl = new double[u2max];
  double *E2mn = new double[v2max];
  calcEijt(E2ij, nx3, nx4, q, muq, xq - x3, xq - x4);
  calcEijt(E2kl, ny3, ny4, q, muq, yq - y3, yq - y4);
  calcEijt(E2mn, nz3, nz4, q, muq, zq - z3, zq - z4);

  const int tmax = t1max + t2max - 1;
  const int umax = u1max + u2max - 1;
  const int vmax = v1max + v2max - 1;

  double *Rntuv = new double[tmax * umax * vmax * (tmax + umax + vmax - 2)];

  calcRntuv(Rntuv, tmax - 1, umax - 1, vmax - 1, al, xp - xq, yp - yq, zp - zq);

  double vt = 0;
  for (int t1 = 0; t1 < t1max; t1++) {
    for (int u1 = 0; u1 < u1max; u1++) {
      for (int v1 = 0; v1 < v1max; v1++) {
        for (int t2 = 0; t2 < t2max; t2++) {
          for (int u2 = 0; u2 < u2max; u2++) {
            for (int v2 = 0; v2 < v2max; v2++) {
              vt += E1ij[t1] * E1kl[u1] * E1mn[v1] * E2ij[t2] * E2kl[u2] *
                    E2mn[v2] *
                    Rntuv[(t1 + t2) * umax * vmax + (u2 + u1) * vmax +
                          (v1 + v2)] *
                    pow(-1, v2 + u2 + t2);
            }
          }
        }
      }
    }
  }

  delete[] E1ij;
  delete[] E1kl;
  delete[] E1mn;
  delete[] E2ij;
  delete[] E2kl;
  delete[] E2mn;
  delete[] Rntuv;

  return vt * 8 * acos(0.0) * acos(0.0) / p / q * sqrt(acos(0.0) * 2 / (p + q));
}
