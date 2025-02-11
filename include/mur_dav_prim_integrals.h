#pragma once

class MurDavPrimIntegrals {
public:
  double Sij(int, int, int, double, double, double, double, int, int, int,
             double, double, double, double);

  double Tij(int, int, int, double, double, double, double, int, int, int,
             double, double, double, double);

  double Vij(int, int, int, double, double, double, double, int, int, int,
             double, double, double, double, double, double, double);

  double Vijkl(int, int, int, double, double, double, double, int, int, int,
               double, double, double, double, int, int, int, double, double,
               double, double, int, int, int, double, double, double, double);

private:
  int calcEijt(double *, int, int, double, double, double, double);
  int calcEij3(double &, double &, double &, int, int, double, double, double,
               double);
  int calcRntuv(double *, int, int, int, double, double, double, double);
};
