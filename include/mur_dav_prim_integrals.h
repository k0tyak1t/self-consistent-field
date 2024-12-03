#ifndef MURDAVPRIMINTEGRALS_H
#define MURDAVPRIMINTEGRALS_H

class MurDavPrimIntegrals {
public:
    double Sij(const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&,
        const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&);

    double Tij(const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&,
        const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&);

    double Vij(const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&,
        const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&,
        const double&,
        const double&,
        const double&);

    double Vijkl(const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&,
        const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&,
        const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&,
        const int&,
        const int&,
        const int&,
        const double&,
        const double&,
        const double&,
        const double&);

private:
    int calcEijt(double*, const int&, const int&, const double&, const double&, const double&, const double&);
    int calcEij3(double&, double&, double&, const int&, const int&, const double&, const double&, const double&, const double&);
    int calcRntuv(double*, const int&, const int&, const int&, const double&, const double&, const double&, const double&);
};

#endif // MURDAVPRIMINTEGRALS_H
