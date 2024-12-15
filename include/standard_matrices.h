#ifndef STANDART_MATRICES_H
#define STANDART_MATRICES_H
#include "matrix.h"

class standard_matrices {
public:
    standard_matrices(const unsigned int);
    standard_matrices();
    ~standard_matrices();
    void init(const unsigned int);
    const double get_Vee(int, int, int, int) const;
    void set_Vee(int, int, int, int, double);
    matrix S, T, H, Ven, X;
    const int get_num_el() const;
    void set_num_el(const unsigned int);
    void set_total_Vnn(const double);
    double get_total_Vnn() const;
    const int get_nAO() const;

private:
    int nAO;
    int num_el;
    double total_Vnn;
    double* Vee;
};

#endif // STANDARD_MATRICES_H
