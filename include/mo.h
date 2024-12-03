#ifndef MOS_H
#define MOS_H
#include "matrix.h"

class MO {
public:
    MO();
    ~MO();
    matrix C; // MO to AO expansion. MO in columns
    int init(const int);
    double get_mo_energy(const int) const;
    int get_size();
    int set_mo_energy(const int, const double);
    int set_total_energy(const double);
    double get_total_energy() const;
    int get_irrep_characters(const int&) const;
    bool set_c2v(int* symmAO, const double& limit);
    void set_mo_energies(const double*);

private:
    int n; // basis dimentionality
    double* mo_energies;
    double total_energy;
    int* irrep;
};

#endif
