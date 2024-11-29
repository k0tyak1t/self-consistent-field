#ifndef MOS_H
#define MOS_H
#include "matrix.h"
class MOs{
public:
	MOs();
	~MOs();
	matrix C; //матрица разложения МО по АО, МО в столбцах
	int init(const int);
	double get_mo_energy(const int) const;
	int get_size();
	int set_mo_energy(const int, const double);
	int set_total_energy(const double);
	double get_total_energy() const;
	double check_ij(const int, const int, const char*) const;
	int get_irrep(const int&) const;
	bool set_c2v(int* symmAO, const double& limit);
private:
	int n;   //размерность базиса данной задачи
	double* mo_energies; //массив из энергий молекулярных орбиталей
	double total_energy; //энергия всей системы
	int* irrep; // массив из неприводимых представлений
};
#endif
