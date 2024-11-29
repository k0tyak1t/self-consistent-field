#ifndef STANDART_MATRICES_H
#define STANDART_MATRICES_H
#include"matrix.h"

class standard_matrices{
public:
	~standard_matrices();
	standard_matrices();
	standard_matrices(int);
	int init(int);
	double get_Vee(int,int,int,int) const;
	int set_Vee(int,int,int,int,double);
	matrix S, T, H, Ven, X;
	int get_num_el() const;
	int set_num_el(int);
	int set_total_Vnn(const double);
	double get_total_Vnn() const;
	int get_nAO() const;
private:
	int nAO;
	int num_el;
	double total_Vnn;
	double* Vee;
};

#endif
