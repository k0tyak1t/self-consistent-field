#include"standard_matrices.h"
#include"matrix.h"
#include<iostream>
standard_matrices::standard_matrices(){nAO=0;}

standard_matrices::~standard_matrices(){if (nAO > 0) {delete[] Vee;} }

standard_matrices::standard_matrices(int new_nAO){
	if (new_nAO <= 0) {
		std::cerr<<"size of nAO <=0 \n";
	}
	nAO = new_nAO;
	S.init(new_nAO);
	T.init(new_nAO);
	H.init(new_nAO);
	Ven.init(new_nAO);
	Vee = new double[new_nAO*new_nAO*new_nAO*new_nAO];
}

int standard_matrices::init(int n_new){
	if (nAO>0) {
		std::cerr<<"nAO yet exist\n";
		return 1;
	}
	if (n_new<=0) {
		std::cerr<<"size of nAO <=0 \n";
		return 1;
	}
	nAO=n_new;
	S.init(n_new);
	T.init(n_new);
	H.init(n_new);
	Ven.init(n_new); 
	Vee=new double [n_new*n_new*n_new*n_new];
	return 0;
}

int standard_matrices::set_total_Vnn(const double a) {
	total_Vnn=a;
	return 0;
}

double standard_matrices::get_total_Vnn() const {
	return total_Vnn;
}

int standard_matrices::get_nAO() const {
	return nAO;
}
	
double standard_matrices::get_Vee(int i, int j, int k, int l) const {  // (ij|kl)
	if (((i<0)||(i>=nAO)) || ((j<0)||(j>=nAO)) || ((k<0)||(k>=nAO)) || ((l<0)||(l>=nAO))) {
		std::cerr<<" **ERROR*** numeration of indexes is wrong"<<"  get_Vee "<<'\n';
		return 1;
	};
	return Vee[i + nAO*j + nAO*nAO*k + nAO*nAO*nAO*l];
}

int standard_matrices::set_Vee(int i, int j, int k, int l, double a) { // (ij|kl)
	if (((i<0)||(i>=nAO)) || ((j<0)||(j>=nAO)) || ((k<0)||(k>=nAO)) || ((l<0)||(l>=nAO))) {
		std::cerr<<" **ERROR*** numeration of indexes is wrong"<<"  set_Vee "<<'\n';
		return 1;
	};
	Vee[i + nAO*j + nAO*nAO*k + nAO*nAO*nAO*l] = a;
	return 0;
}

int standard_matrices::get_num_el() const{
	return num_el;
}

int standard_matrices::set_num_el(int new_num_el) {
	if ((new_num_el != 1) && (new_num_el % 2) !=0) {
		std::cerr<<"** ERROR *** wrong number of electrons "<<'\n';
		return 1;
	};
	num_el=new_num_el;
	return 0;
}



