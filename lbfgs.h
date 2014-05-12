/*
  lbfgs.h 2007-03-29

  Modified from CRF++ toolkit http://sourceforge.net/projects/crfpp/
  
  Author taku
*/


#ifndef LBFGS_H
#define LBFGS_H
#include <vector>
#include <stdio.h>
using namespace std;

class LBFGS{
private:
	int n;
	int m;
	int iflag;
	double beta;
	vector<double> alpha;
	vector<double> rho;
	int iter;
	bool finish;
	vector<double *> s;
	vector<double *> y;
	double *q;
	double stp1;
	double stp;
	int info;
	int bound;
	int nfev;
	void save(double *buf, char *file_name, int index=-1);
	void load(double *buf, char *file_name, int index=-1);
	int prior;//0: speed prior, 1: mem prior
	double *diag;
	double *w;
public:
	explicit LBFGS(): n(0), m(5), iflag(0) , finish(false), iter(0), prior(0), w(NULL), diag(NULL){}
	~LBFGS();
	bool init(int parameter_num, int depth, int pri);
	int optimize(double *x, double *f, double *g, int orthant, double *w0, double *w1);
};


#endif
