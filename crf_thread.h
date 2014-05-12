#ifndef CRF_THREAD_H
#define CRF_THREAD_H
#include "thread.h"
#include "crf.h"
class CRF;
struct sequence;
struct sequence1;
class crf_thread: public thread	//all data are borrowed from crf
{
public:
	CRF *c;//provide data
	unsigned int start_i;
	void build_lattice(sequence &seq);//calculate all path, i.e. \sum \lambda * f for each path
	void build_lattice(sequence1 &seq1);
	double path_cost(sequence &seq);//calculate \sum_{f \in y} \lambda * f(x,y)
	double path_cost(sequence1 &seq1);
	double seq_fx_gx(sequence &seq);//calculate z(\lambda) - \sum_{f \in y} \lambda * f(x,y) of current seq
	double seq_fx_gx(sequence1 &seq1);
	void forward_backward(sequence &seq,double &z);
	void forward_backward(sequence1 &seq1,double &z);
	//calculate alpha beta z(\lambda)
	void calculate_gradient(sequence &seq,double &z);
	void calculate_gradient(sequence1 &seq1,double &z);
	//calculate d L(\lambda) / d \lambda for current seq
	void node_margin(sequence &seq, vector<vector<double> >&node_p,double &z);
	void assign_tag(sequence &seq, vector<int> &node_tag);
	void viterbi(sequence &seq);
	void viterbi(sequence &seq, vector<int> &con_pos, vector<int> &con_y);
	void viterbi(sequence1 &seq1);
	void ap_update(sequence &seq);
	void ap_update(sequence1 &seq1);
	void pa_update(sequence &seq);
	void pa_update(sequence1 &seq1);
	//temp working space
	double obj;
	double *gradient;
	vector<double> alpha;
	vector<double> beta;
	vector<double> path;
	vector<int> first_cal;
	vector<double> margin;
	int *fmap;
	int times;//for ap
	void run();

	vector<vector<int> > best_path;
	vector<int> bst_path;
	vector<vector<double> >last_best;
	vector<vector<double> >cur_best;
	vector<vector<vector< pair<int, int> > > >best_prev;
	vector<vector<vector< int > > >best_link;
	vector<double> final_path;
	double (crf_thread::*_cal_cost)(double w, double *fvalue, int k);
	int (crf_thread::*_findex)(int index);
	inline double cal_cost(double w, double *fvalue, int k)
	{
		return (this->*_cal_cost)(w,fvalue,k);
	}
	inline int findex(int index)
	{
		return (this->*_findex)(index);
	}
	inline double bool_cost(double w, double *fvalue, int k)
	{
		return w;
	}

	inline double real_cost(double w, double *fvalue, int k)
	{
		return w*fvalue[k];
	}
	inline int map_findex(int index)
	{
		return fmap[index];
	}
	inline int dir_findex(int index){
		return index;
	}

};
#endif