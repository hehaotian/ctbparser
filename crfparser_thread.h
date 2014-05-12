#ifndef CRFPARSER_THREAD_H
#define CRFPARSER_THREAD_H
#include "thread.h"
#include "crfparser.h"
#include <vector>
using namespace std;
struct sentence;
struct sentence_tmp;
class crfparser;

typedef struct span{
	int head;//head < tail
	int tail;
	int y;//-2: invalid, -1: no label
	double score;
	span *left;//if left == NULL, right!=NULL, span signature is 4 <-
	span *right;//if left != NULL, right==NULL, span signature is 3 ->
}span;

class span_cmp{
	public:
	bool operator()(const span &s, const span &t) const{
		if(s.y==-2)
			return false;
		if(t.y==-2)
			return true;
		return s.score>t.score;
	}
};

/*
signature=
0						F  F  (T, if |i-j|=1; F else)
1						F  T  F
2						T  F  F
3						F  T  T
4						T  F  T
5						arc(i->j)
6						arc(i<-j)
*/
class crfparser_thread: public thread{
public:
	crfparser *_c;
	unsigned int _start_i;
	void build_lattice(sentence &sen);
	void collapse_lattice(sentence &sen);
	void build_lattice_max(sentence &sen);
	double sen_fx_gx(sentence &sen);
	void inside_outside(sentence &sen, double &z);
	void cyk(sentence &sen);
	void cyk_trace(span &s);
	int cell_index(int i,int j, int l);
	//temp working space
	double _obj;
	double *_gradient;
	vector<double> _beta;
	vector<double> _alpha;
	vector<double> _lattice;
	vector<bool> _has_alpha;
	vector<bool> _has_beta;
	vector<double> _margin;
	vector<double> _edgep;
	vector<double> _par_margin;
	vector<double> _child_margin;
	vector<span> _left_opt;
	vector<span> _right_opt;
	void run();
	void join_span(int index1,int index2,int index3);
	bool merge_span_nbest(int i,int k,int j,int l);
	void split_span(int index1,int index2,int index3, int index4, int index5);// beta[index5]-> beta[index1]+alpha[index4]
	void close_span(int i,int j,int l);
	void close_span_nbest(int i,int j,int l);
	void disclose_span(int i,int j,int l);
	void calculate_gradient(sentence &sen, double &z);
	void calculate_gradient_factor(sentence &sen, double &z);
	void calculate_margin(sentence &sen, double &z);
	double tree_cost(sentence &sen);
	void ap_update(sentence &sen);
	void pa_update(sentence &sen);
	vector<vector<int> > _par;
	vector<vector<int> > _label;
	vector<int> _cur_par;
	vector<int> _cur_label;
	vector<int> _best_label;//best relation type on each edge
	vector<int> _best_alpha_label;//best relation type on each edge
	int _times;
	vector<span> _spans;
	void assign_label(sentence &s, int *par, int *label);
};

#endif