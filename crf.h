#ifndef CRF_H
#define CRF_H

#include <map>
#include <vector>
#include <string>
#include <stdio.h>
#include "freelist.h"
#include "crf_thread.h"
#include "dat.h"
#include "const.h"
#include "fun.h"
using namespace std;


typedef freelist<char, str_length> charlist;



typedef struct templet
{
	vector<pair<int,int> > x;
	vector<int> y;
	bool end_of_group;
	int groupid;
}templet;
/*
define the templet "%x[1,0]%x[-1,0]y[0]" in templet file
here,

x={[1,0],[-1,0]}
y={0}

templets share the same y are clustered into one group:

%x[0,0]%y[0]			groupid=0	end_of_group=true
%x[0,0]%y[-1]%y[0]		groupid=1	end_of_group=false
x[-2,0]%y[-1]%y[0]		groupid=1	end_of_group=true
*/


class templet_cmp
{
public:
	bool operator()(const templet &s, const templet &t) const
	{
		for(int i=0;i<s.y.size();i++)
		{
			if(s.y[i]<t.y[i])
				return false;
			else if(s.y[i]>t.y[i])
				return true;
		}
		return false;
	}
};


struct node;

typedef struct clique{
	int *fvector;
	double *fvalue;
	int feature_num;
	node **nodes;
	int node_num;
	int groupid;
	int key;
}clique;
/*
for a concrete x, CRF use all templets for matching, and creates one clique for each templet group
fvector[i] ~ fvector[i] + ysize^node_num is the lambda space for the i th templet in the group
(ysize is the tag number, for "BIO" tagging, ysize=3)
feature_num is the group size.
nodes is the tokens linked by the clique, e.g. for "y[-1]y[0]" clique, the previous and current tokens are linked
node_num is the number of nodes, e.g. it is 2 for "y[-1]y[0]" clique
key: if the nodes are tagged, then clique.key is the index of the corresponding lambda in the clique's lambda space
e.g. for "y[-1]y[0]" clique, nodes are tagged as "IB", then
key = 2*3^1+0*3^0 = 6, i.e. r*f_i = lambda[fvector[i]+6] (f_i=1 here)
*/

typedef struct node{
public:
	clique **cliques;
	int clique_num;
	int key;
}node;

/*
cliques are the cliques linked with current node
*/

typedef struct sequence{
	node* nodes;
	int node_num;
}sequence;

/*
node is the structure of token, path is the edge.
For "BIO" tagging, a sequence "x1...xn" contains n nodes,
if the order of CRF is 1, i.e. uses bigram features at most, then
the path number is n * 3 * 3

B	  .			.


I	  .			.


O	  .			.

	node_i	node_{i+1}

  corresponding paths: (i_B,{i+1}_B),(i_B,{i+1}_I),...(i_O,{i+1}_O)
*/

/*
following structures used for 1st order markov chain only
*/
typedef struct vertex{
	int key;
	int feature_num;
	int *fvector;
	double *fvalue;
}vertex;

typedef struct edge{
	int feature_num;
	int *fvector;
	double *fvalue;
}edge;

typedef struct sequence1{
	int vertex_num;
	vertex *vertexes;
	edge *edges;
}sequence1;


class crf_thread;
class CRF{
friend class crf_thread;
public:
	CRF();
	bool set_para(char *para_name, char *para_value);
	~CRF();
	bool learn(char* templet_file, char *training_file, char *model_file);
	bool load_model(char *model_file);
	void tag(vector<vector<vector<string> > > &ext_table, vector<vector<vector<double> > > &val_table, vector<vector<string> > &best_tag,vector<double> &sequencep, vector<vector<double> > &nodep);
	void tag(vector<vector<vector<string> > > &ext_table, vector<vector<vector<double> > > &val_table, vector<vector<string> > &best_tag,vector<double> &sequencep, vector<vector<double> > &nodep, vector<int> &con_pos, vector<string> &con_tag);
	//see readme.html
	bool is_real;
	vector<char *>tags;//sorted in alphabet order
private:
	bool load_templet(char *templet_file);
	//load templets from file.
	bool add_templet(char *line);
	//add templets, line is the use defined templet string, cur_group is the number of group currently formed.
	bool set_order();
	//set groupid,end_of_group,order
	void set_group();
	//after templets loading, set templet_group
	void set_chain_type();
	//set chain_type
	bool check_training(char *training_file);
	//check the training file , get x cols, tag number
	bool generate_features(char *training_file);
	//generate features according to templates and training file
	void compress();
	void adjust_data();
	//adjust fmap, lambda_size
	vector<templet> templets;
	//templates
	vector<vector<vector<int> > > templet_group;
	//templet_group[i] is the i th group
	//templet_group[i].size() is the number of parameters number of i th group (size of the parameter space)
	//e.g. for "y[-1]y[0]" group, its size is 3*3=9
	//templet_group[i][j] is the indexes of the paths that affect the i th group, j th parameter
	//for 2 order CRF, ("y[-2]y[-1]y[0]"), group[i] is "y[-2]y[0]" group, tagged with IO, 
	//group[i][1*3^1+2*3^0]=group[i][5]
	//then paths that affected the 5th parameter are the {5,14,23} th paths
	//group[i][5]={5,14,23}
	int gsize;//templet_group.size()
	//x
	map<char *, int, str_cmp> xindex;//<"1:Confidence", 132> > 132 th x
	vector<int> x_freq;//x_freq[i] is the frequency of the i th x
	charlist x_str;//space storing x

	
	charlist tag_str;//space storing tags
	
	double *lambda;//parameters
	size_t lambda_size;//parameter size

	int cols;//cols of training file
	int ysize;//tags.size()
	int order;//order of CRF

	double sigma;//sigma||\lambda||^2
	int freq_thresh;//see readme.html
	int max_iter;//max iteration number
	double eta;//controls the iteration accuracy

	double* gradient;//gradient=d L(\lambda) / d \lambda
	double fvalue;//L(\lambda)
	//L(\lambda)=\sum_{sequences} [ z(\lambda) - \sum_{f \in y} \lambda * f(x,y) ]
	//z(\lambda)=\log{\sum_y exp(\sum_{f \in y} \lambda * f(x,y)) }


	vector<sequence> sequences_tmp;
	sequence *sequences;
	sequence1 *sequence1s;
	int sequence_num;
	freelist<node> nodes;//storing all nodes
	freelist<clique> cliques;//storing all cliques
	freelist<node*> clique_node;//storing clique linked nodes' addresses
	freelist<clique*> node_clique;//storing nodes linked cliques' addresses
	freelist<int> clique_feature;//clique->feature
	freelist<double> clique_fvalue;//clique->fvalue
	bool add_x(vector<char *> &table);//add x from table to xindex
	bool insert_x(char *x, int & index);//insert x to xindex, return its index
	void shrink_feature();//eliminate the x with frequency less than freq_thresh
	void write_model(char *model_file,bool first_part);
	//test code
	void generate_sequence(vector<vector<vector<string> > >&ext_table, vector<vector<vector<double> > >&val_table, sequence &seq);
	//generate seq according to user's input table
	bool margin;
	//true if margin probabilities for each node calculation is required
	bool seqp;
	int nbest;
	//output nbest results
	//global temporary variables
	int path_num;//=pow((double)ysize,order+1);
	int node_anum;//pow((double)ysize,order),alpha(beta) number of each node
	double head_offset;//=log((double)ysize)*order; node_{-1}'s prior possibility


	vector<int> fmap_tmp;
	int *fmap;
	int fmap_size;
	int version;

	vector<crf_thread> threads;
	unsigned int thread_num;

	char *work_space;
	unsigned long work_size;
	void unload();
	void load();

	int chain_type;
	int transit;
	int max_seq_size;
	double *transit_buff;
	int algorithm;
	int total_nodes;
	int depth;
	int prior;
	double_array_trie *dat;
	vector<string> table_head;
	vector<string> table_tail;
};
#endif
