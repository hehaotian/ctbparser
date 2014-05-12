#ifndef CRF_PARSER_H
#define CRF_PARSER_H
#include "base_feature.h"
#include "crfparser_thread.h"
#include "freelist.h"

#include <vector>
using namespace std;
typedef struct feature{
	int y;//dependency label
	int index;//feature index
}feature;
typedef struct sentence{
	int l;//length
	int *parent;//parent[i] is the parent of i
	int *label;//
	int **fnum;//fnum[i] feature number in grid[i], (l+2)*l grids
	//i*l+j i!=j, i->j
	//      i==j, i<-root
	//l*l+i i act as parent
	//l*l+l+i act as child
	int ***fvector;
}sentence;
class feature_pool;
class crfparser_thread;
class switch_thread;


class crfparser{
friend class crfparser_thread;
friend class switch_thread;
public:
	bool learn_api(char *training_file);
	bool model_api(bool **fmap,char *training_file, int action);//action=0: refine, 1:generate, 2:load

	crfparser();
	~crfparser();
	bool learn(char *training_file,char *model_file);//return whether learning success
	void tag(vector<vector<string> > &table, vector<vector<int> > &parent, vector<vector<string> > &label,vector<double> &treep, vector<vector<vector<double> > > &edgep);
	
	bool load_model(char *model_file);
	bool set_para(char *para_name,char *para_value);
private:
	char *_work_space;
	long _work_bytes;
	long _max_sen_bytes;
	char *_total_space;
	vector<sentence> _sentences;
	double *_lambda;
	double *_gradient;
	feature_pool *_fpool;

	long _total_bytes;
	freelist<char> _temp_space;


	int _sen_num;
	void copy_sentence();
	int _lambda_size;
	int _thread_num;
	vector<crfparser_thread> _threads;
	int _max_iter;
	double _sigma;
	double _eta;
	bool write_model(char *model_file, int part=0);
	void convert2sentence(sentence_tmp &st, sentence &sen);
	bool load_sentence();
	int _version;
	int _depth;

	void unload();
	void load();
	int _margin;
	
	int _prior;
	int _algorithm;
	int _training_model;
	int _nbest;
	
	
	bool _need_scale;
	double _factor_gradient;
	double _factor;
	int _start;
};
#endif
