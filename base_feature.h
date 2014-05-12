#ifndef BASE_FEATURE_H
#define BASE_FEATURE_H

#include <vector>
#include <string>
#include "crfparser.h"
using namespace std;
class base_feature;

struct feature;
struct sentence;

typedef struct sentence_tmp{
	vector<int> parent;
	vector<int> label;
	vector<vector<feature> > fvector;
}sentence_tmp;

class feature_pool{
public:
	friend class crfparser;
	friend class base_feature;
	bool check_training(char *training_file);
	void push_back(base_feature *bf);
	void pop();//pop and clear all features
	feature_pool();
	~feature_pool();
	void generate_feature(char *training_file);
	
private:
	int _feature_num;
	int _cols;//column number of training file.
	int _sen_num;
	
	void generate_sentence(vector<vector<string> >&table, sentence_tmp & st);
	void write_sentence(sentence_tmp & st,char *filename);
	void save_sentence(sentence_tmp & st,char *filename);
	bool load_sentence(char *cur_space, char *&new_space, char *filename,long &rest_bytes);
	bool is_projective(vector<vector<string> > &table);
	void generate_feature();
	int root_num(vector<vector<string> > &table);
	vector<int> _base_feature_ids;
	vector<base_feature *> _features;
	//working varables
	vector<sentence_tmp> _sentence_tmps;
	int _total_words;
	long _total_bytes;
	long _max_sen_bytes;


};
//each feature class should be independent, ie, contains its own member
class base_feature{
public:
	friend class feature_pool;
	friend class crfparser;
	base_feature(){_feature_num=0;_base_feature_id=0;};
	virtual ~base_feature(){};
	virtual void generate_feature(vector<vector<string> > & table, vector<vector<feature> > &vf) = 0;
	virtual void refine_feature(bool *fmap){};
	static int get_label_index(string &label);//shared in train and test
	static vector<string> _labels;
protected:
	int _base_feature_id;//all features generate by templet arranges from _base_feature_id.
	int _feature_num;
};
#endif
