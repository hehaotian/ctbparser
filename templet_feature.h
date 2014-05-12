/*
This is an implementation of 2D Trie Structure for fast template based feature extraction
Reference:
Xian Qian, Qi Zhang, Xuanjing Huang and Lide Wu,
2D Trie for fast parsing,
Proceedings of COLING 2010
*/
#ifndef TEMPLET_FEATURE_H
#define TEMPLET_FEATURE_H
#include <string>
#include <vector>
#include <map>
#include "freelist.h"
#include "fun.h"
#include "base_feature.h"
#include "dat.h"
using namespace std;

typedef struct graph_node{
	int row;
	int col;
	int type;
	//type: 0 p, 1 b, 2 c, 3 r, 4 dir, 5 dis, 6 type
	int par;
	int sib;		//number of its brothers
	int id;
	int cluster_id;
	vector<int> child;	//index of children in the tree, child.size()=pos.size()
	vector<int> follow;	//index of children in the graph, follow.size()>=child.size()
}graph_node;

typedef struct trie_node_tmp{
	int graph_id;//index of the node in graph
	vector<int> index;
	vector<vector<trie_node_tmp *> > child;	//child.size()=_graph_nodes[graph_id].pos.size()
				//child[i].size()=_code[_graph_nodes[graph_id].col_id].size();
				//child[i][j] is a trie_node*
				//if has no child, child=NULL;
}trie_node_tmp;

typedef struct trie_node{
	int graph_id;//index of the node in graph
	int *index;
	int length;
	trie_node ***child;	//child.size()=_graph_nodes[graph_id].pos.size()
				//child[i].size()=_code[_graph_nodes[graph_id].col_id].size();
				//child[i][j] is a trie_node*
				//if has no child, child=NULL;
}trie_node;

typedef int subscript;
typedef struct tunit{
	subscript *base;
	subscript check;
}tunit;


typedef struct scan_unit{
	int child_size;	//-1 means feature
	int cluster_id;	//
	int child_id;	//index of its child
	int gid;
	bool par_child;
	bool par_only;
	bool child_only;
	bool has_root;
}scan_unit;



class templet_feature: public base_feature{
public:
	templet_feature(int x_freq, int ysize);//cols of train file
	~templet_feature();
	bool load_templet(char *templet_file);
	bool construct_index(char *train_file);
	void generate_feature_candidate(char *training_file);
	void generate_feature(vector<vector<string> > & table, vector<vector<feature> > &vf);
	bool write_model(char *model_file);
	bool load_model(char *model_file,bool require_features,bool require_trie);
	void refine_feature(bool *fmap);
private:
	int _ysize;
	int _x_freq;
	map<char *, int, str_cmp> _feature2id;
	freelist<char> _feature_str;
	vector<string> _templet_str;
	vector<bool> _has_par;
	vector<bool> _has_child;
	vector<bool> _has_root;
	vector<bool> _has_between;
	vector<bool> _has_distance;
	vector<bool> _has_direction;
	vector<bool> _has_label;
	vector<bool> _par_only;
	vector<bool> _child_only;
	vector<bool> _par_child;
	vector<int> _ysz;

	
	bool insert_feature(string &fs);
	bool search_feature(string &fs, int &index);
	bool get_feature_string(vector<vector<string> > & table, int templet_id, int par_pos,int child_pos,int y, vector<string> &s);
	void construct_templet_graph(vector<string> &templets);
	trie_node *_root;
	int get_col_id(char *,int col);
	vector<vector<int> > _templet_path;
	vector<graph_node> _gnode;
	int _tsize;

	
	freelist<char> _data;

	int get_csize(graph_node &gn);
	bool add_feature(vector<vector<int> > &ids, int par_pos,int child_pos, vector<vector<feature> > &vf);

	bool add_feature(vector<vector<string> > &table);

	vector<vector<int> > _table_head;//_table_head[i][2]=index of "B_2"
	vector<vector<int> > _table_tail;//_table_tail[i][2]="E_2"

	void table2id(vector<vector<string> > & table, vector<vector<int> > & index_table);
	bool insert_trie(vector<int> &path, vector<int> &ids, trie_node_tmp *root, int y, int ysz, int index);




	
	vector<int> _scan_order;
	vector<graph_node> _gcluster;

	void generate_ids(vector<vector<int> > &index_table, vector<vector<int> > &ids);

	tunit *_dat;
	bool *_used;
	vector<subscript> _bas;
	subscript _dat_size;
	subscript dat_resize(subscript);
	void insert_dat(trie_node_tmp *s, vector<subscript> &bases);
	subscript _dat_next_pos;
	vector<subscript> _bases;
	vector<int> _child_size;
	vector<int> _cluster_id;
	vector<int> _que;
	vector<scan_unit> _su;//
	vector<double_array_trie> _da;
};







#endif
