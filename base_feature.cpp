#include "base_feature.h"
#include "fun.h"
#include <iostream>
#include <fstream>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#include "const.h"

using namespace std;

vector<string> base_feature::_labels;

feature_pool::feature_pool(){
	_feature_num=0;
	_cols=0;
	_sen_num=0;
	_total_words=0;
	_total_bytes=0;
	_max_sen_bytes=0;
}

feature_pool::~feature_pool(){
	int i;
	for(i=0;i<_features.size();i++){
		delete _features[i];
	}
}

bool feature_pool::check_training(char *training_file){
	ifstream fin;
	int i;
	fin.open(training_file);
	if(!fin.is_open()){
		cout<<"can not open file: "<<training_file<<endl;
		return false;
	}
	
	int lines=0;
	_cols=0;
	vector<vector<string> > table;
	bool has_next;
	
	do{
		has_next=build_table(fin,table,lines);
		if(table.size()){
			for(i=0;i<table.size();i++){
				if(_cols && _cols!=table[i].size()){//incompatible
					lines=lines-table.size()+i;
					if(!has_next)
						lines++;
					cout<<"line "<<lines<<": columns incompatible"<<endl;
					fin.close();
					return false;
				}
				_cols=table[i].size();
				if(_cols<3){
					lines=lines-table.size()+i;
					if(!has_next)
						lines++;
					cout<<"line "<<lines<<": columns should be no less than 3"<<endl;
					fin.close();
					return false;
				}
				string t=table[i].back();//label
				int par_pos=atoi(table[i][table[i].size()-2].c_str());
				if(par_pos!=table.size()){
					int index,insert_pos;
					if(!vector_search(base_feature::_labels,t,index,insert_pos)){
						vector_insert(base_feature::_labels,t,insert_pos);
					}
				}
			}
			//check projective
			if(!is_projective(table)){
				lines=lines-table.size();
				if(!has_next)
					lines++;
				cout<<"sentence beginning at line "<<lines<<" is not projective"<<endl;
				fin.close();
				return false;
			}
			int rn=root_num(table);
			if(rn>1){
				lines=lines-table.size();
				if(!has_next)
					lines++;
				cout<<"sentence beginning at line "<<lines<<" has multiple roots"<<endl;
				fin.close();
				return false;
			}else if(rn==0){
				lines=lines-table.size();
				if(!has_next)
					lines++;
				cout<<"sentence beginning at line "<<lines<<" has no roots"<<endl;
				fin.close();
				return false;
			}
		}
	}while(has_next);
	fin.close();
	return true;
}


int base_feature::get_label_index(string &t){
	int index,insert_pos;
	if(!vector_search(base_feature::_labels,t,index,insert_pos))
		return -1;//not found
	return index;
}

void feature_pool::push_back(base_feature *bf){
	bf->_base_feature_id=_feature_num;
	_feature_num=_feature_num+bf->_feature_num;
	_features.push_back(bf);
	_base_feature_ids.push_back(bf->_base_feature_id);
}



void feature_pool::generate_feature(char *training_file){
	int lines=0;
	ifstream fin;
	vector<vector<string> > table;
	bool has_next;
	int i=0;
	cout<<"generate sentences"<<endl;
	fin.open(training_file);
	do{
		has_next=build_table(fin,table,lines);
		if(table.size()){
			sentence_tmp st;
			generate_sentence(table,st);
			char fn[100];
			sprintf(fn,"__sentence%d",i);
			save_sentence(st,fn);
			i++;
			cout<<i<<" ";
		}
	}while(has_next);
	cout<<endl;
	fin.close();
}

void feature_pool::pop(){
	for(int i=0;i<_features.size();i++){
		base_feature *bf=_features[i];
		//_feature_num=_feature_num-bf->_feature_num;
		delete bf;
	}
	_features.clear();
}

void feature_pool::generate_feature(){
	long cur_bytes;
	sentence sen;
	for(int i=0;1;i++){
		char fn[100];
		sprintf(fn,"__sentence%d",i);
		if(!access(fn,0)){//has file
			_sen_num++;
			FILE *fp;
			fp=fopen(fn,"rb");
			fread(&cur_bytes,sizeof(long),1,fp);
			if(_max_sen_bytes<cur_bytes)
				_max_sen_bytes=cur_bytes;
			fread(&sen,sizeof(sentence),1,fp);
			_total_words+=sen.l;
			fclose(fp);
		}else{
			break;
		}
	}
}
void feature_pool::generate_sentence(vector<vector<string> > & table, sentence_tmp & st){
	int i,j,k;
	int rows=table.size();
	_sen_num++;
	_total_words+=rows;

	st.label.clear();
	st.parent.clear();
	st.fvector.clear();

	st.label.resize(rows);
	st.parent.resize(rows);
	st.fvector.resize(rows*2+rows*rows);//w_n act as child: len edges
	for(i=0;i<table.size();i++){
		int index;
		int par=atoi(table[i][table[i].size()-2].c_str());
		if(par<table.size())
			index=base_feature::get_label_index(table[i].back());
		else
			index=-1;
		st.label[i]=index;
		st.parent[i]=par;
	}
	for(i=0;i<_features.size();i++){
		vector<vector<feature> > vf(rows*2+rows*rows);
		_features[i]->generate_feature(table,vf);
		//fill to sentence
		for(j=0;j<vf.size();j++){
			for(k=0;k<vf[j].size();k++){
				feature f;
				f.index=_features[i]->_base_feature_id+vf[j][k].index;
				f.y=vf[j][k].y;
				st.fvector[j].push_back(f);
			}
		}
	}
}


void feature_pool::save_sentence(sentence_tmp & st,char *fn){
	//memory format:
	int i,j,k,ii;
	int ysize=base_feature::_labels.size();
	
	

	long cur_bytes=0;
	cur_bytes+=sizeof(sentence);
	int l=st.label.size();
	cur_bytes+=l*sizeof(int)*2;//parent, label
	cur_bytes+=l*(l+2)*sizeof(int *);//fnum
	cur_bytes+=l*(l+2)*sizeof(int)*ysize;//fnum
	cur_bytes+=l*(l+2)*sizeof(int **);//fvector
	cur_bytes+=l*(l+2)*sizeof(int *)*ysize;
	int total_fnum=0;
	for(j=0;j<st.fvector.size();j++){
		total_fnum+=st.fvector[j].size();
	}
	cur_bytes+=sizeof(int)*total_fnum;
	

	FILE *fp;
	fp=fopen(fn,"wb");
	//memory format: bytes,sentence,parents,labels,
	//	fnum[0],...,fnum[N],fnum[0][0],...fnum[N][y],
	//	fvector[0],...,fvector[N],fvector[0][0],...fvector[N][y],fvector[0][0][0],...
	//put cur_bytes
	fwrite(&cur_bytes,sizeof(long),1,fp);
	//put sentence
	sentence sen;
	int len=st.label.size();
	int sz=st.fvector.size();
	sen.l=len;
	fwrite(&sen,sizeof(sentence),1,fp);
	//put parent,label
	fwrite(&st.parent[0],sizeof(int),len,fp);
	fwrite(&st.label[0],sizeof(int),len,fp);
	//calculate fnum,fvector
	vector<vector<int> > fnum(sz);
	vector<vector<vector<int> > >fvector(sz);
	for(i=0;i<sz;i++){
		fnum[i].resize(ysize,0);
		fvector[i].resize(ysize);
		for(j=0;j<st.fvector[i].size();j++){
			int y=st.fvector[i][j].y;
			int index=st.fvector[i][j].index;
			if(y>=0){
				//cout<<i<<"\t"<<y<<endl;
				fnum[i][y]++;
				fvector[i][y].push_back(index);
			}else{//y=0, root
				fnum[i][0]++;
				fvector[i][0].push_back(index);
			}
		}
	}
	//put fnum
	fwrite(&fnum[0],sizeof(int*),sz,fp);
	for(i=0;i<sz;i++)
		fwrite(&fnum[i][0],sizeof(int),ysize,fp);
	//put fvector
	fwrite(&fvector[0],sizeof(int**),sz,fp);//&len is no use
	for(i=0;i<sz;i++)
		fwrite(&fvector[i][0],sizeof(int*),ysize,fp);//len is no use
	for(i=0;i<sz;i++)
		for(j=0;j<ysize;j++)
			if(fnum[i][j])
				fwrite(&fvector[i][j][0],sizeof(int),fnum[i][j],fp);//len is no use
	
	_total_bytes+=cur_bytes;
	if(_max_sen_bytes<cur_bytes)
		_max_sen_bytes=cur_bytes;
	fclose(fp);
}


bool feature_pool::load_sentence(char *cur_space,char *&new_space, char *fn,long &rest_bytes){
	//memory format:
	int i,j,k,ii;
	int ysize=base_feature::_labels.size();
	FILE *fp;
	fp=fopen(fn,"rb");
	long cur_bytes;
	fread(&cur_bytes,sizeof(long),1,fp);
	if(cur_bytes>rest_bytes){
		fclose(fp);
		return false;
	}
	fread(cur_space,sizeof(char),cur_bytes,fp);
	fclose(fp);
	new_space+=cur_bytes;
	rest_bytes-=cur_bytes;
	//calculate sen.parent, sen.label etc
	//memory format: bytes,sentence,parents,labels,
	//	fnum[0],...,fnum[N],fnum[0][0],...fnum[N][y],
	//	fvector[0],...,fvector[N],fvector[0][0],...fvector[N][y],fvector[0][0][0],...
	//put cur_bytes
	sentence *sen=(sentence *)cur_space;
	int l=sen->l;
	char *tmp_space=(char *)sen+sizeof(sentence);
	sen->parent=(int*)tmp_space;
	tmp_space+=sizeof(int)*l;
	sen->label=(int*)tmp_space;
	tmp_space+=sizeof(int)*l;
	sen->fnum=(int**)tmp_space;
	int sz=l*(l+2);
	tmp_space+=sizeof(int *)*sz;
	for(i=0;i<sz;i++){
		sen->fnum[i]=(int*)tmp_space;
		tmp_space+=sizeof(int)*ysize;
	}
	sen->fvector=(int***)tmp_space;
	tmp_space+=sizeof(int **)*sz;
	for(i=0;i<sz;i++){
		sen->fvector[i]=(int**)tmp_space;
		tmp_space+=sizeof(int*)*ysize;
	}
	for(i=0;i<sz;i++){
		for(j=0;j<ysize;j++){
			sen->fvector[i][j]=(int*)tmp_space;
			tmp_space+=sen->fnum[i][j]*sizeof(int);
		}
	}
	return true;
}


void feature_pool::write_sentence(sentence_tmp & st,char *fn){
	int i,j,k,ii;
	int ysize=base_feature::_labels.size();
	vector<int> fnum(ysize);
	vector<vector<int> > fvector(ysize);
	ofstream fout;
	fout.open(fn);
	fout<<st.label.size()<<endl;
	for(j=0;j<st.label.size();j++)
		fout<<st.label[j]<<endl;
	for(j=0;j<st.parent.size();j++)
		fout<<st.parent[j]<<endl;
	fout<<st.fvector.size()<<endl;
	for(j=0;j<st.fvector.size();j++){
		fill(fnum.begin(),fnum.end(),0);
		fvector.clear();
		fvector.resize(ysize);
		for(k=0;k<st.fvector[j].size();k++){
			int y=st.fvector[j][k].y;
			if(y>=0){
				fnum[y]++;
				fvector[y].push_back(st.fvector[j][k].index);
			}else{//y=1, root
				fnum[0]++;
				fvector[0].push_back(st.fvector[j][k].index);
			}
		}
		for(k=0;k<ysize;k++)
			fout<<fnum[k]<<endl;
		for(k=0;k<ysize;k++){
			for(ii=0;ii<fnum[k];ii++)
				fout<<fvector[k][ii]<<endl;
		}
	}
	long cur_bytes=0;

	cur_bytes+=sizeof(sentence);
	
	int l=st.label.size();
	cur_bytes+=l*sizeof(int)*2;//parent, label
	cur_bytes+=l*(l+2)*sizeof(int *);//fnum
	cur_bytes+=l*(l+2)*sizeof(int)*ysize;//fnum
	cur_bytes+=l*(l+2)*sizeof(int **);//fvector
	cur_bytes+=l*(l+2)*sizeof(int *)*ysize;
	int total_fnum=0;
	for(j=0;j<st.fvector.size();j++){
		total_fnum+=st.fvector[j].size();
	}
	cur_bytes+=sizeof(int)*total_fnum;
	_total_bytes+=cur_bytes;
	if(_max_sen_bytes<cur_bytes)
		_max_sen_bytes=cur_bytes;
	fout.close();
}


//is projective
bool feature_pool::is_projective(vector<vector<string> > &table){
	vector<int> pos1,pos2;
	int i,j;
	//pos1[i] is the former position of edge i
	//pos2[i] is the latter position of edge i
	//pos1.size()=pos2.size(), pos1[i]<pos2[i]
	//complexity: n^2, n=pos1.size()
	for(i=0;i<table.size();i++){
		j=atoi(table[i][table[i].size()-2].c_str());
		if(i>j){
			pos1.push_back(j);
			pos2.push_back(i);
		}else{
			pos2.push_back(j);
			pos1.push_back(i);
		}
	}
	for(i=1;i<pos1.size();i++){
		for(j=0;j<i;j++){
			if(pos1[i]==pos1[j] || pos1[i]==pos2[j] || pos2[i]==pos1[j] || pos2[i]==pos2[j])
				continue;
			int rel_pos1,rel_pos2;
			if(pos1[i]<pos1[j] || pos1[i]>pos2[j])
				rel_pos1=1;//outside
			else
				rel_pos1=0;//inside
			if(pos2[i]<pos1[j] || pos2[i]>pos2[j])
				rel_pos2=1;//outside
			else
				rel_pos2=0;//inside
			if(rel_pos1!=rel_pos2)//cross
				return false;
		}
	}
	return true;
}

int feature_pool::root_num(vector<vector<string> > &table){
	int rows=table.size();
	int n=0;
	int i,j;
	for(i=0;i<rows;i++){
		j=atoi(table[i][table[i].size()-2].c_str());
		if(j==rows)
			n++;
	}
	return n;
}
