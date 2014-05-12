#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <set>
#include <map>
#include <assert.h>
#include <queue>
#include <ctime>
#include "const.h"
#include "templet_feature.h"
#include "fun.h"
const char *CUT="\001";


using namespace std;





class sscmp{//set string compare
public:
	bool operator()(const set<string>& s1, const set<string>& s2)const{
		if (s1.size()<s2.size())
			return true;
		else if (s1.size()>s2.size())
			return false;
		else if (s1<s2) 
			return true;
		else return false;
	}
};



templet_feature::templet_feature(int x_freq, int ysize){
	_x_freq=x_freq;
	_ysize=ysize;
	_tsize=0;
	_data.set_size(1024*1024*256);//256M
	_feature_str.set_size(PAGESIZE);//
	trie_node s;
	s.child=0;
	s.graph_id=-1;
	_root=(trie_node *)_data.push_back((char *)&s,sizeof(trie_node));
	_dat=0;
	_used=0;
	_dat_size=0;
}


bool templet_feature::load_templet(char *templet_file){
	ifstream fin;
	//read template
	fin.open(templet_file);
	if (!fin.is_open()){
		cout<<"template file: "<<templet_file<<" not found"<<endl;
		return false;
	}
	char line[MAXSTRLEN];
	_templet_str.clear();
	while(!fin.eof()){
		fin.getline(line,MAXSTRLEN-1);
		if(line[0] && line[0]=='%')
			_templet_str.push_back(line);
	}
	fin.close();
	_tsize=_templet_str.size();
	construct_templet_graph(_templet_str);
	return true;
}



void templet_feature::construct_templet_graph(vector<string> &templets){
	char s[1024];
	set<set<string> ,sscmp> tset;
	int i,j;
	_tsize=templets.size();
	for(i=0;i<_tsize;i++){
		strcpy(s,templets[i].c_str());
		vector<char *> row;
		set<string> t;
		split_string(s,"%",row);
		for(j=0;j<row.size();j++){
			if(strlen(row[j]))
				t.insert(row[j]);
		}
		tset.insert(t);
	}
	vector<set<string> > tvec(_tsize);

	//update templet_str
	i=0;
	for(set<set<string> ,sscmp>::iterator it=tset.begin();it!=tset.end();it++,i++){
		tvec[i]=*it;
		_templet_str[i]="";
		for(set<string>::iterator it1=it->begin();it1!=it->end();it1++){
			_templet_str[i]+="%";
			_templet_str[i]+=*it1;
		}
	}
	//templets loaded
	//get parents and children
	vector<set<int> > parent(_tsize),child(_tsize);
	for (i=0;i<_tsize;i++){
		for (j=i+1;j<_tsize;j++){
			if (is_subset(tvec[i],tvec[j])){
				parent[j].insert(i);
				child[i].insert(j);
			}
		}
	}
	//remove intersection
	vector<int> dirpar(_tsize);//direct parent
	vector<set<int> > dirch(_tsize);//direct child
	vector<set<int> > dirpars(_tsize);//temp variable
	for(i=0;i<_tsize;i++){
		for (j=i+1;j<_tsize;j++){
			if (is_subset(tvec[i],tvec[j]) && !is_inter(child[i],parent[j]))
				dirpars[j].insert(i);
		}

	}
	for(i=0;i<_tsize;i++){
		if(dirpars[i].empty())
			dirpar[i]=-1;
		else
			dirpar[i]=*dirpars[i].begin();
	}
	//get direct child in tree
	for (i=0;i<_tsize;i++){
		if(dirpar[i]!=-1){
			dirch[dirpar[i]].insert(i);
		}
	}
	//add empty node
	set<string> enode;
	tvec.push_back(enode);
	parent.resize(_tsize+1);
	child.resize(_tsize+1);
	dirch.resize(_tsize+1);
	dirpar.resize(_tsize+1);
	dirpar[_tsize]=-1;

	for(i=0;i<_tsize;i++){
		parent[i].insert(_tsize);
		child[_tsize].insert(i);
		if(dirpar[i]==-1){
			dirpar[i]=_tsize;
			dirch[_tsize].insert(i);
		}
	}
	//add virtual node
	for(i=0;i<_tsize;){
		int par=dirpar[i];
		set<string> subr=set_sub(tvec[i],tvec[par]);
		if(subr.size()>1){
			set<int> sib=dirch[par];
			map<string,int> count;
			for(set<int>::iterator it1=sib.begin();it1!=sib.end();it1++){
				set<string> t=set_sub(tvec[*it1],tvec[par]);
				if(t.size()>1){
					for(set<string>::iterator it2=t.begin();it2!=t.end();it2++){
						count[*it2]++;
					}
				}
			}
			int max_count=0;
			string max_str;
			for(map<string,int>::iterator it3=count.begin();it3!=count.end();it3++){
				if(it3->second>max_count){
					max_str=it3->first;
					max_count=it3->second;
				}
			}
			set<string> vnode=tvec[par];//new virtual node
			vnode.insert(max_str);
			tvec.push_back(vnode);
			//set parent child dirch dirpar
			int vid=tvec.size()-1;
			parent.resize(tvec.size());
			child.resize(tvec.size());
			dirch.resize(tvec.size());
			dirpar.resize(tvec.size());
			for(j=0;j<vid;j++){
				if(is_subset(tvec[j],vnode)){
					parent[vid].insert(j);
					child[j].insert(vid);
				}else if(is_subset(vnode,tvec[j])){
					child[vid].insert(j);
					parent[j].insert(vid);
				}
			}
			//direct
			dirch[par].insert(vid);
			dirpar[vid]=par;
			for(set<int>::iterator it1=sib.begin();it1!=sib.end();it1++){
				if(is_subset(vnode,tvec[*it1])){
					dirpar[*it1]=vid;
					dirch[par].erase(*it1);
					dirch[vid].insert(*it1);
				}
			}
		}else{
			i++;
		}
	}
	//get diff
	vector<set<string> > new_tvec(tvec.size());
	for(i=0;i<tvec.size();i++){
		if(!tvec[i].empty()){
			new_tvec[i]=set_sub(tvec[i],tvec[dirpar[i]]);
		}
	}
	tvec=new_tvec;
	//generate templet graph, templet tree, templet path
	_gnode.resize(tvec.size());
	for(i=0;i<tvec.size();i++){
		graph_node &gn=_gnode[i];
		gn.id=i;
		gn.child.resize(dirch[i].size());
		gn.follow.resize(child[i].size());
		gn.par=dirpar[i];
		gn.type=-1;
		if(!tvec[i].empty()){
			char add[8192];
			strcpy(add,(tvec[i].begin())->c_str());
			if(add[0]=='p'){
				gn.type=0;
				char *p=strstr(add,",");*p=0;
				gn.row=atoi(add+2);
				gn.col=atoi(p+1);
			}else if(add[0]=='b'){
				gn.type=1;
				char *p=strstr(add,",");*p=0;
				gn.row=atoi(add+2);
				gn.col=atoi(p+1);
			}else if(add[0]=='c'){
				gn.type=2;
				char *p=strstr(add,",");*p=0;
				gn.row=atoi(add+2);
				gn.col=atoi(p+1);
			}else if(add[0]=='r'){
				gn.type=3;
				char *p=strstr(add,",");*p=0;
				gn.row=atoi(add+2);
				gn.col=atoi(p+1);
			}else if(add[0]=='d'){
				gn.type=4;
			}else if(add[0]=='l'){
				gn.type=5;
			}else if(add[0]=='y'){
				gn.type=6;
			}
		}

		j=0;
		for(set<int>::iterator it1=dirch[i].begin();it1!=dirch[i].end();it1++,j++){
			gn.child[j]=*it1;
			_gnode[gn.child[j]].sib=j;
		}
		j=0;
		for(set<int>::iterator it1=child[i].begin();it1!=child[i].end();it1++,j++){
			gn.follow[j]=*it1;
		}
	}
	//set templet sig

	_has_par.clear();
	_has_child.clear();
	_has_root.clear();
	_has_between.clear();
	_has_distance.clear();
	_has_direction.clear();
	_has_label.clear();
	_par_only.clear();
	_child_only.clear();
	_par_child.clear();
	_ysz.clear();

	_has_par.resize(_tsize,false);
	_has_child.resize(_tsize,false);
	_has_root.resize(_tsize,false);
	_has_between.resize(_tsize,false);
	_has_distance.resize(_tsize,false);
	_has_direction.resize(_tsize,false);
	_has_label.resize(_tsize,false);
	_par_only.resize(_tsize,false);
	_child_only.resize(_tsize,false);
	_par_child.resize(_tsize,false);
	_ysz.resize(_tsize,base_feature::_labels.size());

	for(i=0;i<_tsize;i++){
		j=i;
		//type: 0 p, 1 b, 2 c, 3 r, 4 dir, 5 dis, 6 type
		while(j!=_tsize){
			switch(_gnode[j].type){
				case 0:	_has_par[i]=true; break;
				case 1:	_has_between[i]=true; break;
				case 2:	_has_child[i]=true; break;
				case 3:	_has_root[i]=true; _ysz[i]=1; break;
				case 4:	_has_direction[i]=true; break;
				case 5:	_has_distance[i]=true; break;
				case 6:	_has_label[i]=true; break;
			}
			j=_gnode[j].par;
		}
		if(_has_par[i] && _has_child[i] || _has_distance[i] || _has_direction[i])
			_par_child[i]=true;
		else if(_has_par[i] && !_has_child[i] && !_has_distance[i] && !_has_direction[i])
			_par_only[i]=true;
		else if(_has_child[i] && !_has_par[i] && !_has_distance[i] && !_has_direction[i])
			_child_only[i]=true;
	}

	//get templet path
	_templet_path.resize(_tsize);
	for(i=0;i<_tsize;i++){
		j=i;
		while(_gnode[j].par!=-1){		
			_templet_path[i].insert(_templet_path[i].begin(),_gnode[j].sib);
			j=_gnode[j].par;
		}
	}
	//get scan order
	_scan_order.clear();
	queue<int> gq;
	gq.push(_tsize);
	while(!gq.empty()){
		i=gq.front();
		_scan_order.push_back(i);
		gq.pop();
		for(j=0;j<_gnode[i].child.size();j++){
			gq.push(_gnode[i].child[j]);
		}
	}
	//get gnode cluster
	_gcluster.clear();
	for(i=0;i<_gnode.size();i++){
		int cluster_id=-1;
		for(j=0;j<_gcluster.size();j++){
			if((_gnode[i].type==-1 || _gnode[i].type==4 || _gnode[i].type==5 || _gnode[i].type==6) && _gnode[i].type==_gcluster[j].type)
				cluster_id=j;
			else if(_gcluster[j].col==_gnode[i].col && _gcluster[j].row==_gnode[i].row && _gcluster[j].type==_gnode[i].type)
				cluster_id=j;
		}
		if(cluster_id==-1){
			cluster_id=_gcluster.size();
			_gcluster.resize(cluster_id+1);
			_gcluster.back().type=_gnode[i].type;
			if(!(_gnode[i].type==4 || _gnode[i].type==5 || _gnode[i].type==6)){
				_gcluster.back().row=_gnode[i].row;
				_gcluster.back().col=_gnode[i].col;
			}
		}
		_gnode[i].cluster_id=cluster_id;
	}

	//set scan unit
	scan_unit su;
	int child_id=_gnode[_tsize].child.size();
	for(i=0;i<_scan_order.size();i++){
		int gid=_scan_order[i];
		graph_node &gn=_gnode[gid];
		for(j=0;j<gn.child.size();j++){
			//scan_unit su;
			su.child_id=child_id;
			su.gid=gid;
			su.child_size=_gnode[gn.child[j]].child.size();
			if(gn.child[j]<_tsize)
				su.child_size++;
			child_id+=su.child_size;
			su.cluster_id=gn.cluster_id;
			su.child_only=false;
			su.par_child=false;
			su.has_root=false;
			su.par_only=false;
		
			_su.push_back(su);
		}
		if(gid<_tsize){
			su.child_id=0;
			su.gid=gid;
			su.child_size=0;
			su.cluster_id=gn.cluster_id;

			su.child_only=_child_only[gid];
			su.par_child=_par_child[gid];
			su.has_root=_has_root[gid];
			su.par_only=_par_only[gid];
			_su.push_back(su);

		}
	}
	_bases.resize(_su.size());
	_child_size.resize(_su.size());
	_cluster_id.resize(_su.size());
	_que.resize(_su.size());
	_root->graph_id=_tsize;
}

void templet_feature::table2id(vector<vector<string> > & table, vector<vector<int> > & index_table){//for test
	int index,insert_pos;
	index_table.resize(table.size());
	for(int i=0;i<table.size();i++){
		index_table[i].resize(_da.size());
		for(int j=0;j<_da.size();j++){
			const char *p=table[i][j].c_str();
			if(table[i][j]=="")
				index_table[i][j]=-1;
			else if(!_da[j].search((unsigned char *)p,index_table[i][j]))
				index_table[i][j]=-1;
		}
	}
}


int templet_feature::get_csize(graph_node &gn){
	int csize=0;
	if(gn.type==-1)
		csize=1;
	else if(gn.type<4)
		csize=_da[gn.col].size();
	else if(gn.type==4)//%d direction
		csize=2;
	else if(gn.type==5)//%l
		csize=6;//1,2,3,4 : 0,1,2,3    5-9: 4   10- : 5
	else if(gn.type==6)//%y
		csize=base_feature::_labels.size();
	return csize;
}



bool templet_feature::add_feature(vector<vector<int> > &tids, int par_pos, int child_pos, vector<vector<feature> > &vf){
	int i,j,k;
	int rows=tids.size();
	int min_pos=_min(par_pos,child_pos);
	int max_pos=_max(par_pos,child_pos);
	//scan templet_tree
	//get ids
	
	vector<int> ids(_gcluster.size());
	for(i=0;i<_gcluster.size();i++){
		graph_node &gc=_gcluster[i];
		int &id=ids[i];
		if( (gc.type==4||gc.type==5) && par_pos==child_pos)
			continue;
		int row,col;
		if(gc.type==-1){//%
			id=0;
		}else if(gc.type==0){//%p
			id=tids[par_pos][i];
		}else if(gc.type==1){//%b
		}else if(gc.type==2){//%c
			id=tids[child_pos][i];
		}else if(gc.type==3){//%r
			id=tids[child_pos][i];
		}else if(gc.type==4){//%d
			if(child_pos<par_pos)
				id=0;
			else
				id=1;
		}else if(gc.type==5){//%l
			if(max_pos-min_pos<5)
				id=max_pos-min_pos-1;
			else if(max_pos-min_pos<10)
				id=4;
			else
				id=5;
		}else if(gc.type==6){//%y
			//id=_yset;
		}
	}
	//scan
	
	bool par_child= (par_pos==child_pos);
	int h=0;//get here
	int t=_gnode[_tsize].child.size();//put here
	while(h<t){
		scan_unit &su=_su[_que[h]];
		int child_size=su.child_size;
		int cluster_id=su.cluster_id;
		subscript base=_bases[h];
		if(!base){
			h++;
			continue;
		}
		int id=ids[cluster_id];
		if(id!=-1 && base==_dat[base+id].check){
			if(child_size){
				for(i=0;i<child_size;i++){

					if(_dat[base+id].base[i]){
						_bases[t]=_dat[base+id].base[i];
						_que[t]=su.child_id+i;
						t++;
					}
				}
			}else{//feature index
				feature f;
				if(su.par_child && !par_child){
					int grid=par_pos * rows + child_pos;
					for(k=0;k<_ysz[su.gid];){
						f.index=((int *)_dat[base+id].base)[k];
						if(f.index>=0){
							f.y=k;
							vf[grid].push_back(f);
							k++;
						}else{
							k+=-f.index;
						}
					}
				}else if(par_child){
					if(su.par_only){
						int grid=par_pos+rows * rows;
						for(k=0;k<_ysz[su.gid];){
							f.index=((int *)_dat[base+id].base)[k];
							if(f.index>=0){
								f.y=k;
								vf[grid].push_back(f);
								k++;
							}else{
								k+=-f.index;
							}
						}
					}else if(su.child_only){
						int grid=child_pos + rows + rows*rows;
						for(k=0;k<_ysz[su.gid];){
							f.index=((int *)_dat[base+id].base)[k];
							if(f.index>=0){
								f.y=k;
								vf[grid].push_back(f);
								k++;
							}else{
								k+=-f.index;
							}
						}
						
					}else if(su.has_root){
						int grid=child_pos * rows + child_pos;
						f.index=((int *)_dat[base+id].base)[0];
						if(f.index>=0){
							f.y=-1;
							vf[grid].push_back(f);
						}
					}
				}
			}
		}
		h++;
	}
	return true;
}

templet_feature::~templet_feature(){
	if(_dat)	delete [] _dat;
	if(_used)	delete [] _used;
}


void templet_feature::generate_feature_candidate(char *training_file){
	vector<vector<string> > table;
	int sens=0;
	int lines=0;
	int i,j;
	ifstream fin;
	fin.open(training_file);
	bool has_next;
	do{
		has_next=build_table(fin,table,lines);
		if(table.size()){
			add_feature(table);
			table.clear();
			sens++;
			if(!(sens%100))
				cout<<sens<<".. ";
		}
	}while(has_next);
	fin.close();
	//remove features occurs less than _x_freq
	cout<<"shrinking features"<<endl;
	map<char *, int, str_cmp>::iterator p,q;
	_feature_num=0;
	for(p=_feature2id.begin();p!=_feature2id.end();){
		if(p->second<_x_freq){
			q=p;
			p++;
			_feature2id.erase(q);
		}else{
			p->second=_feature_num;
			_feature_num++;
			p++;
		}
	}
}

bool templet_feature::add_feature(vector<vector<string> > &table){
	int i,j,k;
	int rows=table.size();
	for(i=0;i<rows;i++){
		for(j=0;j<_tsize;j++){
			int par_pos=atoi(table[i][table[i].size()-2].c_str());
			int child_pos=i;
			int y=-1;
			if(par_pos<rows)
				y=base_feature::get_label_index(table[child_pos].back());
			vector<string> fs;
			if(get_feature_string(table,j, par_pos,child_pos,y, fs)){
				for(k=0;k<fs.size();k++)
					insert_feature(fs[k]);
			}
		}
	}
	return true;
}

bool templet_feature::insert_feature(string &fs){
	map<char *, int , str_cmp>::iterator p;
	char t[1024];
	strcpy(t,fs.c_str());
	p=_feature2id.find(t);
	if(p==_feature2id.end()){
		char *s=(char *)_feature_str.push_back((char *)t,strlen(t)+1);
		_feature2id.insert(make_pair(s,1));
		return true;
	}else{
		p->second++;
		return false;
	}
}




bool templet_feature::search_feature(string &fs, int &index){
	map<char *, int , str_cmp>::iterator p;
	char t[1024];
	strcpy(t,fs.c_str());
	p=_feature2id.find(t);
	if(p==_feature2id.end())
		return false;
	index=p->second;
	return true;
}

bool templet_feature::get_feature_string(vector<vector<string> > & table, int templet_id, int par_pos, int child_pos, int y, vector<string> &fs){
	int i,j,k;
	int rows=table.size();
	char s[8192];
	int index1,index2;
	fs.clear();
	vector<int> &path=_templet_path[templet_id];
	int max_pos=_max<int>(child_pos,par_pos);
	int min_pos=_min<int>(child_pos,par_pos);

	if(par_pos==rows && !_has_root[templet_id] || _has_root[templet_id] && par_pos!=rows)
		return false;
	//sprintf(s,"%d:",templet_id);
	fs.push_back("");
	for(i=templet_id;i!=_tsize;){
		graph_node &gn=_gnode[i];
//		assert(gn.type!=-1);
		vector<string> adds;
		if(gn.type==0){//%p
			index2=gn.col;
			index1=gn.row+par_pos;
			if(index1<0){
				index1=-index1-1;
				sprintf(s,"B_%d",index1);
				adds.push_back(s);
			}else if(index1>=rows){
				index1-=rows;
				sprintf(s,"E_%d",index1);
				adds.push_back(s);
			}else if(table[index1][index2].size()){
				adds.push_back(table[index1][index2]);
			}
		}else if(gn.type==1){//%b
			index2=gn.col;
			for(j=min_pos+1;j<max_pos;j++){
				index1=gn.row+j;
				index2=gn.col;
				if(index1<0){
					index1=-index1-1;
					sprintf(s,"B_%d",index1);
					adds.push_back(s);
				}else if(index1>=rows){
					index1-=rows;
					sprintf(s,"E_%d",index1);
					adds.push_back(s);
				}else if(table[index1][index2].size()){
					adds.push_back(table[index1][index2]);
				}
			}
		}else if(gn.type==2 || gn.type==3){//%c,%r
			index2=gn.col;
			index1=gn.row+child_pos;
			if(index1<0){
				index1=-index1-1;
				sprintf(s,"B_%d",index1);
				adds.push_back(s);
			}else if(index1>=rows){
				index1-=rows;
				sprintf(s,"E_%d",index1);
				adds.push_back(s);
			}else if(table[index1][index2].size()){
				adds.push_back(table[index1][index2]);
			}
		}else if(gn.type==4){//%d
			if(child_pos<par_pos)// <-
				adds.push_back("0");
			else
				adds.push_back("1");
		}else if(gn.type==5){//%l
			int dis;
			if(max_pos-min_pos<5)
				dis=max_pos-min_pos-1;
			else if(max_pos-min_pos<10)
				dis=4;
			else
				dis=5;
			sprintf(s,"%d",dis);
			adds.push_back(s);
		}else if(gn.type==6){//%y
			adds.push_back(base_feature::_labels[y]);
		}
		//string catch
		if(adds.empty())
			return false;
		vector<string> new_fs;
		for(j=0;j<adds.size();j++){
			for(k=0;k<fs.size();k++){
				if(fs[k].size())
					new_fs.push_back(adds[j]+CUT+fs[k]);
				else
					new_fs.push_back(adds[j]);
			}
		}
		fs=new_fs;
		i=gn.par;
	}
	//add templet id and label
	for(i=0;i<fs.size();i++){
		if(!_has_root[templet_id])
			sprintf(s,"%d:%s%s%s",templet_id,fs[i].c_str(),CUT,base_feature::_labels[y].c_str());
		else
			sprintf(s,"%d:%s",templet_id,fs[i].c_str());
		fs[i]=s;
	}
	return true;
}

bool templet_feature::write_model(char *model_file){
	ofstream fout;
	fout.open(model_file,ios::binary);
	if(!fout.is_open()){
		cout<<"can not open model file: "<<model_file<<endl;
		return false;
	}
	int i,j,k;
	int l;
	char inbit=sizeof(char *);
	fout.write((char *)&inbit,1);
	//write templets
	fout.write((char *)&_tsize,sizeof(int));
	for(i=0;i<_tsize;i++){
		l=_templet_str[i].length();
		fout.write((char *)&l,sizeof(int));
		fout.write(_templet_str[i].c_str(),l);
	}
	fout.write((char *)&_feature_num,sizeof(int));
	map<char *, int, str_cmp>::iterator it;
	for(it = _feature2id.begin(); it != _feature2id.end(); it++){
		l=strlen(it->first);
		fout.write((char *)&l,sizeof(int));
		fout.write(it->first,l);
		fout.write((char *)&(it->second),sizeof(int));
	}

	//construct index
	char line[8192];
	int index;
	vector<set<string> > vss;
	for(it = _feature2id.begin(); it != _feature2id.end(); it++){
		index=it->second;
		strcpy(line,it->first);
		char *t=strstr(line,":");
		*t=0;
		int templet_id=atoi(line);
		t++;
		vector<char *> columns;
		split_string(t,CUT,columns);
		k=templet_id;
		if(_has_root[templet_id])
			j=columns.size()-1;
		else
			j=columns.size()-2;
		for(;j>=0;j--){
			graph_node &gn=_gnode[k];
			k=gn.par;
			if(gn.type>3)
				continue;
			int col=gn.col;
			int row=gn.row;
			if(vss.size()<col+1)
				vss.resize(col+1);
			vss[col].insert(columns[j]);
		}
	}
	
	vector<int> sz(vss.size());
	fill(sz.begin(),sz.end(),0);
	for(i=0;i<vss.size();i++)
		sz[i]=vss[i].size();

	l=vss.size();
	_da.resize(l);
	for(i=0;i<vss.size();i++){
		vector<char *> keys(sz[i]);
		freelist<char> key_str;
		key_str.set_size(PAGESIZE);
		vector<int> values(sz[i]);
		set<string>::iterator it;
		j=0;
		for(it=vss[i].begin();it!=vss[i].end();it++,j++){
			keys[j]=key_str.push_back((char *)it->c_str(),strlen(it->c_str())+1);
			values[j]=j;
		}
		if(keys.size())
			_da[i].build(keys.size(), (unsigned char **)(&keys[0]), &values[0]); 
		else
			_da[i].build(keys.size(), NULL, NULL); 

	}
	vss.clear();
	fout.write((char *)&l,sizeof(int));
	for(i=0;i<l;i++)
		_da[i].write(fout);
	//convert model file to trie tree
	for(i=0;i<_gnode.size();i++){
		graph_node &gn=_gnode[i];
		if(gn.type>3 || gn.type==-1)
			continue;
		if(gn.row<0){
			if(gn.col+1>_table_head.size())
				_table_head.resize(gn.col+1);
			if(_table_head[gn.col].size()<-gn.row)
				_table_head[gn.col].resize(-gn.row,-1);
		}else if(gn.row>0){
			if(gn.col+1>_table_tail.size())
				_table_tail.resize(gn.col+1);
			if(_table_tail[gn.col].size()<gn.row)
				_table_tail[gn.col].resize(gn.row,-1);
		}
	}

	cout<<"Building 2D Trie, this may take several minutes."<<endl;
	
	trie_node_tmp *root=new trie_node_tmp;
	//root->start=0;//initial,
	root->graph_id=_tsize;
	for(it = _feature2id.begin(); it != _feature2id.end(); it++){
		strcpy(line,it->first);
		index=it->second;
		vector<char *> columns;
		char *t=strstr(line,":");
		*t=0;
		int templet_id=atoi(line);
		t++;
		split_string(t,CUT,columns);
		k=templet_id;
		vector<int> &path=_templet_path[templet_id];
		vector<int> ids(path.size(),0);
		int y;//offset to index probe
		if(_has_root[templet_id]){
			y=0;
			j=columns.size()-1;
		}else{
			j=columns.size()-2;
			string ystr=columns.back();
			y=base_feature::get_label_index(ystr);
		}
		for(;j>=0;j--){
			graph_node &gn=_gnode[k];
			k=gn.par;
			int pos,insert_pos;
			if(gn.type<4){
				if(!_da[gn.col].search((unsigned char *)(columns[j]),pos))
					pos=-1;//this should not appear
				//get offset
				if(gn.row<0 && !strncmp(columns[j],"B_",2)){
					int minus=atoi(columns[j]+2);
					_table_head[gn.col][minus]=pos;
				}else if(gn.row>0 && !strncmp(columns[j],"E_",2)){
					int addon=atoi(columns[j]+2);
					_table_tail[gn.col][addon]=pos;
				}
			}else if(gn.type==6){//%y
				string ystr=columns[j];
				pos=base_feature::get_label_index(ystr);
			}else{
				pos=atoi(columns[j]);//%l,%d
			}
			ids[j]=pos;
		}
		insert_trie(path,ids,root,y,_ysz[templet_id],index);
	}
	_feature2id.clear();
	_feature_str.clear();
	//convert to dat
	_dat_size=0;
	dat_resize(PAGESIZE);
	_dat_next_pos=1;
	vector<subscript> vi;
	insert_dat(root,vi);
	for(i=0;i<vi.size();i++){
		_bases[i]=vi[i];
		_que[i]=i;
	}
	int max_code_size=0;
	for(i=0;i<_da.size();i++)
		if(max_code_size<_da[i].size())
			max_code_size=_da[i].size();
	dat_resize(_dat_size+(subscript)max_code_size);
	delete [] _used;
	_used=0;
	fout.write((char *)&_dat_size,sizeof(subscript));
	for(i=0;i<_dat_size;i++){
		fout.write((char *)&(_dat[i].base),sizeof(subscript *));
		fout.write((char *)&(_dat[i].check),sizeof(subscript));
	}
	l=_bas.size();
	fout.write((char *)&l,sizeof(int));
	fout.write((char *)&_bas[0],_bas.size()*sizeof(subscript));
	//_bases,_que
	l=_bases.size();
	fout.write((char *)&l,sizeof(int));
	fout.write((char *)&_bases[0],_bases.size()*sizeof(subscript));
	l=_que.size();
	fout.write((char *)&l,sizeof(int));
	fout.write((char *)&_que[0],_que.size()*sizeof(subscript));
	l=_table_head.size();
	fout.write((char *)&l,sizeof(int));
	for(i=0;i<l;i++){
		j=_table_head[i].size();
		fout.write((char *)&j,sizeof(int));
		if(j>0)
			fout.write((char *)&_table_head[i][0],sizeof(int)*j);
	}
	l=_table_tail.size();
	fout.write((char *)&l,sizeof(int));
	for(i=0;i<l;i++){
		j=_table_tail[i].size();
		fout.write((char *)&j,sizeof(int));
		if(j>0)
			fout.write((char *)&_table_tail[i][0],sizeof(int)*j);
	}

	fout.close();


	delete [] _dat;
	_dat=0;
	_bas.clear();
	cout<<"template feature saved"<<endl;
	return true;
}

void templet_feature::refine_feature(bool *fmap){
	map<char *, int , str_cmp>::iterator it,iter;
	int new_index=0;
	for(it=_feature2id.begin();it!=_feature2id.end();){
		if(fmap[it->second]){
			it->second=new_index++;
			it++;
		}else{
			iter=it;
			it++;
			_feature2id.erase(iter);
		}
	}
	_feature_num=new_index;
}



void templet_feature::generate_ids(vector<vector<int> > &index_table, vector<vector<int> > &ids){
	ids.clear();
	int rows=index_table.size();
	ids.resize(rows);
	int i,j,k;
	for(i=0;i<ids.size();i++){
		ids[i].resize(_gcluster.size());
		for(j=0;j<_gcluster.size();j++){
			graph_node &gc=_gcluster[j];

			int &id=ids[i][j];
			int row,col;
			if(gc.type<4 && gc.type!=-1){//%p,b,c,r
				row=i+gc.row;
				col=gc.col;
				if(row>=0 && row<rows)
					id=index_table[row][col];
				else if(row<0)
					id=_table_head[col][-row-1];
				else
					id=_table_tail[col][row-rows];
			}
		}
	}
}

void templet_feature::generate_feature(vector<vector<string> > & table, vector<vector<feature> > &vf){
	int i,j,k,y;
	int rows=table.size();
	int ysize=base_feature::_labels.size();
	vector<vector<int> > index_table;
	vector<vector<int> > ids;
	table2id(table,index_table);
	generate_ids(index_table,ids);
	index_table.clear();
	int child_pos=-1,par_pos=-1;
	for(par_pos=0;par_pos<rows;par_pos++){
		for(child_pos=0;child_pos<rows;child_pos++){
			add_feature(ids,par_pos,child_pos,vf);
		}
	}
	ids.clear();
	
}



bool templet_feature::load_model(char *model_file,bool require_features,bool require_trie){
	ifstream fin;
	int i,j,k;
	char line[MAXSTRLEN];
	fin.open(model_file,ios::binary);
	if(!fin.is_open()){
		cout<<"can not open model file: "<<model_file<<endl;
		return false;
	}
	//load templets
	_templet_str.clear();
	int l;
	char curbit=sizeof(char *);
	char inbit;
	fin.read((char *)&inbit,1);
	fin.read((char *)&l,sizeof(int));
	for(i=0;i<l;i++){
		int tl;
		fin.read((char *)&tl,sizeof(int));
		fin.read(line,tl);
		line[tl]=0;
		_templet_str.push_back(line);
	}
	construct_templet_graph(_templet_str);
	fin.read((char *)&_feature_num,sizeof(int));
	//skip feature2id
	int index;
	for(i=0;i<_feature_num;i++){
		fin.read((char *)&l,sizeof(int));
		fin.read(line,l);
		line[l]=0;
		fin.read((char *)&index,sizeof(int));
		if(require_features){
			char *s=(char *)_feature_str.push_back(line,l+1);
			_feature2id.insert(make_pair(s,index));
		}
	}

	//trie
	if(!require_trie){
		fin.close();
		return true;
	}
	fin.read((char *)&l,sizeof(int));
	_da.resize(l);
	for(i=0;i<l;i++)
		_da[i].load(fin);
	fin.read((char *)&_dat_size,sizeof(subscript));
	_dat=new tunit[_dat_size];
	if(inbit==curbit){
		for(i=0;i<_dat_size;i++){
			fin.read((char *)&(_dat[i].base),sizeof(subscript *));
			fin.read((char *)&(_dat[i].check),sizeof(subscript));
		}
	}else if(inbit==8){
		char buf[4];
		for(i=0;i<_dat_size;i++){
			fin.read((char *)&(_dat[i].base),4);
			fin.read(buf,4);
			fin.read((char *)&(_dat[i].check),4);
		}
	}else{
		for(i=0;i<_dat_size;i++){
			_dat[i].base=0;
			fin.read((char *)&(_dat[i].base),4);
			fin.read((char *)&(_dat[i].check),4);
		}
	}
	

	fin.read((char *)&l,sizeof(int));
	_bas.resize(l);
	fin.read((char *)&_bas[0],_bas.size()*sizeof(subscript));
	//_bases,_que
	l=_bases.size();
	fin.read((char *)&l,sizeof(int));
	_bases.resize(l);
	fin.read((char *)&_bases[0],_bases.size()*sizeof(subscript));
	l=_que.size();
	fin.read((char *)&l,sizeof(int));
	_que.resize(l);
	fin.read((char *)&_que[0],_que.size()*sizeof(subscript));

	
	fin.read((char *)&l,sizeof(int));
	_table_head.resize(l);
	for(i=0;i<l;i++){
		fin.read((char *)&j,sizeof(int));
		_table_head[i].resize(j);
		if(j>0)
			fin.read((char *)&_table_head[i][0],sizeof(int)*j);
	}
	fin.read((char *)&l,sizeof(int));
	_table_tail.resize(l);
	for(i=0;i<l;i++){
		fin.read((char *)&j,sizeof(int));
		_table_tail[i].resize(j);
		if(j>0)
			fin.read((char *)&_table_tail[i][0],sizeof(int)*j);
	}

	fin.close();

	//adjust base offset
	for(i=0;i<_dat_size;i++){
		memcpy((char *)&j,(char *)&_dat[i].base,sizeof(int));
		if(_dat[i].base)
			_dat[i].base=&_bas[j-1];
	}

	cout<<"template feature loaded, template number: "<<_tsize<<endl;
	return true;
}


bool templet_feature::insert_trie(vector<int> &path, vector<int> &ids, trie_node_tmp *root, int y, int ysz, int index){
	trie_node_tmp* p=root;
	int i,j;
	vector<trie_node_tmp *> empty_child;
	
	for(i=0;i<path.size()+1;i++){
		int col;//id
		int row;
		graph_node &gn=_gnode[p->graph_id];
		if(i==0){
			col=0;
			row=path[i];
		}else if(i==path.size()){
			col=ids[i-1];
			row=gn.child.size();//last row
		}else{
			col=ids[i-1];
			row=path[i];
		}
		int pos,insert_pos;
		if(!vector_search(p->index,col,pos,insert_pos)){
			vector_insert(p->index,col,insert_pos);
			vector_insert(p->child,empty_child,insert_pos);
			pos=insert_pos;
		}
		if(p->child[pos].empty()){
			if(gn.id<_tsize){
				p->child[pos].resize(gn.child.size()+1,0);
			}else
				p->child[pos].resize(gn.child.size(),0);
		}
		if(i!=path.size()){
			if(!p->child[pos][row]){
				p->child[pos][row]=new trie_node_tmp();
				//p->child[pos][row]->index.resize(1,ids[i]);
				p->child[pos][row]->graph_id=gn.child[row];
			}
			p=p->child[pos][row];
		}else{
			
			if(!p->child[pos][row]){
				int *tmp=new int[ysz];
				for(j=0;j<ysz;j++)
					tmp[j]=-1;
				tmp[y]=index;
				p->child[pos][row]=(trie_node_tmp*)tmp;
			}else{
				int *tmp=(int *)(p->child[pos][row]);
				tmp[y]=index;
			}
		}
	}
	return true;
}


subscript templet_feature::dat_resize(subscript size){
	if(_dat_size>=size)
		return _dat_size;
	tunit *tmp=new tunit[size];
	bool *bmp=new bool[size];
	if(!tmp || !bmp)
		assert(0);
	if(!tmp)
		return 0;
	else if(!bmp){
		delete [] tmp;
		return 0;
	}

	int i;
	for(i=0;i<_dat_size;i++){
		tmp[i]=_dat[i];
		bmp[i]=_used[i];
	}
	for(;i<size;i++){
		tmp[i].base=0;
		tmp[i].check=0;
		bmp[i]=false;
	}
	delete [] _dat;
	delete [] _used;
	_dat=tmp;
	_used=bmp;
	_dat_size=size;
	return size;
}
void templet_feature::insert_dat(trie_node_tmp *s, vector<subscript> &begins){//return check of s
	int i,j,k;
	graph_node &gn=_gnode[s->graph_id];
	begins.clear();
	int csz=s->graph_id<_tsize?gn.child.size()+1:gn.child.size();
	begins.resize(csz,0);
	for(i=0;i<gn.child.size();i++){
		vector<subscript> tmpindex;
		vector<trie_node_tmp *> tmptnd;
		for(j=0;j<s->child.size();j++)
			if(s->child[j][i]){
				tmpindex.push_back((subscript)s->index[j]);
				tmptnd.push_back(s->child[j][i]);
			}
		if(tmpindex.size()==0)
			continue;
		subscript range=tmpindex.back();
		subscript pos=_max(_dat_next_pos,tmpindex[0]+1);
		int nonzero=0;
		subscript begin=0;
		bool first=true;
		while(1){
			//pos+max_pos+1 = new size
			subscript new_size=_dat_size;
			while(new_size<pos+range+1){
				new_size*=1.05;
			}
			dat_resize(new_size);
			if(_dat[pos].check){
				nonzero++;
			}else if(first){
				_dat_next_pos=pos;//next pos which check == 0
				first=false;
			}
			begin=pos-tmpindex[0];
			if(!_used[begin]){
				for(j=0;j<tmpindex.size();j++)
					if(_dat[begin+tmpindex[j]].check)
						break;
				if(j==tmpindex.size())//such pos is OK
					break;
			}
			pos++;
		}
		if(((double)nonzero)/(pos - _dat_next_pos + 1) >= 0.95)
			_dat_next_pos=pos;

		begins[i]=begin;

		_used[begin]=true;
		for(j=0;j<tmpindex.size();j++){
			_dat[begin+tmpindex[j]].check=begin;
		}
		for(j=0;j<tmpindex.size();j++){
			vector<subscript> cbase;
			insert_dat(tmptnd[j],cbase);
			int l=_bas.size()+1;
			_bas.insert(_bas.end(),cbase.begin(),cbase.end());
			memcpy(&(_dat[begin+tmpindex[j]].base),&l,sizeof(int));
		}
	}
	if(s->graph_id>=_tsize){
		s->child.clear();
		s->index.clear();
		delete s;
		return;
	}

//s->graph_id<_tsize
	
	vector<subscript> tmpindex;
	vector<int *> tmptnd;
	for(j=0;j<s->child.size();j++){
		if(s->child[j].back()){
			tmpindex.push_back((subscript)s->index[j]);
			tmptnd.push_back((int *)(s->child[j].back()));
		}
	}
	if(tmpindex.size()){
		subscript range=tmpindex.back();
		subscript pos=_max(_dat_next_pos,tmpindex[0]+1);
		int nonzero=0;
		subscript begin=0;
		bool first=true;
		while(1){
			//pos+max_pos+1 = new size
			subscript new_size=_dat_size;
			while(new_size<pos+range+1){
				new_size*=1.05;
			}
			dat_resize(new_size);
			if(_dat[pos].check){
				nonzero++;
			}else if(first){
				_dat_next_pos=pos;//next pos which check == 0
				first=false;
			}
			begin=pos-tmpindex[0];
			if(!_used[begin]){
				for(j=0;j<tmpindex.size();j++)
					if(_dat[begin+tmpindex[j]].check)
						break;
				if(j==tmpindex.size())//such pos is OK
					break;
			}
			pos++;
		}
		if(((double)nonzero)/(pos - _dat_next_pos + 1) >= 0.95)
			_dat_next_pos=pos;

		begins.back()=begin;
		_used[begin]=true;

		

		for(j=0;j<tmpindex.size();j++){
			_dat[begin+tmpindex[j]].check=begin;
			int ysz=_ysz[gn.id];
			int *tmp=new int[ysz];
			memcpy(tmp,tmptnd[j],sizeof(int)*ysz);
			//adjust merge -1,-1 => -2, -1
			for(k=0;k<ysz;k++){
				if(tmp[k]==-1){
					int ii;
					for(ii=k;ii<ysz && tmp[ii]==-1;ii++);
					tmp[k]=k-ii;
				}
			}
			int l=_bas.size()+1;
			for(k=0;k<ysz;k++)
				_bas.push_back(tmp[k]);
			
			memcpy((char *)&(_dat[begin+tmpindex[j]].base),(char *)&l,sizeof(int));
			delete [] tmptnd[j];
			delete [] tmp;
		}
	}
	s->child.clear();
	s->index.clear();
	delete s;
	return;
}
