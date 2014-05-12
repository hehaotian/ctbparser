#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#ifdef _WIN32
#include <io.h>
#endif

#include "const.h"
#include "fun.h"
#include "ctbparser.h"
using namespace std;

bool is_cdigit(const char *p){
	if(!p||strlen(p)==0)
		return false;
	bool is_digit=true;
	while(*p){
		if(!(!strncmp(p,"亿",2)||!strncmp(p,"万",2)||!strncmp(p,"百",2)||!strncmp(p,"千",2)||!strncmp(p,"一",2)||!strncmp(p,"二",2)||!strncmp(p,"三",2)||!strncmp(p,"四",2)||!strncmp(p,"五",2)||!strncmp(p,"六",2)||!strncmp(p,"七",2)||!strncmp(p,"八",2)||!strncmp(p,"九",2)||!strncmp(p,"十",2)||!strncmp(p,"〇",2)))
			return false;
		p+=2;
	}
	return true;
}

bool is_digit(const char *p){
	if(!p||strlen(p)==0)
		return false;
	bool is_digit=true;
	while(*p){
		if(!(!strncmp(p,"１",2)||!strncmp(p,"２",2)||!strncmp(p,"３",2)||!strncmp(p,"４",2)||!strncmp(p,"５",2)||!strncmp(p,"６",2)||!strncmp(p,"７",2)||!strncmp(p,"８",2)||!strncmp(p,"９",2)||!strncmp(p,"０",2)))
			return false;
		p+=2;
	}
	return true;
}


bool is_low_letter(const char *p){
	if(!p||strlen(p)==0)
		return false;
	while(*p){
		unsigned char c1=*p;
		unsigned char c2=*(p+1);
		if(c1!=0xA3)
			return false;
		if(c2<0xC1 || c2>0xDA)
			return false;
		p+=2;
	}
	return true;
}

bool is_upper_letter(const char *p){
	if(!p||strlen(p)<2)
		return false;
	while(*p){
		unsigned char c1=*p;
		unsigned char c2=*(p+1);
		if(c1!=0xA3)
			return false;
		if(c2<0xE1 || c2>0xFA)
			return false;
		p+=2;
	}
	return true;
}

ctbparser::ctbparser(){
	_version=12;
	_segsen=false;
	_full=false;
	_ner_model_file[0]=0;
	_seg_model_file[0]=0;
	_pos_model_file[0]=0;
	_parser_model_file[0]=0;
	_dict_file[0]=0;
	_nbest=1;
	_task=TASK_SEG;

	_ner=new CRF();
	_seg=new CRF();
	_pos=new CRF();
	_dict=new Trie();
	_ns=new normalstr();
	_par=new crfparser();
	
}


ctbparser::~ctbparser(){
	if(_ner){
		delete _ner;
		_ner=0;
	}
	if(_seg){
		delete _seg;
		_seg=0;
	}
	if(_pos){
		delete _pos;
		_pos=0;
	}
	if(_ns){
		delete _ns;
		_ns=0;
	}
	if(_par){
		delete _par;
		_par=0;
	}
}


bool ctbparser::load_config(char *fn){
	cout<<"ctbparser version 0."<<_version<<endl;
	char line[MAXSTRLEN];
	ifstream fin;
	fin.open(fn);
	if(!fin.is_open()){
		cout<<"can not open file: "<<fn<<endl;
		return false;
	}
		
	while(!fin.eof()){
		fin.getline(line,MAXSTRLEN-1);
		if(!line[0] || line[0]=='#')
			continue;
		char *p=strchr(line,'\t');
		*p=0;
		p++;
		//set para, seg
		if(!strcmp(line,"nbest")){
			_nbest=atoi(p);
			_seg->set_para("nbest",p);
			_pos->set_para("nbest",p);
			_par->set_para("nbest",p);
			if(_nbest>1){
				_seg->set_para("seqp","1");
				_pos->set_para("seqp","1");
				_par->set_para("margin","1");
			}
		}else if(!strcmp(line,"task")){
			if(!strcmp(p,"seg"))
				_task=TASK_SEG;
			if(!strcmp(p,"pos"))
				_task=TASK_POS;
			if(!strcmp(p,"parse"))
				_task=TASK_PARSE;
		}else if(!strcmp(line,"ner_model_file"))
			strcpy(_ner_model_file,p);
		else if(!strcmp(line,"seg_model_file"))
			strcpy(_seg_model_file,p);
		else if(!strcmp(line,"pos_model_file"))
			strcpy(_pos_model_file,p);
		else if(!strcmp(line,"parser_model_file"))
			strcpy(_parser_model_file,p);
		else if(!strcmp(line,"dict_file"))
			strcpy(_dict_file,p);
		else if(!strcmp(line,"full"))
			_full=atoi(p);
		else if(!strcmp(line,"segsen"))
			_segsen=atoi(p);
	}
	fin.close();
	//init
	
	if(!_ns->load("GBK2GB2312")){
		cout<<"initialization failed"<<endl;
		return false;
	}
	
	if(!_ner->load_model(_ner_model_file)){
		cout<<"initialization failed"<<endl;
		return false;
	}
	if(!_seg->load_model(_seg_model_file)){
		cout<<"initialization failed"<<endl;
		return false;
	}
	if(!_pos->load_model(_pos_model_file)){
		cout<<"initialization failed"<<endl;
		return false;
	}
	if(!_par->load_model(_parser_model_file)){
		cout<<"initialization failed"<<endl;
		return false;
	}
	if(!check_dict()){
		cout<<"initialization failed"<<endl;
		return false;
	}
	if(!_dict->Load("__dict")){
		cout<<"initialization failed"<<endl;
		return false;
	}
	if(!access("__dict",0))
		unlink("__dict");
	

	return true;
}


void ctbparser::get_constraint(vector<word> &constraint_words){
	constraint_words.clear();
	if(!_dict->root)
		return;
	int i,j,k,ii;
	int l=_fulls.size();
	vector<vector<word> > matched_words(l+1);
	for(i=0;i<_fulls.size();i++){
		//dict match
		TRIENODE *p=_dict->root;
		for(j=i;j<l &&(p=_dict->NearSearch((char *)(_fulls[j].c_str()),p));j++){
			if(!p->link){//
				int wl=strlen(p->wordptr->word)/2;
				char w[MAXWORDLEN*2+1]="";
				for(k=i;k<l&&k-i<wl;k++)
					strcat(w,_fulls[k].c_str());
				if(k-i==wl && !strcmp(p->wordptr->word,w)){
					//k=i+wl
					//find best pos
					POSFREQCHAIN *pf=p->wordptr->posfreq;
					double posfreq=-INF;
					char pos[100];
					while(pf){
						if(pf->freq>posfreq){
							posfreq=pf->freq;
							strcpy(pos,pf->pos);
						}
						pf=pf->next;
					}

						
					word w;
					w.left=i;
					w.right=k-1;
					w.pos=pos;
					w.parent=-1;
					w.weight=posfreq;
					matched_words[k].push_back(w);
				}
				break;
			}else if(p->link[0]){
				TRIENODE *pt=p->link[0];
				//find best pos
				POSFREQCHAIN *pf=pt->wordptr->posfreq;
				double posfreq=-INF;
				char pos[100];
				while(pf){
					if(pf->freq>posfreq){
						posfreq=pf->freq;
						strcpy(pos,pf->pos);
					}
					pf=pf->next;
				}
				word w;
				w.left=i;
				w.right=j;
				w.pos=pos;
				w.parent=-1;
				w.weight=posfreq;
				matched_words[j+1].push_back(w);
			}
		}
		//rule match
		int type=0;
		for(j=i;j<_fulls.size() && is_digit(_fulls[j].c_str());j++);
		if(j-i>4){
			word w;
			w.left=i;
			w.right=j-1;
			w.pos="CD";
			w.parent=-1;
			w.weight=LOG_INF;
			matched_words[j].push_back(w);
		}
		for(j=i;j<_fulls.size() && (is_upper_letter(_fulls[j].c_str())||is_low_letter(_fulls[j].c_str())||j>i && ((_fulls[j]=="．")||_fulls[j]=="＠"));j++);
		if(j-i>2){
			word w;
			w.left=i;
			w.right=j-1;
			w.pos="FW";
			w.parent=-1;
			w.weight=LOG_INF;
			matched_words[j].push_back(w);
		}
		
	}
	//NER detection
	vector<vector<vector<string> > > ext_table(_fulls.size());
	vector<vector<vector<double> > > value_table;
	vector<vector<string> > tag;
	vector<double> seqp;
	vector<vector<double> > nodep;

	for(i=0;i<_fulls.size();i++){
		ext_table[i].resize(1);
		ext_table[i][0].push_back(_fulls[i]);
	}
	_ner->tag(ext_table,value_table,tag,seqp,nodep);
	

	
	for(i=0;i<tag.size();i++){
		int last_word_end=-1;
		for(j=0;j<tag[i].size();j++){
			if(tag[i][j][0]=='B'||tag[i][j][0]=='S'){
				word w;
				w.parent=-1;
				w.pos="NR";
				w.left=j;
				for(k=j;k<tag[i].size() && tag[i][k][0]!='E' && tag[i][k][0]!='S'; k++);
				if(k==tag[i].size())
					k--;
				w.right=k;
				w.weight=LOG_INF;
				matched_words[k+1].push_back(w);
			}
		}
	}

	//get optimal
	vector<double> opt_w(l+1);
	vector<vector<word> > opt_word(l+1);
	fill(opt_w.begin(),opt_w.end(),0);
	for(i=1;i<l+1;i++){
		//set default
		opt_word[i]=opt_word[i-1];
		opt_w[i]=opt_w[i-1]-1;

		vector<word> &v=matched_words[i];
		for(j=0;j<v.size();j++){
			if(opt_w[v[j].left]+v[j].weight>opt_w[i]){
				opt_w[i]=opt_w[v[j].left]+v[j].weight;
				opt_word[i]=opt_word[v[j].left];
				opt_word[i].push_back(v[j]);
			}
		}
	}
	constraint_words=opt_word.back();
}
void ctbparser::process_sub(){
	int i,j,k;
	int l;
	word w;
	vector<vector<vector<string> > > ext_table(_fulls.size());
	vector<vector<vector<double> > > value_table;
	vector<vector<string> > tag;
	vector<double> seqp;
	vector<vector<double> > nodep;
	_probs.clear();
	//constraint
	vector<word> cw;
	get_constraint(cw);
	//convert to seg constraint
	vector<int> con_pos;
	vector<string> con_tag;
	for(i=0;i<cw.size();i++){
		w=cw[i];
		l=w.right-w.left+1;
		if(l==1){
			con_pos.push_back(w.left);
			con_tag.push_back("S");
		}else{
			con_pos.push_back(w.left);
			con_tag.push_back("B");
			if(l>2){
				con_pos.push_back(w.left+1);
				con_tag.push_back("C");
			}
			if(l>3){
				con_pos.push_back(w.left+2);
				con_tag.push_back("D");
			}
			for(j=w.left+3;j+1<=w.right;j++){
				con_pos.push_back(j);
				con_tag.push_back("I");
			}
			con_pos.push_back(w.right);
			con_tag.push_back("E");
		}
	}
	//seg
	
	for(i=0;i<_fulls.size();i++){
		ext_table[i].resize(2);
		ext_table[i][0].push_back(_fulls[i]);
		if(is_digit(_fulls[i].c_str())){
			ext_table[i][1].push_back("1");
		}else if(is_cdigit(_fulls[i].c_str())){
			ext_table[i][1].push_back("2");
		}else if(is_upper_letter(_fulls[i].c_str())){
			ext_table[i][1].push_back("3");
		}else if(is_low_letter(_fulls[i].c_str())){
			ext_table[i][1].push_back("4");
		}
	}

	_seg->tag(ext_table,value_table,tag,seqp,nodep,con_pos,con_tag);
	
	_probs=seqp;
	
	_words.clear();
	_words.resize(tag.size());

	w.left=-1;
	w.parent=-1;
	w.right=-1;
	w.pos="";

	
	for(i=0;i<tag.size();i++){
		int last_word_end=-1;
		for(j=0;j<tag[i].size();j++){
			if(tag[i][j]=="E"||tag[i][j]=="S"||j+1==tag[i].size()||j+1<tag[i].size() && tag[i][j+1]=="B"||j+1<tag[i].size() && tag[i][j+1]=="S"){
				w.left=last_word_end+1;
				w.right=j;
				_words[i].push_back(w);
				last_word_end=w.right;
			}
		}
	}
	//pos
	if(_task==TASK_SEG)
		return;
	_probs.clear();
	vector<vector<word> > words=_words;
	_words.clear();
	int n=tag.size();
	for(i=0;i<n;i++){
		con_pos.clear();
		con_tag.clear();
		ext_table.clear();
		tag.clear();
		vector<double> seqp1;
		ext_table.resize(words[i].size());
		nodep.clear();
		int con_i=0;
		for(j=0;j<words[i].size();j++){
			ext_table[j].resize(9);//w,c1,c2,c3,b1,c-1,c-2,b-1,l
			string w;
			int start=words[i][j].left;
			int end=words[i][j].right;
			int l=end-start+1;
			for(k=start;k<=end;k++)
				w+=_fulls[k];
			ext_table[j][0].push_back(w);
			//first 3 grams
			ext_table[j][1].push_back(_fulls[start]);
			if(l>1)
				ext_table[j][2].push_back(_fulls[start+1]);
			if(l>2)
				ext_table[j][3].push_back(_fulls[start+2]);
			//first bigram
			if(l>1)
				ext_table[j][4].push_back(_fulls[start+1]+_fulls[start+1]);
			//last 2 grams
			ext_table[j][5].push_back(_fulls[end]);
			if(l>1)
				ext_table[j][6].push_back(_fulls[end-1]);
			//last bigram
			if(l>1)
				ext_table[j][7].push_back(_fulls[end-1]+_fulls[end]);
			string type;
			if(is_digit(w.c_str()))
				type="1";
			else if(is_cdigit(w.c_str()))
				type="2";
			else if(is_upper_letter(w.c_str()))
				type="3";
			else if(is_low_letter(w.c_str()))
				type="4";
			ext_table[j][8].push_back(type);
			if(con_i<cw.size() && cw[con_i].left==words[i][j].left && cw[con_i].right==words[i][j].right){
				if(cw[con_i].pos!="-"){
					con_pos.push_back(j);
					con_tag.push_back(cw[con_i].pos);
				}
				con_i++;				
			}
		}
		_pos->tag(ext_table,value_table,tag,seqp1,nodep,con_pos,con_tag);
		vector<word> ws=words[i];
		for(j=0;j<tag.size();j++){
			for(k=0;k<tag[j].size();k++){
				ws[k].pos=tag[j][k];
			}
			if(seqp1.size()){
				_probs.push_back(seqp1[j]*seqp[i]);
			}
			_words.push_back(ws);
		}
	}
	if(_nbest>1){
		seqp=_probs;
		words=_words;
		vector<int> sorted_index;
		merge_sort(seqp,sorted_index);
		for(i=0;i<_nbest && i<words.size();i++){
			_words[i]=words[sorted_index[i]];
			_probs[i]=seqp[sorted_index[i]];
		}
		if(_words.size()>_nbest){
			_words.resize(_nbest);
			_probs.resize(_nbest);
		}
	}
	if(_task==TASK_POS)
		return;
	
	seqp=_probs;
	_probs.clear();
	words=_words;
	n=_words.size();
	_words.clear();
	for(i=0;i<n;i++){
		vector<double> treep;
		tag.clear();
		vector < vector < int > > par;
		vector<vector<vector < double > > >edgep;
		vector < vector < string > > table(words[i].size());
		//build table
		for(j=0;j<words[i].size();j++){
			table[j].resize(3);
			string w;
			int start=words[i][j].left;
			int end=words[i][j].right;
			int l=end-start+1;
			for(k=start;k<=end;k++)
				w+=_fulls[k];
			char ls[100];
			sprintf(ls,"%d",l);
			int type=0;
			if(is_digit(w.c_str()))
				type=1;
			else if(is_cdigit(w.c_str()))
				type=2;
			else if(is_upper_letter(w.c_str()))
				type=3;
			else if(is_low_letter(w.c_str()))
				type=4;
			table[j][0]=w;
			table[j][1]=words[i][j].pos;
			if(type==0)
				table[j][2]=ls;
			else
				table[j][2]="";
		}
		_par->tag(table,par,tag,treep,edgep);
		vector<word> ws=words[i];
		for(j=0;j<tag.size();j++){
			for(k=0;k<tag[j].size();k++){
				ws[k].dep=tag[j][k];
				ws[k].parent=par[j][k];
			}
			if(treep.size())
				_probs.push_back(treep[j]*seqp[i]);
			_words.push_back(ws);
		}
	}
	if(_nbest>1){
		seqp=_probs;
		words=_words;
		vector<int> sorted_index;
		merge_sort(seqp,sorted_index);
		for(i=0;i<_nbest && i<words.size();i++){
			_words[i]=words[sorted_index[i]];
			_probs[i]=seqp[sorted_index[i]];
		}
		if(_words.size()>_nbest){
			_words.resize(_nbest);
			_probs.resize(_nbest);
		}
	}
}
bool ctbparser::get_char(char *&p, char *c){
	if(*p==0)
		return false;
	c[0]=*p;
	c[1]=0;
	if(c[0]<0){
		p++;
		if(*p==0)
			return false;
		c[1]=*p;
	}
	p++;
	return true;
}


bool ctbparser::get_char(FILE *&fp, char *c){
	if(!fp)
		return false;
	if((c[0]=fgetc(fp))==EOF)
		return false;
	c[1]=0;
	if(c[0]<0){
		if((c[1]=fgetc(fp))==EOF)
			return false;
	}
	return true;
}

void ctbparser::decode_string(char *in, char *out){
	out[0]=0;
	vector<string> words;
	_halfs.clear();
	_fulls.clear();
	char c[3]={0};
	char C[3];
	char *p=in;
	
	while(get_char(p,c)){
		if(!c[1]&&(c[0]=='\r'||c[0]=='\n')){
			if(_segsen)
				continue;
			else{
				if(_fulls.size()){
					process_sub();
					output(out);
					_fulls.clear();
					_halfs.clear();
				}
			}
		}else{
			_ns->convert_string(c,C);
			_halfs.push_back(c);
			_fulls.push_back(C);
			if(_fulls.size()>=MAXSENLEN){//too long, cut
				process_sub();
				output(out);
				_fulls.clear();
				_halfs.clear();
			}else if(_segsen && StrStr(SENSEP,C)){//。。   。）
				while(get_char(p,c)){
					if((c[1]&&StrStr("”’）＇｀＂］｝。？！，、；",c))||strstr("?!;,)]}\042'",c) && _fulls.size()<MAXSENLEN){
						_ns->convert_string(c,C);
						_halfs.push_back(c);
						_fulls.push_back(C);
					}else if(!c[1]&&(c[0]=='\r'||c[0]=='\n')){
						continue;
					}else{
						process_sub();
						output(out);
						_fulls.clear();
						_halfs.clear();
						_ns->convert_string(c,C);
						_halfs.push_back(c);
						_fulls.push_back(C);
						break;
					}
				}
			}
		}
	}
	if(_fulls.size()){
		process_sub();
		output(out);
		_fulls.clear();
		_halfs.clear();
	}
}


void ctbparser::decode_file(char *fn_in, char *fn_out){

	vector<string> words;
	_halfs.clear();
	_fulls.clear();
	char c[3]={0};
	char C[3];
	FILE *p,*out;
	if((p=fopen(fn_in,"r"))==NULL)
	{
		cout<<"can not open file: "<<fn_in<<endl;
		return;
	}
	if((out=fopen(fn_out,"w"))==NULL)
	{
		cout<<"can not open file: "<<fn_out<<endl;
		fclose(p);
		return;
	}
	
	while(get_char(p,c)){
		if(!c[1]&&(c[0]=='\r'||c[0]=='\n')){
			if(_segsen)
				continue;
			else{
				if(_fulls.size()){
					process_sub();
					output(out);
					_fulls.clear();
					_halfs.clear();
				}
			}
		}else{
			_ns->convert_string(c,C);
			_halfs.push_back(c);
			_fulls.push_back(C);
			if(_fulls.size()>=MAXSENLEN){//too long, cut
				process_sub();
				output(out);
				_fulls.clear();
				_halfs.clear();
			}else if(_segsen && StrStr(SENSEP,C)){//。。   。）
				while(get_char(p,c)){
					if((c[1]&&StrStr("”’）＇｀＂］｝。？！，、；",c))||strstr("?!;,)]}\042'",c) && _fulls.size()<MAXSENLEN){
						_ns->convert_string(c,C);
						_halfs.push_back(c);
						_fulls.push_back(C);
					}else if(!c[1]&&(c[0]=='\r'||c[0]=='\n')){
						continue;
					}else{
						process_sub();
						output(out);
						_fulls.clear();
						_halfs.clear();
						_ns->convert_string(c,C);
						_halfs.push_back(c);
						_fulls.push_back(C);
						break;
					}
				}
			}
		}
	}
	if(_fulls.size()){
		process_sub();
		output(out);
		_fulls.clear();
		_halfs.clear();
	}
	fclose(p);
	fclose(out);
}


void ctbparser::output(char *out){	
	int i,j,k;
	char t[100];
	vector<string> &vc=_full?_fulls:_halfs;
	for(i=0;i<_words.size();i++){
		for(j=0;j<_words[i].size();j++){
			for(k=_words[i][j].left;k<=_words[i][j].right;k++)
				strcat(out,vc[k].c_str());
			if(_task==TASK_POS||_task==TASK_PARSE){
				strcat(out,"/");
				strcat(out,_words[i][j].pos.c_str());
			}
			if(_task==TASK_PARSE){
				strcat(out,"/");
				sprintf(t,"%d",_words[i][j].parent);
				strcat(out,t);
				strcat(out,"/");
				strcat(out,_words[i][j].dep.c_str());
			}
			if(j+1<_words[i].size())
				strcat(out,"  ");
			else{
				if(_nbest>1)
					sprintf(out,"%s  %lf\n",out,_probs[i]);
				else
					strcat(out,"\n");
			}
		}
	}
}

void ctbparser::output(FILE *&fout){	
	int i,j,k;
	char t[100];
	vector<string> &vc=_full?_fulls:_halfs;
	for(i=0;i<_words.size();i++){
		for(j=0;j<_words[i].size();j++){
			for(k=_words[i][j].left;k<=_words[i][j].right;k++)
				fprintf(fout,"%s",vc[k].c_str());
			if(_task==TASK_POS||_task==TASK_PARSE)
				fprintf(fout,"/%s",_words[i][j].pos.c_str());
			if(_task==TASK_PARSE)
				fprintf(fout,"/%d/%s",_words[i][j].parent,_words[i][j].dep.c_str());
			if(j+1<_words[i].size())
				fprintf(fout,"  ");
			else{
				if(_nbest>1)
					fprintf(fout,"  %lf\n",_probs[i]);
				else
					fprintf(fout,"\n");
			}
		}
	}
}

bool ctbparser::check_dict(){
	FILE *fp;
	FILE *fout;
	bool ret=true;
	if((fp=fopen(_dict_file,"r"))==NULL)
		return false;//fail while opening the file
	fout=fopen("__dict","w");
	double freq;
	char pos[POS_MAXLENGTH];
	char word[WORD_MAXLENGTH];
	char word_gb2312[WORD_MAXLENGTH];
	int lines=0;
	while(fscanf(fp,"%s\t%s\t%lf",word,pos,&freq)!=EOF){
		lines++;
		int index;
		int insert_pos;
		char *p=pos;
		if(strcmp(p,"-") && (!vector_search(_pos->tags,p,index,insert_pos,str_cmp()))){
			cout<<"FILE "<<_dict_file<<", line "<<lines<<", POS tag incorrect."<<endl;
			ret=false;
		}
		if(!_ns->certify(word)){
			cout<<"warning: FILE "<<_dict_file<<", line "<<lines<<", word is not coded in GB2312."<<endl;
		}
		_ns->convert_string(word,word_gb2312);
		fprintf(fout,"%s\t%s\t%lf\n",word_gb2312,pos,freq);
	}
	fclose(fp);
	fclose(fout);
	return ret;
}
