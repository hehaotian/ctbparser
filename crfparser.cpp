#include <iostream>
#include <fstream>
#include <ctime>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef _WIN32
#include <io.h>
#endif
#include "const.h"
#include "crfparser.h"
#include "lbfgs.h"
const long DEFAULT_WORK_BYTES=1024*1024*256;//256M
using namespace std;
crfparser::crfparser(){
	_start=0;
	_total_space=NULL;
	_work_space=NULL;
	_total_bytes=0;
	_fpool=NULL;
	_lambda=NULL;
	_gradient=NULL;
	_thread_num=1;
	_max_iter=10;
	_sigma=1;
	_eta=0.0001;
	_version=10;
	_fpool=new feature_pool();
	_thread_num=1;
	_threads.resize(1);
	_threads[0]._start_i = 0;
	_threads[0]._obj=0;
	_threads[0]._c=this;
	_threads[0]._gradient=NULL;
	_depth=5;
	_margin=0;
	_nbest=1;
	_prior=SPEED_PRIOR;
	_algorithm=PA_ALGORITHM;
	_temp_space.set_size(MAXSENLEN*MAXSENLEN*2*sizeof(void *));
	_training_model=DIRECT_LEARN;
	_work_bytes=DEFAULT_WORK_BYTES;
	_need_scale=false;

	//change
	
}

crfparser::~crfparser(){
	int i;
	if(_work_space){
		delete [] _work_space;
		_work_space=NULL;
	}

	if(_fpool){
		delete _fpool;
		_fpool=NULL;
	}
	if(_lambda){
		delete [] _lambda;
		_lambda=NULL;
	}
	if(_gradient){
		delete [] _gradient;
		_gradient=NULL;
	}

	for(i=1;i<_thread_num;i++){
		if(_threads[i]._gradient){
			delete [] _threads[i]._gradient;
			_threads[i]._gradient=NULL;
		}
	}
	if(!access("__data1",0))
		unlink("__data1");
	if(!access("__data2",0))
		unlink("__data2");
}

bool crfparser::learn(char *training_file, char *model_file){
	int i,j,k;
	cout<<"crfparser version 0."<<_version<<endl;
	//initialize

	if(_algorithm==CRF_ALGORITHM||_algorithm==L1_ALGORITHM){
		_threads.resize(_thread_num);
		for(i=0;i<_thread_num;i++){
			_threads[i]._c=this;
			_threads[i]._start_i = i;
			_threads[i]._obj=0;
			_threads[i]._gradient=NULL;
		}
	}else if(_algorithm==AP_ALGORITHM||_algorithm==PA_ALGORITHM){
		_thread_num=1;
		_threads.resize(1);
		_threads[0]._start_i = 0;
		_threads[0]._obj=0;
		_threads[0]._c=this;
		_threads[0]._gradient=NULL;
		_prior=SPEED_PRIOR;
	}
	if(!_fpool->check_training(training_file))
		return false;
	cout<<"training data checking complete."<<endl;
	
	if(_training_model==DIRECT_LEARN){//often used, for first learn
	if(!learn_api(training_file))
		return false;
	}else if(_training_model==REFINED_LEARN){//load paramter from a learnt model, and relearn
		if(!model_api(NULL,training_file,1))
			return false;
	}
	
	_sen_num=_fpool->_sen_num;
	if(_work_bytes< _fpool->_max_sen_bytes)
		_work_bytes= _fpool->_max_sen_bytes;
	long senbuf_bytes=_work_bytes;
	_work_bytes*=2;//mirror

	_lambda_size=_fpool->_feature_num;
	if(_prior== MEM_PRIOR && _lambda_size*2*sizeof(double)>_work_bytes){
		_work_bytes=_lambda_size*2*sizeof(double);
	}
	cout<<"allocating "<<_work_bytes<<" bytes"<<endl;
	_work_space=new char[_work_bytes];
	cout<<_work_bytes<<" bytes allocated"<<endl;
	write_model(model_file,0);
	double *w0,*w1;
	if(_prior==SPEED_PRIOR){
		w0=NULL;
		w1=NULL;
	}else{
		w0=(double *)_work_space;
		w1=((double *)_work_space)+_lambda_size;
	}
	_lambda=new double[_lambda_size];
	_gradient=new double[_lambda_size];
	
	memset(_lambda,0,_lambda_size*sizeof(double));
	_threads[0]._gradient=_gradient;
	char disp_buff[1000];
	if(_algorithm==CRF_ALGORITHM){
		cout<<"algorithm: crf\nsentence number: "<<_sen_num<<"\nfeature number: "<<_lambda_size<<"\nsigma: "<<_sigma<<"\nmax_iter: "<<_max_iter<<endl;
		for(i=1;i<_thread_num;i++){
			_threads[i]._gradient=new double[_lambda_size];
		}
		if(_prior==MEM_PRIOR)
			unload();
		int converge=0;
		double old_fvalue=INF;
		clock_t start_time=clock();
		LBFGS l;
		bool init_success=l.init(_lambda_size,_depth,_prior);
		if(!init_success){
			if(_prior==SPEED_PRIOR)
				printf("can not allocate enough memory, use option -m 1 to reduce memory requirement\n");
			else
				printf("can not allocate enough memory, use option -a 1 or -a 2 to train\n");
			exit(1);
		}
		int selected=0;
		for(i=0;i<_max_iter;i++){
			//initial
			bool finish=false;
			for(j=0;j<_thread_num;j++){
				memset(_threads[j]._gradient,0,sizeof(double)*_lambda_size);
				_threads[j]._obj=0;
			}
			double fvalue=0;
			while(!finish){
				finish=load_sentence();

			
			
				for(j=0;j<_thread_num;j++)
					_threads[j].start();
				for(j=0;j<_thread_num;j++)
					_threads[j].join();
			}
			for(j=0;j<_thread_num;j++)
				fvalue += _threads[j]._obj;
			for(j=0;j<_lambda_size;j++)
				for(k=1;k<_thread_num;k++)
					_gradient[j]+=_threads[k]._gradient[j];
			//add _sigma||\_lambda||^2
			for(j=0;j<_lambda_size;j++){
				fvalue+=_lambda[j] * _lambda[j] / (2*_sigma);
				_gradient[j]+=_lambda[j]/_sigma;
			}
			
			double diff=i>0? fabs((old_fvalue-fvalue)/old_fvalue) : 1;
			
			selected=0;
			for(j=0;j<_lambda_size;j++)
				if(_lambda[j]!=0)
					selected++;
			
			sprintf(disp_buff,"iter: %d act: %d diff: %lf obj: %lf", i, selected, diff, fvalue);
			cout<<disp_buff<<endl;

			old_fvalue=fvalue;
			if(diff<_eta)
				converge++;
			else
				converge=0;
			if(i==_max_iter||converge==3) break;//success
			
			int ret=l.optimize(_lambda,&fvalue,_gradient,0,w0,w1);
			if(ret<0){
				cout<<"lbfgs error"<<endl;
				break;
			}else if(ret==0){
				break;//success
			}
			if(_prior==MEM_PRIOR)
				load();
		}
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		cout<<"optimization elapse: "<<elapse<<"s,"<<selected<<" features used"<<endl;
	}if(_algorithm==L1_ALGORITHM){
		cout<<"algorithm: l1 crf\nsentence number: "<<_sen_num<<"\nfeature number: "<<_lambda_size<<"\nl1: "<<_sigma<<"\nmax_iter: "<<_max_iter<<endl;
		for(i=1;i<_thread_num;i++){
			_threads[i]._gradient=new double[_lambda_size];
		}
		if(_prior==MEM_PRIOR)
			unload();
		int converge=0;
		double old_fvalue=INF;
		clock_t start_time=clock();
		LBFGS l;
		bool init_success=l.init(_lambda_size,_depth,_prior);
		if(!init_success){
			if(_prior==SPEED_PRIOR)
				printf("can not allocate enough memory, use option -m 1 to reduce memory requirement\n");
			else
				printf("can not allocate enough memory, use option -a 1 or -a 2 to train\n");
			exit(1);
		}
		
		int selected=0;
		for(i=0;i<_max_iter;i++){
			//initial
			bool finish=false;
			for(j=0;j<_thread_num;j++){
				memset(_threads[j]._gradient,0,sizeof(double)*_lambda_size);
				_threads[j]._obj=0;
			}
			double fvalue=0;
			while(!finish){
				finish=load_sentence();
				for(j=0;j<_thread_num;j++)
					_threads[j].start();
				for(j=0;j<_thread_num;j++)
					_threads[j].join();
			}
			for(j=0;j<_thread_num;j++)
				fvalue += _threads[j]._obj;
			for(j=0;j<_lambda_size;j++)
				for(k=1;k<_thread_num;k++)
					_gradient[j]+=_threads[k]._gradient[j];
			for(j=0;j<_lambda_size;j++){
				fvalue+=fabs(_lambda[j]) / _sigma;
				if(_lambda[j]==0.0){
					if(_gradient[j]<-1/_sigma){
						_gradient[j]+=1/_sigma;
					}else if(_gradient[j]>1/_sigma){
						_gradient[j]-=1/_sigma;
					}else
						_gradient[j]=0;
				}else{
					if(_lambda[j]>0){
						_gradient[j]+=1/_sigma;
					}else{
						_gradient[j]-=1/_sigma;
					}
				}
			}
			
			double diff=i>0? fabs((old_fvalue-fvalue)/old_fvalue) : 1;
			
			selected=0;
			for(j=0;j<_lambda_size;j++)
				if(_lambda[j]!=0)
					selected++;
			
			sprintf(disp_buff,"iter: %d act: %d diff: %lf obj: %lf", i, selected, diff, fvalue);
			cout<<disp_buff<<endl;

			old_fvalue=fvalue;
			if(diff<_eta)
				converge++;
			else
				converge=0;
			if(i==_max_iter||converge==3) break;//success
			
			int ret=l.optimize(_lambda,&fvalue,_gradient,1,w0,w1);
			if(ret<0){
				cout<<"lbfgs error"<<endl;
				break;
			}else if(ret==0){
				break;//success
			}
			if(_prior==MEM_PRIOR)
				load();
		}
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		cout<<"optimization elapse: "<<elapse<<"s,"<<selected<<" features used"<<endl;
	}else if(_algorithm==AP_ALGORITHM){
		
		_nbest=1;
		cout<<"algorithm: ap\nsentence number: "<<_sen_num<<"\nfeature number: "<<_lambda_size<<"\nmax_iter: "<<_max_iter<<endl;
		clock_t start_time=clock();
		memset(_gradient,0,_lambda_size*sizeof(double));
		memset(_lambda,0,_lambda_size*sizeof(double));
		_threads[0]._times=_max_iter*_sen_num;
		for(i=0;i<_max_iter;i++){
			bool finish=false;
			double fvalue=0;
			_threads[0]._obj=0;
			while(!finish){
				finish=load_sentence();
				_threads[0].start();
				_threads[0].join();
			}
			
			int selected=0;
			for(j=0;j<_lambda_size;j++)
				if(_gradient[j]!=0)
					selected++;
			fvalue=_threads[0]._obj/_fpool->_total_words;
			sprintf(disp_buff,"iter: %d act: %d err: %lf", i, selected, fvalue);
			cout<<disp_buff<<endl;
		}
		for(i=0;i<_lambda_size;i++)
			_lambda[i]=_gradient[i]/(_max_iter*_sen_num);
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		cout<<"ap optimization elapse: "<<elapse<<" s"<<endl;
	}else if(_algorithm==PA_ALGORITHM){
		double fvalue=0;
		_nbest=1;
		cout<<"algorithm: pa\nsentence number: "<<_sen_num<<"\nfeature number: "<<_lambda_size<<"\nmax_iter: "<<_max_iter<<"\nC: "<<_sigma<<endl;
		clock_t start_time=clock();
		memset(_gradient,0,_lambda_size*sizeof(double));
		memset(_lambda,0,_lambda_size*sizeof(double));
		_threads[0]._times=_max_iter*_sen_num;
		for(i=0;i<_max_iter;i++){
			bool finish=false;
			double fvalue=0;
			_threads[0]._obj=0;
			while(!finish){
				finish=load_sentence();
				_threads[0].start();
				_threads[0].join();
			}
			int selected=0;
			for(j=0;j<_lambda_size;j++)
				if(_gradient[j]!=0)
					selected++;
			fvalue=_threads[0]._obj/_fpool->_total_words;
			sprintf(disp_buff,"iter: %d act: %d err: %lf", i, selected, fvalue);
			cout<<disp_buff<<endl;
		}
		for(i=0;i<_lambda_size;i++)
			_lambda[i]=_gradient[i]/(_max_iter*_sen_num);
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		cout<<"pa optimization elapse: "<<elapse<<" s"<<endl;
	}
	if(_need_scale && (_algorithm==AP_ALGORITHM || _algorithm==PA_ALGORITHM )){
		//scale
		double max_lambda=0;
		for(i=0;i<_lambda_size;i++){
			if(max_lambda<fabs(_lambda[i])){
				max_lambda=fabs(_lambda[i]);
			}
		}
		for(i=0;i<_lambda_size;i++){
			_lambda[i]=_lambda[i]/max_lambda;//original feature function
		}

		clock_t start_time=clock();
		_factor=1;
		int tmp_algorithm=_algorithm;
		_algorithm=FACTOR_ALGORITHM;
		int converge=0;
		double old_fvalue=1e+37;
		double upper=-1;
		double infer=0;
		
		memcpy(_gradient,_lambda, sizeof(double)*_lambda_size);
		for(i=0;1;i++){
			double fvalue=0;
			_factor_gradient=0;
			_threads[0]._obj=0;
			bool finish=false;
			while(!finish){
				finish=load_sentence();
				_threads[0].start();
				_threads[0].join();
			}
			fvalue=_threads[0]._obj;
			_factor_gradient/=_factor;
			
			double diff=i>0? fabs((old_fvalue-fvalue)/old_fvalue) : 1;
			printf("iter: %d diff: %lf obj: %lf\n", i, diff, fvalue);
			
			old_fvalue=fvalue;
			if(fvalue<_eta)
				converge=3;
			else if(diff<_eta)
					converge++;
				else
					converge=0;
			if(converge==3) break;
			//cout<<_factor<<"\t"<<_factor_gradient<<endl;
			if(_factor_gradient<0)//increase
			{
				infer=max(infer,_factor);
				if(upper<0)
					_factor*=2;
				else
					_factor=(upper+infer)/2;
			}else{//decrease
				if(upper>0)
					upper=min(upper,_factor);
				else
					upper=_factor;
				_factor=(upper+infer)/2;
			}
		}
		//cout<<_factor<<endl;
		for(i=0;i<_lambda_size;i++){
			_lambda[i] *= _factor;
		}
		_algorithm=tmp_algorithm;
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		cout<<"scaling elapse: "<<elapse<<" s"<<endl;
	}
	for(i=1;i<_thread_num;i++){
		delete [] _threads[i]._gradient;
		_threads[i]._gradient=NULL;
	}
	write_model(model_file,1);
	delete _fpool;
	_fpool=NULL;
	for(i=0;i<_sen_num;i++){
		char fn[100];
		sprintf(fn,"__sentence%d",i);
		if(!access(fn,0))
			unlink(fn);
	}
	
	return true;
}

bool crfparser::load_sentence(){
	int i;
	long rest=_work_bytes;
	_sentences.clear();
	char *zone=_work_space;
	char *new_zone=_work_space;
	for(i=_start;i<_sen_num;i++){
		char fn[100];
		sprintf(fn,"__sentence%d",i);
		if(!_fpool->load_sentence(zone,new_zone,fn,rest))
			break;
		
		_sentences.push_back(*((sentence *)zone));
		zone=new_zone;
	}
	
	_start=i;
	if(_start==_sen_num){
		_start=0;
		return true;
	}
	return false;
}

bool crfparser::load_model(char *model_file){
	ifstream fin;
	int i;
	char line[PAGESIZE];
	fin.open(model_file);
	if(!fin.is_open()){
		cout<<"can not open model file: "<<model_file<<endl;
		return false;
	}
	fin.getline(line,PAGESIZE-1);
	_version=atoi(line);
	cout<<"model version: 0."<<_version<<endl;
	fin.getline(line,PAGESIZE-1);
	int ysize=atoi(line);
	base_feature::_labels.resize(ysize);
	for(i=0;i<ysize;i++){
		fin.getline(line,PAGESIZE-1);
		base_feature::_labels[i]=line;
	}
	fin.getline(line,PAGESIZE-1);
	fin.getline(line,PAGESIZE-1);
	_lambda_size=atoi(line);
	_lambda=new double [_lambda_size];
	for(i=0;i<_lambda_size;i++){
		fin.getline(line,PAGESIZE-1);
		_lambda[i]=atof(line);
	}
	fin.close();
	if(!model_api(NULL,NULL,2))
		return false;
	cout<<"model loaded, feature number:"<<_lambda_size<<endl;
	return true;
}


bool crfparser::write_model(char *model_file, int part){
	ofstream fout;
	int i;
	char line[PAGESIZE];
	if(part==0){
		fout.open(model_file);
		if(!fout.is_open()){
			cout<<"can not open model file: "<<model_file<<endl;
			return false;
		}
		fout<<_version<<endl;
	}else if(part==1){
		//refine feature, remove feature with zero weight
		bool *fmap=new bool[_lambda_size];
		int new_lambda_size=0;
		for(i=0;i<_lambda_size;i++){
			if(_lambda[i]>EPS || _lambda[i]<-EPS){
				fmap[i]=true;
				new_lambda_size++;
			}else{
				fmap[i]=false;
			}
		}
		double *tmp=new double[_lambda_size];
		memcpy(tmp,_lambda,_lambda_size*sizeof(double));
		delete [] _lambda;
		_lambda=new double[new_lambda_size];
		int index=0;
		for(i=0;i<_lambda_size;i++){
			if(fmap[i]){
				_lambda[index++]=tmp[i];
			}
		}
		delete [] tmp;
		_lambda_size=new_lambda_size;
		vector<int> base_feature_ids=_fpool->_base_feature_ids;
		_fpool->_base_feature_ids.clear();
		vector<bool *> fmaps(base_feature_ids.size());
		for(i=0;i<fmaps.size();i++)
			fmaps[i]=&fmap[base_feature_ids[i]];
		if(!model_api(&fmaps[0],NULL,0)){
			delete [] fmap;
			fmap=NULL;
			cout<<"feature refinement failed"<<endl;
			return false;
		}else{
			delete [] fmap;
			fmap=NULL;
		}
		//end refinement
		fout.open(model_file,ios::app);
		if(!fout.is_open()){
			cout<<"can not open model file: "<<model_file<<endl;
			return false;
		}

		fout<<base_feature::_labels.size()<<endl;
		for(i=0;i<base_feature::_labels.size();i++)
			fout<<base_feature::_labels[i]<<endl;
		fout<<endl;

		fout<<_lambda_size<<endl;
		for(i=0;i<_lambda_size;i++){
			sprintf(line,"%lf",_lambda[i]);
			fout<<line<<endl;
		}
	}
	fout.close();
	return true;
}

void crfparser::tag(vector<vector<string> > &table, vector<vector<int> > &parent, vector<vector<string> > &label, vector<double> &treep, vector<vector<vector<double> > >&edgep){
	sentence_tmp st;
	int len=table.size();
	int i,j,k;
	_fpool->generate_sentence(table,st);


	sentence sen;
	convert2sentence(st,sen);
	
	crfparser_thread &cur_thread=_threads.front();
	cur_thread.build_lattice_max(sen);
	cur_thread.cyk(sen);
	vector<string> labels(len);
	for(i=0;i<cur_thread._par.size();i++){
		parent.push_back(cur_thread._par[i]);
		for(j=0;j<len;j++){
			if(cur_thread._label[i][j]!=-1)
				labels[j]=base_feature::_labels[cur_thread._label[i][j]];
			else
				labels[j]="";
		}
		label.push_back(labels);
	}
	if(_margin){
		double z;
		double s;
		treep.clear();
		int ysize=base_feature::_labels.size();
		cur_thread.build_lattice(sen);
		cur_thread.collapse_lattice(sen);

		cur_thread.inside_outside(sen,z);
		if(_margin==2 || _margin==3){
			edgep.clear();
			edgep.resize(len);
			for(i=0;i<len;i++){
				edgep[i].resize(len);
				for(j=0;j<len;j++){
					edgep[i][j].resize(ysize,0);
				}
			}
			cur_thread.calculate_margin(sen,z);
			for(i=0;i<len;i++)
				for(j=0;j<len;j++)
					for(k=0;k<ysize;k++)
						edgep[i][j][k]=cur_thread._margin[ysize*(i*len+j)+k];
		}
		if(_margin==1 || _margin==3){
			for(i=0;i<cur_thread._par.size();i++){
				cur_thread.assign_label(sen,&(cur_thread._par[i][0]),&(cur_thread._label[i][0]));
				s=cur_thread.tree_cost(sen);
				treep.push_back(exp(s-z));
			}
		}
	}
}

void crfparser::convert2sentence(sentence_tmp &st, sentence &sen){
//used in test
	int i,j;
	int ysize=base_feature::_labels.size();
	vector<int> fnum(ysize);
	_temp_space.free();
	sen.l=st.label.size();
	sen.label=(int *)_temp_space.push_back((char*)&st.label[0],sizeof(int)*st.label.size());
	sen.parent=(int *)_temp_space.push_back((char*)&st.parent[0],sizeof(int)*st.parent.size());
	sen.fnum=(int **)_temp_space.alloc(sizeof(int*)*st.fvector.size());
	sen.fvector=(int ***)_temp_space.alloc(sizeof(int **)*st.fvector.size());
	for(i=0;i<st.fvector.size();i++){
		sen.fnum[i]=(int *)_temp_space.alloc(sizeof(int)*ysize);
		sen.fvector[i]=(int **)_temp_space.alloc(sizeof(int*)*ysize);
		fill(fnum.begin(),fnum.end(),0);
		for(j=0;j<st.fvector[i].size();j++){
			int y=st.fvector[i][j].y;
			if(y<0)
				y=0;
			fnum[y]++;
		}
		for(j=0;j<ysize;j++){
			sen.fnum[i][j]=fnum[j];
			if(sen.fnum[i][j]>0)
				sen.fvector[i][j]=(int *)_temp_space.alloc(sizeof(int)*fnum[j]);
			else
				sen.fvector[i][j]=NULL;
		}
		
		fill(fnum.begin(),fnum.end(),0);
		for(j=0;j<st.fvector[i].size();j++){
			feature &f=st.fvector[i][j];
			int y=f.y;
			if(y<0)
				y=0;
			sen.fvector[i][y][fnum[y]++]=f.index;
		}
	}
}
void crfparser::unload(){
	FILE *fp=fopen("__data1","wb");
	fwrite(_total_space,sizeof(char),_total_bytes/2,fp);
	fclose(fp);
	fp=fopen("__data2","wb");
	fwrite(_total_space+_total_bytes/2,sizeof(char),_total_bytes-_total_bytes/2,fp);
	fclose(fp);
}

void crfparser::load(){
	FILE *fp=fopen("__data1","rb");
	fread(_total_space,sizeof(char),_total_bytes/2,fp);
	fclose(fp);
	fp=fopen("__data2","rb");
	fread(_total_space+_total_bytes/2,sizeof(char),_total_bytes-_total_bytes/2,fp);
	fclose(fp);
}


bool crfparser::set_para(char *para_name, char *para_value){
	int i;
	if(!strcmp(para_name,"sigma"))
		_sigma=atof(para_value);
	else if(!strcmp(para_name,"margin"))
		_margin=atoi(para_value);
	else if(!strcmp(para_name,"work_bytes"))
		_work_bytes=atoi(para_value)*1024*1024;
	else if(!strcmp(para_name,"nbest"))
		_nbest=atoi(para_value);
	else if(!strcmp(para_name,"max_iter"))
		_max_iter=atoi(para_value);
	else if(!strcmp(para_name,"depth"))
		_depth=atoi(para_value);
	else if(!strcmp(para_name,"eta"))
		_eta=atof(para_value);
	else if(!strcmp(para_name,"training_model"))
		_training_model=atoi(para_value);
	else if(!strcmp(para_name,"need_scale"))
		_need_scale=atoi(para_value)==1;
	else if(!strcmp(para_name,"thread_num"))
		_thread_num=atoi(para_value);
	else if(!strcmp(para_name,"algorithm"))
		_algorithm=atoi(para_value);
	else if(!strcmp(para_name,"mem_model"))
		_prior=atoi(para_value);
	return true;
}
