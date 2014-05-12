#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctime>
#ifdef _WIN32
#include <io.h>
#endif
#include "crf.h"
#include "fun.h"
#include "lbfgs.h"
#include "const.h"



bool CRF::learn(char* templet_file, char *training_file,char *model_file)
{
	//initialize
	if(algorithm==CRF_ALGORITHM||algorithm==L1_ALGORITHM)
	{
		threads.resize(thread_num);
		for(int i=0;i<thread_num;i++)
		{
			threads[i].c=this;
			threads[i].start_i = i;
			threads[i].obj=0;
			threads[i].gradient=NULL;
		}
	}else if(algorithm==AP_ALGORITHM || algorithm==PA_ALGORITHM){
		thread_num=1;
		threads.resize(1);
		threads[0].start_i = 0;
		threads[0].obj=0;
		threads[0].c=this;
		threads[0].gradient=NULL;
	}
	printf("pocket crf, version 0.%d\n",version);
	if(!load_templet(templet_file))
		return false;
	printf("templates loaded\n");
	set_chain_type();
	if(!check_training(training_file))
		return false;
	if(!generate_features(training_file))
		return false;
	printf("training data loaded\n");
	shrink_feature();
	printf("features shrinked\n");
	if(fmap_tmp.size())
		fmap=&fmap_tmp[0];			//temporary set
	else
		fmap=NULL;
	fmap_size=fmap_tmp.size();
	sequences=&sequences_tmp[0];	//temporary set
	sequence_num=sequences_tmp.size();
	write_model(model_file,true);
	fmap=NULL;
	sequences=NULL;
	//free memory
	tags.clear();
	vector<char *>(tags).swap(tags);
	tag_str.clear();
	xindex.clear();
	map<char *,int,str_cmp>(xindex).swap(xindex);
	x_freq.clear();
	vector<int>(x_freq).swap(x_freq);
	x_str.clear();
	templets.clear();
	vector<templet>(templets).swap(templets);
	vector<vector<vector<int> > >(templet_group).swap(templet_group);
	compress();

		
	//begin learning

	//set global temporary variables

	path_num=pow((double)ysize,order+1);
	node_anum=pow((double)ysize,order);//alpha(beta) number of each node
	head_offset=-log((double)node_anum);
	//initialize

	int i,j,k;

	lambda=new double[lambda_size];
	gradient=new double[lambda_size];
	threads[0].gradient=gradient;
	for(i=0;i<thread_num;i++){
		threads[i].fmap=fmap;
		if(is_real)
			threads[i]._cal_cost=&crf_thread::real_cost;
		else
			threads[i]._cal_cost=&crf_thread::bool_cost;
		if(fmap)
			threads[i]._findex=&crf_thread::map_findex;
		else
			threads[i]._findex=&crf_thread::dir_findex;
	}
	
	if(algorithm==CRF_ALGORITHM)
	{
		if(prior==MEM_PRIOR)
			unload();
		printf("algorithm: crf\nsequence number: %d\nparameter number: %d\nsigma: %g\nfeature frequence threshold: %d\n",sequence_num,lambda_size,sigma,freq_thresh);
		for(i=1;i<thread_num;i++)
			threads[i].gradient=new double[lambda_size];
		
		
		int converge=0;
		double old_fvalue=1e+37;
		clock_t start_time=clock();
		LBFGS l;
		bool init_success=l.init(lambda_size,depth,prior);
		if(!init_success)
		{
			if(prior==SPEED_PRIOR)
			{
				printf("can not allocate enough memory, use option -m 1 to reduce memory requirement\n");
			}else{
				printf("can not allocate enough memory, use option -f to remove less frequent features\n");
			}
			exit(1);
		}
		memset(lambda,0,lambda_size*sizeof(double));
		
		for(i=0;i<max_iter;i++)
		{
			fvalue=0;
			if(chain_type==SIMPLE_CHAIN)
			{
				memset(transit_buff,0,sizeof(double)*(ysize*max_seq_size+path_num*(max_seq_size-1)));
				for(j=0;j<path_num;j++)
				{
					if(fmap){
						for(k=0;k<max_seq_size-1;k++)
							transit_buff[ysize+(ysize+path_num)*k+ysize+j]=lambda[fmap[transit+j]];
					}else{
						for(k=0;k<max_seq_size-1;k++)
							transit_buff[ysize+(ysize+path_num)*k+ysize+j]=lambda[transit+j];
					}
				}
			}
			for(j=0;j<thread_num;j++)
				threads[j].start();
			for(j=0;j<thread_num;j++)
				threads[j].join();

			for(j=0;j<thread_num;j++)
			{
				fvalue += threads[j].obj;
			}
			for(j=0;j<lambda_size;j++)
				for(k=1;k<thread_num;k++)
					gradient[j]+=threads[k].gradient[j];
			
			//add sigma||\lambda||^2
			for(j=0;j<lambda_size;j++)
			{
				fvalue+=lambda[j] * lambda[j] / (2*sigma);
				gradient[j]+=lambda[j]/sigma;
			}
		
			double diff=i>0? fabs((old_fvalue-fvalue)/old_fvalue) : 1;
			printf("iter: %d diff: %lf obj: %lf\n", i, diff, fvalue);
			old_fvalue=fvalue;
			if(diff<eta)
				converge++;
			else
				converge=0;
			if(i==max_iter||converge==3) break;//success
			double *w0,*w1;
			if(prior==SPEED_PRIOR){
				w0=NULL;
				w1=NULL;
			}else{
				w0=(double *)work_space;
				w1=((double *)work_space)+lambda_size;
			}
			
			if(l.optimize(lambda,&fvalue,gradient,0,w0,w1)<=0)
			{
				printf("lbfgs error\n");
				break;//error
			}
			if(prior==MEM_PRIOR)
				load();
		}
		adjust_data();
		//if(prior==MEM_PRIOR)
		//	unload();//save data structure
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		printf("crf optimization elapse: %g s\n",elapse);
	}else if(algorithm==L1_ALGORITHM){
		if(prior==MEM_PRIOR)
			unload();
		printf("algorithm: crf\nsequence number: %d\nparameter number: %d\nl1: %g\nfeature frequence threshold: %d\n",sequence_num,lambda_size,sigma,freq_thresh);
		for(i=1;i<thread_num;i++)
			threads[i].gradient=new double[lambda_size];
		
		int converge=0;
		double old_fvalue=1e+37;
		clock_t start_time=clock();
		LBFGS l;
		bool init_success=l.init(lambda_size,depth,prior);
		if(!init_success)
		{
			if(prior==SPEED_PRIOR)
			{
				printf("can not allocate enough memory, use option -m 1 to reduce memory requirement\n");
			}else{
				printf("can not allocate enough memory, use option -f to remove less frequent features\n");
			}
			exit(1);
		}
		memset(lambda,0,lambda_size*sizeof(double));
		
		for(i=0;i<max_iter;i++)
		{
			fvalue=0;
			if(chain_type==SIMPLE_CHAIN)
			{
				memset(transit_buff,0,sizeof(double)*(ysize*max_seq_size+path_num*(max_seq_size-1)));
				for(j=0;j<path_num;j++)
				{
					if(fmap){
						for(k=0;k<max_seq_size-1;k++)
							transit_buff[ysize+(ysize+path_num)*k+ysize+j]=lambda[fmap[transit+j]];
					}else{
						for(k=0;k<max_seq_size-1;k++)
							transit_buff[ysize+(ysize+path_num)*k+ysize+j]=lambda[transit+j];
					}
				}
			}
			for(j=0;j<thread_num;j++)
				threads[j].start();
			for(j=0;j<thread_num;j++)
				threads[j].join();

			for(j=0;j<thread_num;j++)
			{
				fvalue += threads[j].obj;
			}
			for(j=0;j<lambda_size;j++)
				for(k=1;k<thread_num;k++)
					gradient[j]+=threads[k].gradient[j];
			for(j=0;j<lambda_size;j++)
			{
				fvalue+=fabs( lambda[j]) / sigma;
				if(lambda[j]==0.0)
				{
					if(gradient[j]<-1/sigma)
					{
						gradient[j]+=1/sigma;
					}else if(gradient[j]>1/sigma){
						gradient[j]-=1/sigma;
					}else
						gradient[j]=0;
				}else{
					if(lambda[j]>0)
					{
						gradient[j]+=1/sigma;
					}else{
						gradient[j]-=1/sigma;
					}
				}
			}
			
			double diff=i>0? fabs((old_fvalue-fvalue)/old_fvalue) : 1;
			
			int selected=0;
			for(j=0;j<lambda_size;j++)
			{
				if(lambda[j]!=0)
					selected++;
			}
			printf("iter: %d act: %d diff: %lf obj: %lf\n", i, selected, diff, fvalue);
			old_fvalue=fvalue;
			if(diff<eta)
				converge++;
			else
				converge=0;
			if(i==max_iter||converge==3) break;//success
			double *w0,*w1;
			if(prior==SPEED_PRIOR){
				w0=NULL;
				w1=NULL;
			}else{
				w0=(double *)work_space;
				w1=((double *)work_space)+lambda_size;
			}
			if(l.optimize(lambda,&fvalue,gradient,1,w0,w1)<=0)
			{
				printf("lbfgs error\n");
				break;//error
			}
			if(prior==MEM_PRIOR)
				load();
		}

		adjust_data();
		//if(prior==MEM_PRIOR)
		//	unload();//save data structure
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		printf("crf feature selection elapse: %g s, %d features selected\n",elapse,lambda_size);
	}else if(algorithm==AP_ALGORITHM){
		printf("algorithm: ap\nsequence number: %d\nparameter number: %d\nfeature frequence threshold: %d\n",sequence_num,lambda_size,freq_thresh);
		clock_t start_time=clock();
		memset(gradient,0,lambda_size*sizeof(double));
		memset(lambda,0,lambda_size*sizeof(double));
		threads[0].times=max_iter*sequence_num;
		for(i=0;i<max_iter;i++)
		{
			threads[0].start();
			threads[0].join();
			fvalue=threads[0].obj/total_nodes;
			printf("iter: %d err: %lf \n",i,fvalue);
		}
		for(i=0;i<lambda_size;i++)
			lambda[i]=gradient[i]/(max_iter*sequence_num);
		adjust_data();
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		printf("ap optimization elapse: %g s\n",elapse);
	}else if(algorithm==PA_ALGORITHM){
		printf("algorithm: pa\nsequence number: %d\nparameter number: %d\nC: %g\nfeature frequence threshold: %d\n",sequence_num,lambda_size,sigma,freq_thresh);
		clock_t start_time=clock();
		memset(gradient,0,lambda_size*sizeof(double));
		memset(lambda,0,lambda_size*sizeof(double));
		threads[0].times=max_iter*sequence_num;
		for(i=0;i<max_iter;i++)
		{
			threads[0].start();
			threads[0].join();
			fvalue=threads[0].obj/total_nodes;
			printf("iter: %d err: %lf \n",i,fvalue);
		}
		for(i=0;i<lambda_size;i++)
			lambda[i]=gradient[i]/(max_iter*sequence_num);
		adjust_data();
		double elapse = static_cast<double>(clock() - start_time) / CLOCKS_PER_SEC;
		printf("pa optimization elapse: %g s\n",elapse);
	}

	write_model(model_file,false);
	
	return true;
}

CRF::CRF()
{
	depth=5;
	algorithm=CRF_ALGORITHM;
	transit_buff=NULL;
	work_space=NULL;
	sequences=NULL;
	sequence1s=NULL;
	fmap=NULL;
	sigma=1;
	prior=SPEED_PRIOR;
	freq_thresh=0;
	margin=false;
	seqp=false;
	nbest=1;
	x_str.set_size(PAGESIZE * 16);//each time alloc 16 pages
	tag_str.set_size(256);//each time alloc 256 bytes
	nodes.set_size(PAGESIZE);
	cliques.set_size(PAGESIZE*16);
	clique_node.set_size(PAGESIZE*16);
	node_clique.set_size(PAGESIZE*16);
	clique_feature.set_size(PAGESIZE*16);
	clique_fvalue.set_size(PAGESIZE*16);
	order=-1;
	max_iter=10000;
	eta=0.0001;
	lambda=NULL;
	gradient=NULL;

	is_real=true;
	version=47;//0.47
	thread_num=1;
	threads.resize(thread_num);
	threads.front().start_i = 0;
	threads.front().obj=0;
	threads.front().c=this;
	dat=new double_array_trie();
}

CRF::~CRF()
{
	if(transit_buff)
	{
		delete [] transit_buff;
		transit_buff=NULL;
	}
	if(lambda)
	{
		delete [] lambda;
		lambda=NULL;
	}
	if(gradient)
	{
		delete [] gradient;
		gradient=NULL;
	}
	if(work_space)
	{
		delete [] work_space;
		work_space=NULL;
		sequences=NULL;
		sequence1s=NULL;
		fmap=NULL;
	}
	if(sequences)
	{
		delete [] sequences;
		sequences=NULL;
	}
	if(sequence1s)
	{
		delete [] sequence1s;
		sequence1s=NULL;
	}
	if(fmap)
	{
		delete [] fmap;
		fmap=NULL;
	}
	if(!access("__data1",0))
		unlink("__data1");
	if(!access("__data2",0))
		unlink("__data2");
	if(dat)
	{
		delete dat;
		dat=NULL;
	}
}
void CRF::set_chain_type()
{
	if(order!=1)
	{
		chain_type=GENERAL_CHAIN;
		return;
	}else{
		for(int i=0;i<templets.size();i++)
		{
			
			if(templets[i].y.size()>1 && templets[i].x.size())
			{
				chain_type=FIRST_CHAIN;
				return;
			}
		}
	}
	chain_type=SIMPLE_CHAIN;
}
bool CRF::set_order()
{
	//set groupid,end_of_group,order
	if(!templets.size())// no templets
		return false;
	//set table_head, table_tail
	for(int i=0;i<templets.size();i++)
	{
		for(int j=0;j<templets[i].x.size();j++){
			if(templets[i].x[j].first<0 && -templets[i].x[j].first>table_head.size())
				table_head.resize(-templets[i].x[j].first);
			if(templets[i].x[j].first>0 && templets[i].x[j].first>table_tail.size())
				table_tail.resize(templets[i].x[j].first);
		}
	}

	for(int i=0;i<table_head.size();i++){
		char t[10]={0};
		sprintf(t,"B_%d",i);
		table_head[i]=t;
	}

	for(int i=0;i<table_tail.size();i++){
		char t[10]={0};
		sprintf(t,"E_%d",i);
		table_tail[i]=t;
	}


	order=-templets[0].y[0];
	templets[0].groupid=0;
	templets[0].end_of_group=false;
	for(int i=1;i<templets.size();i++)
	{
		templet &n=templets[i];
		templet &last_n=templets[i-1];
		int j;
		for(j=0;j<last_n.y.size();j++)
		{
			if(last_n.y[j]!=n.y[j])
			{//start of group
				last_n.end_of_group=true;
				n.groupid=last_n.groupid+1;
				if(-n.y[0]>order)		order=-n.y[0];
				break;
			}
		}
		if(j==last_n.y.size())//not start of group
		{
			n.groupid=last_n.groupid;
			last_n.end_of_group=false;
		}
	}
	templets.back().end_of_group=true;
	gsize=templets.back().groupid+1;
	templet_group.resize(gsize);
	return true;
}
bool CRF::set_para(char *para_name, char *para_value)
{
	if(!strcmp(para_name,"sigma"))
		sigma=atof(para_value);
	else if(!strcmp(para_name,"freq_thresh"))
		freq_thresh=atoi(para_value);
	else if(!strcmp(para_name,"margin"))
		margin=atoi(para_value);
	else if(!strcmp(para_name,"seqp"))
		seqp=atoi(para_value);
	else if(!strcmp(para_name,"nbest"))
		nbest=atoi(para_value);
	else if(!strcmp(para_name,"max_iter"))
		max_iter=atoi(para_value);
	else if(!strcmp(para_name,"depth"))
		depth=atoi(para_value);
	else if(!strcmp(para_name,"eta"))
		eta=atof(para_value);
	else if(!strcmp(para_name,"prior"))
		prior=atoi(para_value);
	else if(!strcmp(para_name,"thread_num"))
		thread_num=atoi(para_value);
	else if(!strcmp(para_name,"algorithm"))
		algorithm=atoi(para_value);
	return true;
}

bool CRF::add_templet(char *line)
{
	if(!line[0]||line[0]=='#') 
		return false;

	templet n;
	char *p=line,*q;
	char word[1000];
	

	char index_str[1000];
	int index1,index2;
	while(q=catch_string(p,"%x[",word))
	{
		p=q;
		p=catch_string(p,",",index_str);
		index1=atoi(index_str);
		p=catch_string(p,"]",index_str);
		index2=atoi(index_str);
		n.x.push_back(make_pair(index1,index2));

	}
	q=catch_string(p,"%y[",word);
	if(!q)
	{
		printf("templet: %s incorrect\n",line);
		return false;
	}

	p=q-3;
	while(p=catch_string(p,"%y[","]",index_str))
	{
		index1=atoi(index_str);
		n.y.push_back(index1);
	}
	int insert_pos;
	vector_search(templets,n,index1,insert_pos,templet_cmp());
	vector_insert(templets,n,insert_pos);
	return true;
}
bool CRF::load_templet(char *templet_file)
{
	FILE *fp;
	//read template
	fp=fopen(templet_file,"r");
	if (!fp) 
	{
		printf("template file: %s not found\n",templet_file);
		return false;
	}
	char line[MAXSTRLEN];
	while(fgets(line,MAXSTRLEN-1,fp))
	{
		trim_line(line);
		add_templet(line);
	}
	fclose(fp);
	return set_order();
}





bool CRF::check_training(char *training_file)
{
	FILE *fp;
	if((fp=fopen(training_file,"r"))==NULL)	return false;
	char line[MAXSTRLEN];
	int lines=0;
	cols=0;
	while(fgets(line,MAXSTRLEN-1,fp))
	{
		lines++;
		if(strlen(line)==1) continue;
		trim_line(line);
		vector<char *>columns;
		if(!split_string(line,"\t",columns))
		{
			printf("columns should be greater than 1\n");
			fclose(fp);
			return false;//columns should greater than 1
		}
		if(cols && cols!=columns.size())
		{//incompatible
			printf("line %d: columns incompatible\n",lines);
			fclose(fp);
			return false;
		}
		cols=columns.size();
		if(cols<2){
			printf("columns should be greater than 1\n");
			fclose(fp);
			return false;//columns should greater than 1
		}
		//check if all features are real
		if(is_real){
			for(int i=0;i+1<cols && is_real;i++)
			{
				vector<char *> units;
				split_string(columns[i]," ",units);
				for(int j=0;j<units.size();j++)
				{
					if(units[j][0] && !strrchr(units[j],':'))
					{
							is_real=false;
							break;
					}
				}
			}
		}

		char *t=columns.back();//tag
		int index,insert_pos;
		if(!vector_search(tags,t,index,insert_pos,str_cmp()))
		{
			char *p=tag_str.push_back(t);//copy string
			vector_insert(tags,p,insert_pos);
		}
	}
	fclose(fp);
	ysize=tags.size();
	set_group();
	return true;
}

bool CRF::generate_features(char *filename)
{
	char line[MAXSTRLEN];
	max_seq_size=0;
	total_nodes=0;
	vector<char *>table;// table[i,j] = table[i*cols+j]
	charlist table_str;
	table_str.set_size(PAGESIZE);//1 page memory
	lambda_size=0;//lambda size
	int lines=0;
	int i;
	FILE *fp=fopen(filename,"r");
	while(fgets(line,MAXSTRLEN-1,fp))
	{
		trim_line(line);
		if(line[0])
		{
			vector<char *>columns;
			split_string(line,"\t",columns);
			for(i=0;i<cols;i++)
			{
				char *p=table_str.push_back(columns[i]);
				table.push_back(p);
			}
		}else{//get one sequence
			if(table.size())//non-empty line
			{
				if(max_seq_size<table.size())
					max_seq_size=table.size();
				total_nodes+=table.size();
				add_x(table);
				table.clear();//prepare for new table
				table_str.free();
				lines++;
				if(!(lines%100))
					printf("%d.. ",lines);
			}
		}
	}
	if(table.size())//non-empty line
	{
		if(max_seq_size<table.size())
			max_seq_size=table.size();
		total_nodes+=table.size();
		add_x(table);
		table.clear();//prepare for new table
		table_str.free();
		lines++;
		if(!(lines%100))
			printf("%d.. ",lines);
	}
	fclose(fp);
	//set transit_buff
	if(chain_type==SIMPLE_CHAIN && max_seq_size && (algorithm==CRF_ALGORITHM||algorithm==L1_ALGORITHM))
		transit_buff=new double[ysize*max_seq_size+ysize*ysize*(max_seq_size-1)];
	return true;
}

void CRF::set_group()
{
	//calculate templet_group
	//set the size of each group
	//order and tags.size() of CRF must be known
	int i,j;
	for(i=0,j=0;i<templets.size();i++)
	{

		if(templets[i].end_of_group)
		{
			int n=pow((double)ysize,(int)templets[i].y.size());
			templet_group[j++].resize(n);//group j has n offsets
		}
	}
	vector<int> path_index(order+1,0);
	int path_size=pow((double)ysize,order+1);
	for(i=0;i<path_size;i++)//assosiate path i with templet_group
	{
		int cur_group=0;
		for(j=0;j<templets.size();j++)
		{
			if(templets[j].end_of_group)
			{
				vector<int> &ytemp=templets[j].y;
				int k,offset;
				vector<int> temp;
				for(k=0,offset=0;k<ytemp.size();k++)
					offset=offset*ysize+path_index[-ytemp[k]];

				templet_group[cur_group++][offset].push_back(i);//path i added to current group's offset
			}
		}
		for(j=0;j<order+1 && path_index[j]==ysize-1;j++);
		if(j==order+1) break;
		path_index[j]++;
		for(j--;j>=0;j--)	path_index[j]=0;
	}
	for(i=0;i<templet_group.size();i++)
		for(j=0;j<templet_group[i].size();j++)
		{
			((vector<int>)(templet_group[i][j])).swap(templet_group[i][j]);
		}
}
bool CRF::add_x(vector<char *> &table)
{
	int i,j,k,c;
	int rows=table.size()/cols;
	char s[1024];

	char s1[1024];
	char s2[1024];
	
	sequence seq;
	vector<int> y;
	node* nod=nodes.alloc(rows);
	seq.node_num=rows;
	seq.nodes=nod;
	sequences_tmp.push_back(seq);

	vector<vector<char *> > ext_table(table.size());//split each unit by " "
	vector<vector<double> > val_table(table.size());
	for(i=0;i<table.size();i++){
		vector<char *> units;
		split_string(table[i]," ",units);
		for(j=0;j<units.size();)
		{
			if(!units[j][0])
				units.erase(units.begin()+j);
			else
			{
				if(is_real)
				{
					char *p=strrchr(units[j],':');
					if(p)
					{
						if(p==units[j])//	:0.5
						{
							units.erase(units.begin()+j);
							continue;
						}
						*p=0;
						p++;
						val_table[i].push_back(atof(p));
					}else{//last column
						val_table[i].push_back(0);
					}
				}
				j++;
			}
		}
		ext_table[i]=units;
	}
	//below, using ext_table
	for(i=0;i<rows;i++)
	{
		y.resize(y.size()+1);
		vector_search(tags,ext_table[(i+1)*cols-1][0],y.back(),j,str_cmp());//get the tag of current node,j is invalid
		nod[i].key=0;
		for(j=i-order;j<=i;j++)
			if(j>=0)
				nod[i].key=nod[i].key*ysize+y[j];
		vector<clique*> clisp;//features that affect on current nodes
		
		vector<int> feature_vector;
		vector<double> fvalue_vector;
		for(j=0;j<templets.size();j++)
		{
			//get first y's offset
			templet &pat=templets[j];
			if(pat.y[0]+i<0)
				continue;
			double fval=1;
			if(pat.x.size()>0){//if has x
				int index1,index2;
				bool has_xstring=true;//false, if unit=""
				vector<int> xid(pat.x.size(),0);//xid=(0,0): 0 th units + 0 th units
				vector<int> xtop(pat.x.size(),0);
				//get xtop
				
				for(k=0;k<pat.x.size();k++)
				{
					index1=pat.x[k].first+i;
					index2=pat.x[k].second;
					if(index1<0)
					{
						xtop[k]=1;
					}else if(index1>=rows){
						xtop[k]=1;
					}else if(!ext_table[index1*cols+index2].size()){//no string here
						xtop[k]=0;
						has_xstring=false;
						break;
					}else{
						xtop[k]=ext_table[index1*cols+index2].size();
					}
				}
				if(has_xstring)
				{
					xid.back()=-1;
					while(1)
					{
						//increase and check whether stop
						for(k=xid.size()-1;k>=0 && xid[k]+1 == xtop[k];k--)
							xid[k]=0;
						if(k<0)
							break;//stop
						xid[k]++;
						sprintf(s, "%d", j);
						strcat(s,":");
						
						//get x
						for(k=0;k<pat.x.size();k++)
						{
							index1=pat.x[k].first+i;
							index2=pat.x[k].second;
							assert(index2>=0 && index2<cols-1);
							if(index1<0)
								strcpy(s1,table_head[-index1-1].c_str());
							else if(index1>=rows)
								strcpy(s1,table_tail[index1-rows].c_str());
							else{
								assert(ext_table[index1*cols+index2].size());
								strcpy(s1,ext_table[index1*cols+index2][xid[k]]);
								if(is_real)
									fval*=val_table[index1*cols+index2][xid[k]];
							}
							strcat(s,s1);
							strcat(s,"\001");
						}
						//x obtained, insert x
						int index;//index of feature s
						if(insert_x(s,index))
						{
							c=pow((double)ysize,(int)pat.y.size());
							lambda_size+=c;
							//fmap
							if(freq_thresh)
								fmap_tmp.resize(lambda_size,0);
						}
						//get clique
						feature_vector.push_back(index);
						if(is_real)
							fvalue_vector.push_back(fval);
						fval=1;
					}
				}
			}else{//else , no x
				sprintf(s, "%d", j);
				strcat(s,":");
				//x obtained, insert x
				int index;//index of feature s
				if(insert_x(s,index))
				{
					c=pow((double)ysize,(int)pat.y.size());
					lambda_size+=c;
					if(freq_thresh)
						fmap_tmp.resize(lambda_size,0);
				}
				//get clique
				feature_vector.push_back(index);
				if(is_real)
					fvalue_vector.push_back(fval);
				fval=1;
			}

			if(pat.end_of_group)
			{//creat new clique
				clique cli;
				
				vector<node*> ns;
				int key=0;
				for(k=0;k<pat.y.size();k++)
				{
					ns.push_back(nod+i+pat.y[k]);
					key=key*ysize+ y[i+pat.y[k]];
				}
				//set feature count
				if(freq_thresh)
					for(k=0;k<feature_vector.size();k++)
						fmap_tmp[feature_vector[k]+key]++;
				node ** np=clique_node.push_back(&ns[0],ns.size());
				cli.nodes=np;
				cli.node_num=ns.size();
				cli.key=key;
				int *f=NULL;
				double *fv=NULL;
				if(feature_vector.size()){
					f=clique_feature.push_back(&feature_vector[0],feature_vector.size());
					if(is_real)
						fv=clique_fvalue.push_back(&fvalue_vector[0],fvalue_vector.size());
				}
				cli.fvector=f;
				cli.fvalue=fv;
				cli.feature_num=feature_vector.size();
				cli.groupid=templets[j].groupid;
				clique *new_clique=cliques.push_back(&cli,1);
				clisp.push_back(new_clique);
				feature_vector.clear();
				fvalue_vector.clear();
			}
		}
		//set node -> clique
		if(clisp.size())
			nod[i].cliques = node_clique.push_back(&clisp[0],clisp.size());
		else
			nod[i].cliques = NULL;
		nod[i].clique_num =clisp.size();
	}
	return true;
}

bool CRF::insert_x(char *target, int &index)
{
	map<char *, int , str_cmp>::iterator p;
	p=xindex.find(target);
	if(p!=xindex.end())
	{
		index=p->second;
		x_freq[index]++;
		return false;
	}else{
		char *q=x_str.push_back(target);
		xindex.insert(make_pair(q,lambda_size));
		index=lambda_size;
		x_freq.resize(index+1,0);
		x_freq[index]=1;
		return true;
	}
}



void CRF::compress()
{
	//count bytes
	unsigned long bytes=0;
	int i,j,k,ii;
	if(chain_type==GENERAL_CHAIN){
		//sequence data
		for(i=0;i<sequences_tmp.size();i++)
		{
			sequence &seq=sequences_tmp[i];
			for(j=0;j<seq.node_num;j++)
			{
				node &nod=seq.nodes[j];
				for(k=0;k<nod.clique_num;k++)
				{
					if(!nod.cliques[k])
						continue;
					clique &cli=*nod.cliques[k];
					bytes+=cli.feature_num*sizeof(int);
					if(is_real)
						bytes+=cli.feature_num*sizeof(double);
					bytes+=cli.node_num*sizeof(node*);
					bytes+=sizeof(clique);
				}
				bytes+=nod.clique_num*sizeof(clique*);
			}
			bytes+=seq.node_num*sizeof(node);
			bytes+=sizeof(sequence);
		}
	}else if(chain_type==FIRST_CHAIN){
		//sequence data
		for(i=0;i<sequences_tmp.size();i++)
		{
			sequence &seq=sequences_tmp[i];
			for(j=0;j<seq.node_num;j++)
			{
				node &nod=seq.nodes[j];
				for(k=0;k<nod.clique_num;k++)
				{
					if(!nod.cliques[k])
						continue;
					clique &cli=*nod.cliques[k];
					bytes+=cli.feature_num*sizeof(int);
					if(is_real)bytes+=cli.feature_num*sizeof(double);
				}
			}
			bytes+=seq.node_num*sizeof(vertex);
			bytes+=(seq.node_num-1)*sizeof(edge);
			bytes+=sizeof(sequence1);
		}
	}else if(chain_type==SIMPLE_CHAIN){
		//sequence data
		for(i=0;i<sequences_tmp.size();i++)
		{
			sequence &seq=sequences_tmp[i];
			for(j=0;j<seq.node_num;j++)
			{
				node &nod=seq.nodes[j];
				for(k=0;k<nod.clique_num;k++)
				{
					if(!nod.cliques[k])
						continue;
					clique &cli=*nod.cliques[k];
					if(templet_group[cli.groupid].size()==ysize)//vertex feature
					{
						bytes+=cli.feature_num*sizeof(int);
						if(is_real)bytes+=cli.feature_num*sizeof(double);
					}
				}
			}
			bytes+=seq.node_num*sizeof(vertex);
			bytes+=sizeof(sequence1);
		}
	}
	//fmap
	bytes+=fmap_tmp.size()*sizeof(int);

	if(prior== MEM_PRIOR && (algorithm==CRF_ALGORITHM||algorithm==L1_ALGORITHM) && bytes<lambda_size*sizeof(double)*2)
		bytes=lambda_size*sizeof(double)*2;
	work_size=bytes;
	//allocate
	work_space=new char[work_size];
	char *p=work_space;
	//copy

	//sequence data
	if(chain_type==GENERAL_CHAIN)
	{
		memcpy(p,&sequences_tmp[0],sequences_tmp.size()*sizeof(sequence));
		sequence_num=sequences_tmp.size();
		sequences_tmp.clear();
		vector<sequence>(sequences_tmp).swap(sequences_tmp);
		sequences=(sequence *)p;
		p+=sequence_num*sizeof(sequence);
		for(i=0;i<sequence_num;i++)
		{
			sequence &seq=sequences[i];
			memcpy(p,seq.nodes,sizeof(node)*seq.node_num);
			node *tmp_nodes=seq.nodes;
			seq.nodes=(node *)p;
			p+=sizeof(node)*seq.node_num;

			for(j=0;j<seq.node_num;j++)
			{
				node &nod=seq.nodes[j];
				vector<clique*> clisp(nod.clique_num);
				for(k=0;k<nod.clique_num;k++)
				{
					if(!nod.cliques[k])
					{
						clisp[k]=NULL;
						continue;
					}
					clique &cli=*nod.cliques[k];
					
					if(cli.feature_num)
					{
						memcpy(p,cli.fvector,cli.feature_num*sizeof(int));
						cli.fvector=(int*) p;
						p+=cli.feature_num*sizeof(int);
						if(is_real)
						{
							memcpy(p,cli.fvalue,cli.feature_num*sizeof(double));
							cli.fvalue=(double*) p;
							p+=cli.feature_num*sizeof(double);
						}else{
							cli.fvalue=NULL;
						}
					}//else cli.fvector=NULL;
					

					vector<node*> ns(cli.node_num);
					for(ii=0;ii<cli.node_num;ii++)
						ns[ii]=cli.nodes[ii]-tmp_nodes+seq.nodes;

					memcpy(p,&ns[0],cli.node_num*sizeof(node*));
					cli.nodes=(node**) p;
					p+=cli.node_num*sizeof(node*);

					memcpy(p,&cli,sizeof(clique));
					clisp[k]=(clique*)p;
					p+=sizeof(clique);
				}
				nod.cliques=NULL;
				if(nod.clique_num)
				{
					memcpy(p,&clisp[0],sizeof(clique*)*nod.clique_num);
					nod.cliques=(clique**)p;
					p+=sizeof(clique*)*nod.clique_num;
				}
			}
		}
	}else if(chain_type==FIRST_CHAIN){
		sequence_num=sequences_tmp.size();
		vector<sequence1> seq1s(sequence_num);
		for(i=0;i<sequence_num;i++)
		{
			sequence &seq=sequences_tmp[i];
			seq1s[i].vertex_num=seq.node_num;
			vector<vertex> vertexes(seq.node_num);
			vector<edge> edges(seq.node_num-1);
			for(j=0;j<seq.node_num;j++)
			{
				node &nod=seq.nodes[j];
				for(k=0;k<nod.clique_num;k++)
				{
					if(!nod.cliques[k])
						continue;
					clique &cli=*nod.cliques[k];
					if(templet_group[cli.groupid].size()==ysize)//vertex feature
					{
						if(cli.feature_num)
						{
							memcpy(p,cli.fvector,cli.feature_num*sizeof(int));
							vertexes[j].fvector=(int *)p;
							p+=cli.feature_num*sizeof(int);
							if(is_real)
							{
								memcpy(p,cli.fvalue,cli.feature_num*sizeof(double));
								vertexes[j].fvalue=(double *)p;
								p+=cli.feature_num*sizeof(double);
							}else{
								vertexes[j].fvalue=NULL;
							}
						}else{
							vertexes[j].fvector=NULL;
							vertexes[j].fvalue=NULL;
						}

						vertexes[j].feature_num=cli.feature_num;
						vertexes[j].key=nod.key;
					}else{//edge feature
						if(cli.feature_num)
						{
							memcpy(p,cli.fvector,cli.feature_num*sizeof(int));
							edges[j-1].fvector=(int *)p;
							p+=cli.feature_num*sizeof(int);
							if(is_real)
							{
								memcpy(p,cli.fvalue,cli.feature_num*sizeof(double));
								edges[j-1].fvalue=(double *)p;
								p+=cli.feature_num*sizeof(double);
							}else{
								edges[j-1].fvalue=NULL;
							}
						}else{
							edges[j-1].fvector=NULL;
							edges[j-1].fvalue=NULL;
						}
						edges[j-1].feature_num=cli.feature_num;
					}
				}
			}
			memcpy(p,&vertexes[0],seq.node_num*sizeof(vertex));
			seq1s[i].vertexes=(vertex*)p;
			p+=seq.node_num*sizeof(vertex);

			if(seq.node_num>1){
				memcpy(p,&edges[0],edges.size()*sizeof(edge));
				seq1s[i].edges=(edge*)p;
				p+=edges.size()*sizeof(edge);
			}else{
				seq1s[i].edges=NULL;
			}
		}
		memcpy(p,&seq1s[0],seq1s.size()*sizeof(sequence1));
		sequence1s=(sequence1*)p;
		p+=seq1s.size()*sizeof(sequence1);
		sequences_tmp.clear();
		vector<sequence>(sequences_tmp).swap(sequences_tmp);
	}else if(chain_type==SIMPLE_CHAIN){
		transit=-1;
		sequence_num=sequences_tmp.size();
		vector<sequence1> seq1s(sequence_num);
		for(i=0;i<sequence_num;i++)
		{
			sequence &seq=sequences_tmp[i];
			seq1s[i].vertex_num=seq.node_num;
			vector<vertex> vertexes(seq.node_num);
			for(j=0;j<seq.node_num;j++)
			{
				node &nod=seq.nodes[j];
				for(k=0;k<nod.clique_num;k++)
				{
					if(!nod.cliques[k])
						continue;
					clique &cli=*nod.cliques[k];
					if(templet_group[cli.groupid].size()==ysize)//vertex feature
					{
						if(cli.feature_num)
						{
							memcpy(p,cli.fvector,cli.feature_num*sizeof(int));
							vertexes[j].fvector=(int *)p;
							p+=cli.feature_num*sizeof(int);
							if(is_real)
							{
								memcpy(p,cli.fvalue,cli.feature_num*sizeof(double));
								vertexes[j].fvalue=(double *)p;
								p+=cli.feature_num*sizeof(double);
							}else{
								vertexes[j].fvalue=NULL;
							}
						}else{
							vertexes[j].fvector=NULL;
							vertexes[j].fvalue=NULL;
						}
						vertexes[j].feature_num=cli.feature_num;
						vertexes[j].key=nod.key;
					}else if(transit<0){//edge feature
						transit=cli.fvector[0];
					}
				}
			}
			memcpy(p,&vertexes[0],seq.node_num*sizeof(vertex));
			seq1s[i].vertexes=(vertex*)p;
			p+=seq.node_num*sizeof(vertex);
			seq1s[i].edges=NULL;
		}
		memcpy(p,&seq1s[0],seq1s.size()*sizeof(sequence1));
		sequence1s=(sequence1*)p;
		p+=seq1s.size()*sizeof(sequence1);
		sequences_tmp.clear();
		vector<sequence>(sequences_tmp).swap(sequences_tmp);
	}
	nodes.clear();
	cliques.clear();
	clique_node.clear();//clique affected nodes
	node_clique.clear();//node->clique
	clique_feature.clear();//clique->feature

	
	//fmap
	if(fmap_tmp.size())
	{
		memcpy(p,&fmap_tmp[0],sizeof(int)*fmap_tmp.size());
		fmap=(int *)p;
	}
	fmap_size=fmap_tmp.size();
	fmap_tmp.clear();
	vector<int>(fmap_tmp).swap(fmap_tmp);
	p+=sizeof(int)*fmap_size;
}


void CRF::write_model(char *model_file, bool first_part)
{
	/*
	version
	templates

	is_real

	ysize
	y

	cols

	xnum
	x

	*/

	FILE *fout;
	if(first_part)
	{
		fout=fopen(model_file,"wb");
		if(!fout)
		{
			printf("can not open model file: %s\n",model_file);
			return;
		}
		int i,j;
		//write version
		fwrite(&version, sizeof(int), 1, fout);
		//write templets
		i=templets.size();
		fwrite(&i, sizeof(int), 1, fout);
		char line[MAXSTRLEN];
		for(i=0;i<templets.size();i++)
		{
			templet &cur_templet=templets[i];
			line[0]=0;
			for(j=0;j<cur_templet.x.size();j++)
				sprintf(line,"%s%%x[%d,%d]",line,cur_templet.x[j].first,cur_templet.x[j].second);
			for(j=0;j<cur_templet.y.size();j++)
				sprintf(line,"%s%%y[%d]",line,cur_templet.y[j]);
			j=strlen(line);
			fwrite(&j, sizeof(int), 1, fout);
			fwrite(line, sizeof(char), j, fout);
		}
		
		if(is_real)
			i=1;
		else
			i=0;
		fwrite(&i, sizeof(int), 1, fout);
		//write y
		fwrite(&ysize, sizeof(int), 1, fout);
		for(i=0;i<tags.size();i++)
		{
			j=strlen(tags[i]);
			fwrite(&j, sizeof(int), 1, fout);
			fwrite(tags[i], sizeof(char), j, fout);
		}
		//write x
		fwrite(&cols, sizeof(int), 1, fout);
		i=xindex.size();
		fwrite(&i, sizeof(int), 1, fout);
		map<char*, int, str_cmp>::iterator it;
		for(it = xindex.begin(); it != xindex.end(); it++)
		{
			i=strlen(it->first);
			fwrite(&i, sizeof(int), 1, fout);
			fwrite(it->first, sizeof(char), i, fout);
			fwrite(&it->second, sizeof(int), 1, fout);
		}
		fclose(fout);
	}else{
		int i,j,n;
		
		//compress fmap
		//	index	fmap	ffmap[index]=new_index	new_fmap[new_index]=fmap[index]
		//	0		0		0						0
		//	1		-1		1						-1
		//	2		-2		-1
		//	3		-2		-1
		//	4		1		2						1
		//	5		2		3						2
		//
		//
		vector<int> ffmap(fmap_size);
		vector<int> new_fmap;
		int compress_fmap_index=0;
		bool is_all_valid=true;
		for(i=0;i<fmap_size;i++){
			if(fmap[i]!=-2){
				ffmap[i]=compress_fmap_index++;
				new_fmap.push_back(fmap[i]);
				if(fmap[i]==-1)
					is_all_valid=false;
			}else{
				ffmap[i]=-1;
			}
		}
		
		char line[MAXSTRLEN];
		FILE *fp;
		fp=fopen(model_file,"rb");
		fread(line,sizeof(int),1,fp);//skip version line
		//load templates
		fread(&n,sizeof(int),1,fp);
		for(i=0;i<n;i++){
			fread(&j,sizeof(int),1,fp);
			fread(line,sizeof(char),j,fp);
			line[j]=0;
			add_templet(line);
		}
		set_order();
		//get is_real
		fread(&i,sizeof(int),1,fp);
		if(i)
			is_real=true;
		else
			is_real=false;
		fread(&ysize,sizeof(int),1,fp);
		tags.resize(ysize);
		for(i=0;i<ysize;i++)
		{
			fread(&n,sizeof(int),1,fp);
			fread(line,sizeof(char),n,fp);
			line[n]=0;
			char *q=tag_str.push_back(line);
			tags[i]=q;
		}
		set_group();
		//get cols
		fread(&cols,sizeof(int),1,fp);
		//load x
		int index;
		int x_num;
		fread(&x_num,sizeof(int),1,fp);
		vector<char *> new_x;
		vector<int> new_index;
		for(i=0;i<x_num;i++)
		{
			fread(&n,sizeof(int),1,fp);
			fread(line,sizeof(char),n,fp);
			line[n]=0;
			fread(&index,sizeof(int),1,fp);
			char *q=x_str.push_back(line);
			index=ffmap[index];
			if(index!=-1){
				new_x.push_back(q);
				new_index.push_back(index);
			}
		}
		if(is_all_valid)//no feature removed, 100% dense
		{
			delete [] fmap;
			fmap=NULL;
			fmap_size=0;
		}else{
			memcpy(fmap,&new_fmap[0],sizeof(int)*new_fmap.size());
			fmap_size=new_fmap.size();
		}
		fclose(fp);

		//rewrite first part
		fout=fopen(model_file,"wb");
		if(!fout)
		{
			printf("can not open model file: %s\n",model_file);
			return;
		}
		//write version
		fwrite(&version, sizeof(int), 1, fout);
		//write templets
		i=templets.size();
		fwrite(&i, sizeof(int), 1, fout);
		
		for(i=0;i<templets.size();i++)
		{
			templet &cur_templet=templets[i];
			line[0]=0;
			for(j=0;j<cur_templet.x.size();j++)
				sprintf(line,"%s%%x[%d,%d]",line,cur_templet.x[j].first,cur_templet.x[j].second);
			for(j=0;j<cur_templet.y.size();j++)
				sprintf(line,"%s%%y[%d]",line,cur_templet.y[j]);
			j=strlen(line);
			fwrite(&j, sizeof(int), 1, fout);
			fwrite(line, sizeof(char), j, fout);
		}
		
		if(is_real)
			i=1;
		else
			i=0;
		fwrite(&i, sizeof(int), 1, fout);
		//write y
		fwrite(&ysize, sizeof(int), 1, fout);
		for(i=0;i<tags.size();i++)
		{
			j=strlen(tags[i]);
			fwrite(&j, sizeof(int), 1, fout);
			fwrite(tags[i], sizeof(char), j, fout);
		}
		//write x
		fwrite(&cols, sizeof(int), 1, fout);
		i=new_x.size();
		fwrite(&i, sizeof(int), 1, fout);
		for(i=0;i<new_x.size();i++)
		{
			n=strlen(new_x[i]);
			fwrite(&n, sizeof(int), 1, fout);
			fwrite(new_x[i], sizeof(char), n, fout);
			fwrite(&new_index[i], sizeof(int), 1, fout);
		}

		//write fmap
		fwrite(&fmap_size, sizeof(int), 1, fout);
		fwrite(fmap, sizeof(int), fmap_size, fout);
		//lambda
		fwrite(&lambda_size, sizeof(int), 1, fout);
		fwrite(lambda, sizeof(double), lambda_size, fout);
		//free all members' memory
		templet_group.clear();
		if(lambda)
		{
			delete [] lambda;
			lambda=NULL;
		}
		if(gradient)
		{
			delete [] gradient;
			gradient=NULL;
		}
		for(i=1;i<thread_num;i++)
		{
			delete [] threads[i].gradient;
			threads[i].gradient=NULL;
		}
		//build double array trie
		dat->build(new_x.size(),(unsigned char **)&new_x[0],&new_index[0]);
		dat->write(fout);
		fclose(fout);
	}
}




void CRF::shrink_feature()
{
	int i,j,k,ii;
	if(freq_thresh<=0)
		return;
    map<int, int> old2new;
    int new_lambda_size = 0;
	int new_mapping_size=0;
	char temp[MAXSTRLEN];
	map<char*, int, str_cmp>::iterator it;
	int findex=0;//base feature index
	for(i=0;i<lambda_size;i++)
	{
		if(i<x_freq.size() && x_freq[i]>0)
			findex=i;
		if(fmap_tmp[i]<freq_thresh)
		{
			x_freq[findex]-=fmap_tmp[i];
			fmap_tmp[i]=-1;//	mapping current feature index to null
		}else{
			fmap_tmp[i]=new_lambda_size++;// mapping current feature index to new_lambda_size 
		}
	}

	vector<int> new_fmap;
	lambda_size=new_lambda_size;
    for (it= xindex.begin(); it!= xindex.end();)
	{
		char *key=it->first;
		catch_string(key,":",temp);
		int index=atoi(temp);
		int gram_num=templets[index].y.size();
		if (x_freq[it->second] >0)
		{
			old2new.insert(make_pair<int, int>(it->second, new_mapping_size));
			new_fmap.insert(new_fmap.end(),fmap_tmp.begin()+it->second,fmap_tmp.begin()+it->second+pow(double(ysize),gram_num));
			it->second = new_mapping_size;
			new_mapping_size += pow(double(ysize),gram_num);
			++it;
		}else{
			xindex.erase(it++);
		}
    }
	fmap_tmp=new_fmap;
	new_fmap.clear();
	for(i=0;i<fmap_tmp.size();i++)
		if(fmap_tmp[i]==-1)
			fmap_tmp[i]=lambda_size;
	lambda_size++;

	map<int, int>::iterator iter;
	freelist<int> temp_clique_feature;
	freelist<double> temp_clique_fvalue;
	temp_clique_feature.set_size(PAGESIZE*16);
	temp_clique_fvalue.set_size(PAGESIZE*16);
	clique_feature.free();
	clique_fvalue.free();

	for(i=0;i<sequences_tmp.size();i++)
	{
		sequence &seq=sequences_tmp[i];
		for(j=0;j<seq.node_num;j++)
		{
			node &nod=seq.nodes[j];
			for(k=0;k<nod.clique_num;k++)
			{
				if(!nod.cliques[k])
					continue;
				clique &cli=*nod.cliques[k];
				vector<int> newf;
				vector<double> newfv;
				for(ii=0;ii<cli.feature_num;ii++)
				{
					iter = old2new.find(cli.fvector[ii]);
					if(iter != old2new.end()){
						newf.push_back(iter->second);
						if(is_real)
							newfv.push_back(cli.fvalue[ii]);
					}
				}
				int *f=NULL;
				double *fv=NULL;
				if(newf.size())
				{
					f=temp_clique_feature.push_back(&newf[0],newf.size());
					if(is_real)
						fv=temp_clique_fvalue.push_back(&newfv[0],newfv.size());
				}
				cli.fvector=f;
				cli.fvalue=fv;
				cli.feature_num=newf.size();
			}
		}
	}
	clique_feature.clear();
	for(i=0;i<sequences_tmp.size();i++)
	{
		sequence &seq=sequences_tmp[i];
		for(j=0;j<seq.node_num;j++)
		{
			node &nod=seq.nodes[j];
			for(k=0;k<nod.clique_num;k++)
			{
				if(!nod.cliques[k])
					continue;
				clique &cli=*nod.cliques[k];
				if(cli.feature_num)
				{
					cli.fvector=clique_feature.push_back(cli.fvector,cli.feature_num);
					if(is_real)
						cli.fvalue=clique_fvalue.push_back(cli.fvalue,cli.feature_num);
					else
						cli.fvalue=NULL;
				}else{
					cli.fvector=NULL;
					cli.fvalue=NULL;
				}
			}
		}
	}
}



void CRF::tag(vector<vector<vector<string> > > &ext_table, vector<vector<vector<double> > > &val_table, vector<vector<string> > &best_tag,vector<double> &sequencep, vector<vector<double> > &nodep, vector<int> &con_pos, vector<string> &con_tag)
{
	sequence seq;
	int i,j;
	generate_sequence(ext_table,val_table,seq);
	crf_thread &cur_thread=threads.front();
	vector<int> con_y(con_tag.size());
	for(i=0;i<con_tag.size();i++){
		int insert_pos;
		vector_search(tags,(char *)(con_tag[i].c_str()),con_y[i],insert_pos,str_cmp());
	}
	cur_thread.build_lattice(seq);
	cur_thread.viterbi(seq,con_pos,con_y);
	
	
	for(i=0;i<cur_thread.best_path.size();i++)
	{
		vector<string> cur_tag(seq.node_num);
		for(j=0;j<seq.node_num;j++)
			cur_tag[j]=tags[cur_thread.best_path[i][j]];
		best_tag.push_back(cur_tag);
	}
	if(margin||seqp)
	{
		double z;
		cur_thread.forward_backward(seq,z);
		if(margin)
			cur_thread.node_margin(seq,nodep,z);
		if(seqp)
		{
			sequencep.resize(cur_thread.best_path.size());
			for(i=0;i<cur_thread.best_path.size();i++)
			{
				cur_thread.assign_tag(seq,cur_thread.best_path[i]);
				double c=cur_thread.path_cost(seq);
				sequencep[i]=exp(c-z);
			}
		}
	}
}


void CRF::tag(vector<vector<vector<string> > >&ext_table, vector<vector<vector<double> > >&val_table, vector<vector<string> > &best_tag,vector<double> &sequencep, vector<vector<double> > &nodep)
{
	sequence seq;
	generate_sequence(ext_table,val_table,seq);
	crf_thread &cur_thread=threads.front();
	cur_thread.build_lattice(seq);
	cur_thread.viterbi(seq);
	
	int i,j;
	for(i=0;i<cur_thread.best_path.size();i++)
	{
		vector<string> cur_tag(seq.node_num);
		for(j=0;j<seq.node_num;j++)
			cur_tag[j]=tags[cur_thread.best_path[i][j]];
		best_tag.push_back(cur_tag);
	}
	if(margin||seqp)
	{
		double z;
		cur_thread.forward_backward(seq,z);
		if(margin)
			cur_thread.node_margin(seq,nodep,z);
		if(seqp)
		{
			sequencep.resize(cur_thread.best_path.size());
			for(i=0;i<cur_thread.best_path.size();i++)
			{
				cur_thread.assign_tag(seq,cur_thread.best_path[i]);
				double c=cur_thread.path_cost(seq);
				sequencep[i]=exp(c-z);
			}
		}
	}
}


void CRF::generate_sequence(vector<vector<vector<string> > >&ext_table, vector<vector<vector<double> > >&val_table, sequence &seq)
{
	int i,j,k;
	int xcols=cols-1;
	int rows=ext_table.size();
	char s[1024];
	char s1[1024];
	char s2[1024];
	node* nod=nodes.alloc(rows);
	seq.node_num=rows;
	seq.nodes=nod;
	
	for(i=0;i<rows;i++)
	{
		nod[i].key=0;//random initialize
		vector<clique*> clisp;//features that affect on current nodes
		vector<int> feature_vector;
		vector<double> fvalue_vector;
		for(j=0;j<templets.size();j++)
		{
			//get first y's offset
			templet &pat=templets[j];
			if(pat.y[0]+i<0)
				continue;
			double fval=1;
			if(pat.x.size()>0){//if has x
				int index1,index2;
				bool has_xstring=true;//false, if unit=""
				vector<int> xid(pat.x.size(),0);//xid=(0,0): 0 th units + 0 th units
				vector<int> xtop(pat.x.size(),0);
				//get xtop
				for(k=0;k<pat.x.size();k++)
				{
					index1=pat.x[k].first+i;
					index2=pat.x[k].second;
					if(index1<0)
					{
						xtop[k]=1;
					}else if(index1>=rows){
						xtop[k]=1;
					}else if(!ext_table[index1][index2].size()){//no string here
						xtop[k]=0;
						has_xstring=false;
						break;
					}else{
						xtop[k]=ext_table[index1][index2].size();
					}
				}
				if(has_xstring)
				{
					xid.back()=-1;
					while(1)
					{
						//increase and check whether stop
						for(k=xid.size()-1;k>=0 && xid[k]+1 == xtop[k];k--)
							xid[k]=0;
						if(k<0)
							break;//stop
						xid[k]++;
						sprintf(s, "%d", j);
						strcat(s,":");
						//get x
						for(k=0;k<pat.x.size();k++)
						{
							index1=pat.x[k].first+i;
							index2=pat.x[k].second;
							if(index1<0)
								strcpy(s1,table_head[-index1-1].c_str());
							else if(index1>=rows)
								strcpy(s1,table_tail[index1-rows].c_str());
							else{
								strcpy(s1,ext_table[index1][index2][xid[k]].c_str());
								if(is_real)
									fval*=val_table[index1][index2][xid[k]];
							}
							strcat(s,s1);
							strcat(s,"\001");
						}
						//x obtained, insert x
						int index;
						if(dat->search((unsigned char*)s,index)){
							feature_vector.push_back(index);
							if(is_real)
								fvalue_vector.push_back(fval);
						}
						fval=1;
					}
				}
			}else{//else , no x
				sprintf(s, "%d", j);
				strcat(s,":");
				//x obtained, insert x
				int index;
				if(dat->search((unsigned char*)s,index)){
					feature_vector.push_back(index);
					if(is_real)
						fvalue_vector.push_back(fval);
				}
				fval=1;
			}
			if(pat.end_of_group)
			{//creat new clique
				clique cli;
				vector<node*> ns;
				for(k=0;k<pat.y.size();k++)
					ns.push_back(nod+i+pat.y[k]);
				node ** np=clique_node.push_back(&ns[0],ns.size());
				cli.nodes=np;
				cli.node_num=ns.size();
				if(feature_vector.size())
				{
					int *f=clique_feature.push_back(&feature_vector[0],feature_vector.size());
					cli.fvector=f;
					if(is_real)
						cli.fvalue=clique_fvalue.push_back(&fvalue_vector[0],fvalue_vector.size());
					else
						cli.fvalue=NULL;
				}else{
					cli.fvector=NULL;
					cli.fvalue=NULL;
				}
				cli.feature_num=feature_vector.size();
				cli.groupid=templets[j].groupid;
				cli.key=0;//random initialize
				clique *new_clique=cliques.push_back(&cli,1);
				clisp.push_back(new_clique);
				feature_vector.clear();
				fvalue_vector.clear();
			}
		}
		//set node -> clique
		if(clisp.size())
			nod[i].cliques = node_clique.push_back(&clisp[0],clisp.size());
		else
			nod[i].cliques = NULL;
		nod[i].clique_num =clisp.size();
	}
	nodes.free();
	node_clique.free();
	cliques.free();
	clique_node.free();
	clique_feature.free();
}

bool CRF::load_model(char *model_file)
{
	int i,j,n;
	char line[MAXSTRLEN];
	FILE *fp=fopen(model_file,"rb");
	if(!fp)
	{
		printf("model file: %s not found\n",model_file);
		return false;
	}
	//check version
	fread(&version,sizeof(int),1,fp);
	printf("model version: 0.%d\n",version);
	//load templates
	fread(&n,sizeof(int),1,fp);
	for(i=0;i<n;i++){
		fread(&j,sizeof(int),1,fp);
		fread(line,sizeof(char),j,fp);
		line[j]=0;
		add_templet(line);
	}
	set_order();
	printf("template number: %d \n",templets.size());
	//get is_real
	fread(&n,sizeof(int),1,fp);
	if(n==1){
		is_real=true;
		threads[0]._cal_cost=&crf_thread::real_cost;
	}else{
		is_real=false;
		threads[0]._cal_cost=&crf_thread::bool_cost;
	}
	//get ysize
	fread(&ysize,sizeof(int),1,fp);
	tags.resize(ysize);
	for(i=0;i<ysize;i++){
		fread(&j, sizeof(int), 1, fp);
		fread(line,sizeof(char),j,fp);
		line[j]=0;
		char *q=tag_str.push_back(line);
		tags[i]=q;
	}
	printf("tags number: %d \n",ysize);
	set_group();
	//get cols
	fread(&cols,sizeof(int),1,fp);
	//load x
	int index;
	int x_num;
	fread(&x_num,sizeof(int),1,fp);
	for(i=0;i<x_num;i++)
	{
		fread(&n,sizeof(int),1,fp);
		fread(line,sizeof(char),n,fp);
		line[n]=0;
		fread(&index,sizeof(int),1,fp);
	}
	//load fmap
	fread(&fmap_size,sizeof(int),1,fp);
	if(fmap_size)
	{
		fmap=new int[fmap_size];
		for(i=0;i<fmap_size;i++)
			fread(&fmap[i],sizeof(int),1,fp);
	}else{
		fmap=NULL;
	}

	//load lambda
	fread(&lambda_size,sizeof(int),1,fp);
	//fake fmap lambda
	if(fmap)
	{
		lambda=new double[lambda_size+1];
		lambda[lambda_size]=0;
	}
	else
		lambda=new double[lambda_size];
	for(i=0;i<lambda_size;i++)
		fread(&lambda[i],sizeof(double),1,fp);
	//load dat
	dat->load(fp);
	if(fmap)
	{
		for(i=0;i<fmap_size;i++)
			if(fmap[i]==-1)
				fmap[i]=lambda_size;
		lambda_size++;
	}
	if(is_real)
		threads[0]._cal_cost=&crf_thread::real_cost;
	else
		threads[0]._cal_cost=&crf_thread::bool_cost;
	if(fmap)
		threads[0]._findex=&crf_thread::map_findex;
	else
		threads[0]._findex=&crf_thread::dir_findex;
	threads[0].fmap=fmap;
	printf("%d lambda loaded\n",lambda_size);
	fclose(fp);
	path_num=pow((double)ysize,order+1);
	node_anum=pow((double)ysize,order);//alpha(beta) number of each node
	head_offset=-log((double)ysize)*order;
	return true;
}

void CRF::unload()//unload the sequence data to file to set memory free
{
	FILE *fp=fopen("__data1","wb");
	fwrite(work_space,sizeof(char),work_size/2,fp);
	fclose(fp);
	fp=fopen("__data2","wb");
	fwrite(work_space+work_size/2,sizeof(char),work_size-work_size/2,fp);
	fclose(fp);
}

void CRF::load()//unload the sequence data to file to set memory free
{
	FILE *fp=fopen("__data1","rb");
	fread(work_space,sizeof(char),work_size/2,fp);
	fclose(fp);
	fp=fopen("__data2","rb");
	fread(work_space+work_size/2,sizeof(char),work_size-work_size/2,fp);
	fclose(fp);
}

void CRF::adjust_data()
{
	int i,j,k;
	if(!fmap)
	{
		fmap_size=lambda_size;
		fmap=new int[fmap_size];
		lambda_size=0;
		for(i=0;i<fmap_size;i++)
			if(lambda[i]!=0)
			{
				gradient[lambda_size]=lambda[i];//copy lambda
				fmap[i]=lambda_size++;
			}else
				fmap[i]=-1;
	}else{
		lambda_size=0;
		for(i=0;i<fmap_size;i++)
		{
			if(lambda[fmap[i]]!=0)
			{
				gradient[lambda_size]=lambda[fmap[i]];//copy lambda
				fmap[i]=lambda_size++;
			}
			else
				fmap[i]=-1;
		}
	}
	//adjust cliques
	if(chain_type==SIMPLE_CHAIN)
	{
		for(i=0;i<path_num && fmap[transit+i]<0;i++);
		if(i==path_num){
			for(i=0;i<path_num;i++)
				 fmap[transit+i]=-2;
		}
	}
	for(i=0;i<sequence_num;i++){
		if(chain_type==GENERAL_CHAIN){
			sequence &seq=sequences[i];
			for(j=0;j<seq.node_num;j++){
				node &nod=seq.nodes[j];
				for(int k=0;k<nod.clique_num;k++){
					if(!nod.cliques[k])
						continue;
					clique &cli=*nod.cliques[k];
					for(int ii=0;ii<cli.feature_num;){
						int jj;
						for(jj=0;jj<templet_group[cli.groupid].size() && fmap[cli.fvector[ii]+jj]<0;jj++);
						if(jj==templet_group[cli.groupid].size()){
							for(jj=0;jj<templet_group[cli.groupid].size();jj++)
								fmap[cli.fvector[ii]+jj]=-2;//ready to remove corresponding feature string
							for(jj=ii;jj<cli.feature_num-1;jj++)
							{
								cli.fvector[jj]=cli.fvector[jj+1];
								if(is_real)
									cli.fvalue[jj]=cli.fvalue[jj+1];
							}
							cli.feature_num--;
						}else{
							ii++;
						}
					}
					if(cli.feature_num==0)
					{
						cli.fvector=NULL;
						cli.fvalue=NULL;
					}
				}
			}
		}else if(chain_type==FIRST_CHAIN||chain_type==SIMPLE_CHAIN){
			sequence1 &seq1=sequence1s[i];
			for(j=0;j<seq1.vertex_num;j++){
				vertex &vtx=seq1.vertexes[j];
				for(k=0;k<vtx.feature_num;)
				{
					int ii;
					for(ii=0;ii<ysize && fmap[vtx.fvector[k]+ii]<0;ii++);
					if(ii==ysize){
						for(ii=0;ii<ysize;ii++)
							fmap[vtx.fvector[k]+ii]=-2;
						for(ii=k;ii<vtx.feature_num-1;ii++)
						{
							vtx.fvector[ii]=vtx.fvector[ii+1];
							if(is_real)
								vtx.fvalue[ii]=vtx.fvalue[ii+1];
						}
						vtx.feature_num--;
					}else{
						k++;
					}
				}
				if(vtx.feature_num==0)
				{
					vtx.fvector=NULL;
					vtx.fvalue=NULL;
				}
				if(j && chain_type==FIRST_CHAIN){
					edge &e=seq1.edges[j-1];
					for(k=0;k<e.feature_num;)
					{
						int ii;
						for(ii=0;ii<path_num && fmap[e.fvector[k]+ii]<0;ii++);//==-1 || ==-2
						if(ii==path_num){
							for(ii=0;ii<path_num;ii++)
								fmap[e.fvector[k]+ii]=-2;
							for(ii=k;ii<e.feature_num-1;ii++)
							{
								e.fvector[ii]=e.fvector[ii+1];
								if(is_real)
									e.fvalue[ii]=e.fvalue[ii+1];
							}
							e.feature_num--;
						}else{
							k++;
						}
					}
					if(e.feature_num==0)
					{
						e.fvector=NULL;
						e.fvalue=NULL;
					}
				}
			}
		}
	}
	memcpy(lambda,gradient,sizeof(double)*lambda_size);//copy lambda
}
