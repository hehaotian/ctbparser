#include "crf_thread.h"
#include "fun.h"


void crf_thread::build_lattice(sequence &seq)
{
	if(path.size()<c->path_num*seq.node_num)
		path.resize(c->path_num*seq.node_num);//r*f for all path ended with current node
	fill(path.begin(),path.end(),0);
	int i,j,k,ii,jj;
	for(i=0;i<seq.node_num;i++)
	{
		double *cur_path=&path[c->path_num*i];
		node &nod=seq.nodes[i];
		//calculate r*f for paths in cliques order
		for(j=0;j<nod.clique_num;j++)
		{
			if(!nod.cliques[j])
				continue;
			clique &cli=*nod.cliques[j];//get j th clique
			vector<vector<int> > &group=c->templet_group[cli.groupid];
			for(k=0;k<cli.feature_num;k++)	
				for(ii=0;ii<group.size();ii++)
					for(jj=0;jj<group[ii].size();jj++)
						cur_path[group[ii][jj]]+=cal_cost(c->lambda[findex(cli.fvector[k]+ii)],cli.fvalue,k);
		}
	}
}

void crf_thread::build_lattice(sequence1 &seq1)
{
	int i,j,k;
	if(path.size()<c->path_num*(seq1.vertex_num-1) + c->ysize*seq1.vertex_num)
		path.resize(c->path_num*(seq1.vertex_num-1) + c->ysize*seq1.vertex_num);
	if(c->chain_type==FIRST_CHAIN)
		fill(path.begin(),path.end(),0);
	else if(c->chain_type==SIMPLE_CHAIN && (c->algorithm==CRF_ALGORITHM||c->algorithm==L1_ALGORITHM))
		memcpy((double*)&path[0],c->transit_buff,sizeof(double)*path.size());
	else if(c->chain_type==SIMPLE_CHAIN && (c->algorithm==AP_ALGORITHM || c->algorithm==PA_ALGORITHM))
	{
		fill(path.begin(),path.end(),0);
		for(i=0;i<c->path_num;i++)
			for(j=0;j<seq1.vertex_num-1;j++)
				path[2*c->ysize+(c->ysize+c->path_num)*j+i]=c->lambda[findex(c->transit+i)];
	}
	for(i=0;i<seq1.vertex_num;i++)
	{

		int cur_path_offset=i?(c->path_num+ c->ysize)*i-c->path_num:0;
		vertex &vtx=seq1.vertexes[i];
		for(j=0;j<vtx.feature_num;j++)
		{
			for(k=0;k<c->ysize;k++)
				path[cur_path_offset+k]+=cal_cost(c->lambda[findex(vtx.fvector[j]+k)],vtx.fvalue,j);
		}
		if(i && c->chain_type==FIRST_CHAIN)
		{
			edge &e=seq1.edges[i-1];
			for(j=0;j<e.feature_num;j++)
				for(k=0;k<c->path_num;k++)
					path[cur_path_offset+c->ysize+k]+=cal_cost(c->lambda[findex(e.fvector[j]+k)],e.fvalue,j);
		}
	}
}
double crf_thread::path_cost(sequence &seq)
{
	int i;
	double cost=0;
	for(i=0;i<seq.node_num;i++)
		cost+=path[c->path_num*i+seq.nodes[i].key];
	return cost;
}

double crf_thread::path_cost(sequence1 &seq1)
{
	int i;
	double cost=path[seq1.vertexes[0].key];
	for(i=1;i<seq1.vertex_num;i++)
	{
		int offset=(c->path_num+c->ysize)*(i-1)+c->ysize;
		cost+=path[offset+seq1.vertexes[i].key%c->ysize];
		cost+=path[offset+c->ysize+seq1.vertexes[i].key];
	}
	return cost;
}

double crf_thread::seq_fx_gx(sequence &seq)
{
	double s;
	double z;
	build_lattice(seq);
	s=path_cost(seq);
	forward_backward(seq,z);
	calculate_gradient(seq,z);
	return z-s;
}
double crf_thread::seq_fx_gx(sequence1 &seq1)
{
	double s;
	double z;
	build_lattice(seq1);
	s=path_cost(seq1);
	forward_backward(seq1,z);
	calculate_gradient(seq1,z);
	return z-s;
}

void crf_thread::forward_backward(sequence &seq, double &z)
{
	int i,j,k,ii;
	//forward
	int alpha_size=c->node_anum*seq.node_num;
	if(alpha.size()<alpha_size)
	{
		alpha.resize(alpha_size);
		beta.resize(alpha_size);
		first_cal.resize(alpha_size);
	}
	fill(alpha.begin(),alpha.end(),0);
	fill(first_cal.begin(),first_cal.end(),1);
	for(i=0;i<seq.node_num;i++)
	{
		double *cur_path=&path[c->path_num*i];
		//cal alpha of current node
		if(i>0)
		{
			double *cur_alpha=&alpha[i*c->node_anum];
			const double *last_alpha=&alpha[(i-1)*c->node_anum];
			int *cur_first=&first_cal[i*c->node_anum];
			for(j=0;j<c->path_num;j++)
			{
				ii=j % c->node_anum;
				k=j / c->ysize;
				if(!cur_first[ii])
				{
					cur_alpha[ii]=log_sum_exp(last_alpha[k]+cur_path[j],cur_alpha[ii]);
				}else{
					cur_alpha[ii]=last_alpha[k]+cur_path[j];
					cur_first[ii]=0;
				}
			}
		}else{
			double *cur_alpha=&alpha[i*c->node_anum];
			int *cur_first=&first_cal[i*c->node_anum];
			for(j=0;j<c->path_num;j++)
			{
				ii=j % c->node_anum;
				if(!cur_first[ii])
				{
					cur_alpha[ii]=log_sum_exp(cur_path[j]+ c->head_offset ,cur_alpha[ii]);
				}else{
					cur_alpha[ii]=cur_path[j]+ c->head_offset ;
					cur_first[ii]=0;
				}
			}
		}
	}

	//backward
	fill(beta.begin(),beta.end(),0);
	fill(first_cal.begin(),first_cal.end(),true);
	vector<double> last_path(c->path_num,0);
	for(i=seq.node_num-1;i>=0;i--)
	{
		//calculate beta of last node
		if(i<seq.node_num-1)
		{
			double *cur_beta=&beta[i*c->node_anum];
			double *last_beta=&beta[(i+1)*c->node_anum];
			int *cur_first=&first_cal[i*c->node_anum];
			double *last_path=&path[c->path_num*(i+1)];
			for(j=0;j<c->path_num;j++)
			{
				k=j % c->node_anum;
				ii=j / c->ysize;
				if(!cur_first[ii])
				{
					cur_beta[ii]=log_sum_exp(last_beta[k]+last_path[j],cur_beta[ii]);
				}else{
					cur_beta[ii]=last_beta[k]+last_path[j];
					cur_first[ii]=0;
				}
			}
		}else{
			double *cur_beta=&beta[i*c->node_anum];
			for(j=0;j<c->node_anum;j++)
				cur_beta[j]=0;
		}
	}
	//calculate z(x)
	z=alpha[c->node_anum*(seq.node_num-1)];
	for(i=1;i<c->node_anum;i++)
		z=log_sum_exp(z, alpha[c->node_anum*(seq.node_num-1)+i]);
}

void crf_thread::forward_backward(sequence1 &seq1, double &z)
{
	int i,j,k,ii;
	//forward
	int alpha_size=c->ysize*seq1.vertex_num;
	if(alpha.size()<alpha_size)
	{
		alpha.resize(alpha_size);
		beta.resize(alpha_size);
		first_cal.resize(alpha_size);
	}
	fill(alpha.begin(),alpha.end(),0);
	fill(first_cal.begin(),first_cal.end(),1);
	int cur_path_offset=-1;
	int last_path_offset;
	for(i=0;i<seq1.vertex_num;i++)
	{
		last_path_offset=cur_path_offset;
		cur_path_offset=i?(c->path_num+ c->ysize)*i-c->path_num:0;
		double *cur_path=&path[cur_path_offset];
		//cal alpha of current node
		if(i>0)
		{
			double *cur_alpha=&alpha[i*c->ysize];
			const double *last_alpha=&alpha[(i-1)*c->ysize];
			int *cur_first=&first_cal[i*c->ysize];
			for(j=0;j<c->path_num;j++)
			{
				ii=j % c->ysize;
				k=j / c->ysize;
				if(!cur_first[ii])
				{
					cur_alpha[ii]=log_sum_exp(last_alpha[k]+cur_path[j+c->ysize],cur_alpha[ii]);
				}else{
					cur_alpha[ii]=last_alpha[k]+cur_path[j+c->ysize];
					cur_first[ii]=0;
				}
			}
			for(j=0;j<c->ysize;j++)
				cur_alpha[j]+=cur_path[j];
		}else{
			double *cur_alpha=&alpha[0];
			for(j=0;j<c->ysize;j++)
				cur_alpha[j]=cur_path[j];
		}
	}

	//backward
	fill(beta.begin(),beta.end(),0);
	fill(first_cal.begin(),first_cal.end(),1);
	vector<double> last_path(c->path_num,0);
	last_path_offset=-1;
	for(i=seq1.vertex_num-1;i>=0;i--)
	{
		cur_path_offset=i?(c->path_num+ c->ysize)*i-c->path_num:0;
		//calculate beta of last node
		if(i<seq1.vertex_num-1)
		{
			double *cur_beta=&beta[i*c->ysize];
			double *last_beta=&beta[(i+1)*c->ysize];
			int *cur_first=&first_cal[i*c->ysize];
			double *last_path=&path[last_path_offset];
			for(j=0;j<c->path_num;j++)
			{
				k=j % c->ysize;
				ii=j / c->ysize;
				if(!cur_first[ii]){
					cur_beta[ii]=log_sum_exp(last_beta[k]+last_path[j+c->ysize]+last_path[k],cur_beta[ii]);
				}else{
					cur_beta[ii]=last_beta[k]+last_path[j+c->ysize]+last_path[k];
					cur_first[ii]=0;
				}
			}
		}else{
			double *cur_beta=&beta[i*c->ysize];
			for(j=0;j<c->ysize;j++)
				cur_beta[j]=0;
		}
		last_path_offset=cur_path_offset;
	}
	//calculate z(x)
	z=alpha[c->ysize*(seq1.vertex_num-1)];
	for(i=1;i<c->ysize;i++)
		z=log_sum_exp(z, alpha[c->ysize*(seq1.vertex_num-1)+i]);
}


void crf_thread::calculate_gradient(sequence &seq, double &z)
{
	int i,j,k,ii,jj;
	if(!margin.size())
		margin.resize(c->path_num);
	for(i=0;i<seq.node_num;i++)
	{
		double *cur_path=&path[c->path_num*i];
		double *cur_beta=&beta[c->node_anum*i];
		node &nod=seq.nodes[i];
		fill(margin.begin(),margin.end(),0);
		if(i>0)
		{
			double *last_alpha=&alpha[c->node_anum*(i-1)];
			for(j=0;j<c->path_num;j++)
				margin[j]=exp(cur_path[j] + last_alpha[j / c->ysize] + cur_beta[j % c->node_anum] - z);
			for(j=0;j<nod.clique_num;j++)
			{
				if(!nod.cliques[j])
					continue;
				clique &cli=*nod.cliques[j];//get j th clique
				vector<vector<int> > &group=c->templet_group[cli.groupid];
				for(k=0;k<group.size();k++)
					for(ii=0;ii<group[k].size();ii++)
						for(jj=0;jj<cli.feature_num;jj++)
							gradient[findex(cli.fvector[jj]+k)]+=cal_cost(margin[group[k][ii]],cli.fvalue,jj);
				for(k=0;k<cli.feature_num;k++)
					gradient[findex(cli.fvector[k]+cli.key)]-=cal_cost(1,cli.fvalue,k);
			}
		}else{//first node
			for(j=0;j<c->path_num;j++)
				margin[j]=exp(cur_path[j] + c->head_offset + cur_beta[j % c->node_anum] - z);
			for(j=0;j<nod.clique_num;j++)
			{
				if(!nod.cliques[j])
					continue;
				clique &cli=*nod.cliques[j];//get j th clique
				vector<vector<int> > &group=c->templet_group[cli.groupid];
				for(k=0;k<group.size();k++)
					for(ii=0;ii<group[k].size();ii++)
						for(jj=0;jj<cli.feature_num;jj++)
								gradient[findex(cli.fvector[jj]+k)]+=cal_cost(margin[group[k][ii]],cli.fvalue,jj);
				for(k=0;k<cli.feature_num;k++)
					gradient[findex(cli.fvector[k]+cli.key)]-=cal_cost(1,cli.fvalue,k);
			}
		}
	}
}

void crf_thread::calculate_gradient(sequence1 &seq1, double &z)
{
	int i,j,k;
	double margin1;
	int cur_path_offset;
	for(i=0;i<seq1.vertex_num;i++)
	{
		cur_path_offset=i?(c->path_num+ c->ysize)*i-c->path_num:0;
		double *cur_path=&path[cur_path_offset];
		double *cur_beta=&beta[c->ysize*i];
		double *cur_alpha=&alpha[c->ysize*i];
		vertex &vtx=seq1.vertexes[i];
		
		for(j=0;j<c->ysize;j++)
		{
			margin1=exp(cur_alpha[j]+cur_beta[j]-z);
			for(k=0;k<vtx.feature_num;k++)
				gradient[findex(vtx.fvector[k]+j)]+=cal_cost(margin1,vtx.fvalue,k);
		}
		for(j=0;j<vtx.feature_num;j++)
			gradient[findex(vtx.fvector[j] + vtx.key%c->ysize)]-=cal_cost(1,vtx.fvalue,j);
		
		if(i>0)
		{
			double *last_alpha=&alpha[c->ysize*(i-1)];
			if(c->chain_type==FIRST_CHAIN)
			{
				edge &e=seq1.edges[i-1];
				for(j=0;j<c->path_num;j++)
				{
					margin1=exp(cur_path[j+c->ysize] + cur_path[j%c->ysize]+ last_alpha[j / c->ysize] + cur_beta[j % c->ysize] - z);
					for(k=0;k<e.feature_num;k++)
						gradient[findex(e.fvector[k]+j)]+=cal_cost(margin1,e.fvalue,k);
				}
				for(j=0;j<e.feature_num;j++)
					gradient[findex(e.fvector[j]+vtx.key)]-=cal_cost(1,e.fvalue,j);
			}else if(c->chain_type==SIMPLE_CHAIN){
				for(j=0;j<c->path_num;j++)
				{
					margin1=exp(cur_path[j+c->ysize] + cur_path[j%c->ysize]+ last_alpha[j / c->ysize] + cur_beta[j % c->ysize] - z);
					gradient[findex(c->transit+j)]+=margin1;
				}
				gradient[findex(c->transit+vtx.key)]--;
			}
		}
	}
}


void crf_thread::run()
{
	if(c->algorithm==CRF_ALGORITHM|| c->algorithm==L1_ALGORITHM)
	{
		memset(gradient,0,sizeof(double)*c->lambda_size);
		obj=0;
		if(c->chain_type==GENERAL_CHAIN)
		{
			for(int i = start_i; i < c->sequence_num; i += c->thread_num)
			{
				obj += seq_fx_gx (c->sequences[i]);
			}
		}else if(c->chain_type==FIRST_CHAIN || c->chain_type==SIMPLE_CHAIN){
			for(int i = start_i; i < c->sequence_num; i += c->thread_num)
			{
				obj += seq_fx_gx (c->sequence1s[i]);
			}
		}
		if(fmap)
		{
			gradient[c->lambda_size-1]=0;
		}
	}else if(c->algorithm==AP_ALGORITHM){
		obj=0;
		for(int i = 0; i < c->sequence_num; i ++)
		{			
			if(c->chain_type==GENERAL_CHAIN)
			{
				build_lattice(c->sequences[i]);
				viterbi(c->sequences[i]);
				ap_update(c->sequences[i]);
			}else if(c->chain_type==FIRST_CHAIN || c->chain_type==SIMPLE_CHAIN){
				build_lattice(c->sequence1s[i]);
				viterbi(c->sequence1s[i]);
				ap_update(c->sequence1s[i]);
			}
			times--;
			if(fmap)
			{
				c->lambda[c->lambda_size-1]=0;
				gradient[c->lambda_size-1]=0;
			}
		}
	}else if(c->algorithm==PA_ALGORITHM){
		obj=0;
		for(int i = 0; i < c->sequence_num; i ++)
		{			
			if(c->chain_type==GENERAL_CHAIN)
			{
				build_lattice(c->sequences[i]);
				viterbi(c->sequences[i]);
				pa_update(c->sequences[i]);
			}else if(c->chain_type==FIRST_CHAIN || c->chain_type==SIMPLE_CHAIN){
				build_lattice(c->sequence1s[i]);
				viterbi(c->sequence1s[i]);
				pa_update(c->sequence1s[i]);
			}
			times--;
			if(fmap)
			{
				c->lambda[c->lambda_size-1]=0;
				gradient[c->lambda_size-1]=0;
			}
			
		}
	}
}



void crf_thread::assign_tag(sequence &seq, vector<int> &node_tag)
{
	int i,j,k;
	for(i=0;i<seq.node_num;i++){
		seq.nodes[i].key=0;
		for(j=0;j<=c->order;j++){
			if(i+j>=c->order)
				seq.nodes[i].key=seq.nodes[i].key*c->ysize+node_tag[i+j-c->order];
		}
	}
	for(i=0;i<seq.node_num;i++)
	{
		node &nod=seq.nodes[i];
		for(j=0;j<nod.clique_num;j++)
		{
			if(!nod.cliques[j])
				continue;
			clique &cli=*(nod.cliques[j]);
			int key=0;
			for(k=0;k<cli.node_num;k++)
				key= key*c->ysize +cli.nodes[k]->key%c->ysize;
			cli.key=key;
		}
	}
}

void crf_thread::ap_update(sequence &seq)
{
	int i,j,k;
	vector<int> seq_key(seq.node_num);//preserve seq.nodes[i].key
	for(i=0;i<seq.node_num;i++){
		seq_key[i]=seq.nodes[i].key;
		seq.nodes[i].key=0;
		for(j=0;j<=c->order;j++){
			if(i+j>=c->order)
				seq.nodes[i].key=seq.nodes[i].key*c->ysize+bst_path[i+j-c->order];
		}
	}

	for(i=0;i<seq.node_num;i++)
	{
		bool is_err=false;
		node &nod=seq.nodes[i];
		for(j=0;j<nod.clique_num;j++)
		{
			if(!nod.cliques[j])
				continue;
			clique &cli=*(nod.cliques[j]);
			int key=0;
			for(k=0;k<cli.node_num;k++)
				key= key*c->ysize +cli.nodes[k]->key%c->ysize;
			//key = predict cli.key
			if(key!=cli.key)
			{
				is_err=true;
				for(k=0;k<cli.feature_num;k++)
				{
					c->lambda[findex(cli.fvector[k]+key)]-=cal_cost(1,cli.fvalue,k);
					c->gradient[findex(cli.fvector[k]+key)]-=times*cal_cost(1,cli.fvalue,k);
					c->lambda[findex(cli.fvector[k]+cli.key)]+=cal_cost(1,cli.fvalue,k);
					c->gradient[findex(cli.fvector[k]+cli.key)]+=times*cal_cost(1,cli.fvalue,k);
				}
			}
		}
		if(is_err)
			obj++;
	}
	//recover
	for(i=0;i<seq.node_num;i++)
		seq.nodes[i].key=seq_key[i];
}

void crf_thread::pa_update(sequence &seq)
{
	int i,j,k;
	vector<int> seq_key(seq.node_num);//preserve seq.nodes[i].key
	for(i=0;i<seq.node_num;i++){
		seq_key[i]=seq.nodes[i].key;
		seq.nodes[i].key=0;
		for(j=0;j<=c->order;j++){
			if(i+j>=c->order)
				seq.nodes[i].key=seq.nodes[i].key*c->ysize+bst_path[i+j-c->order];
		}
	}

	//calculate \delta f(x,y) = f(x,y) - f(x,y~)
	map<int,double> res;
	double loss=0;
	for(i=0;i<seq.node_num;i++)
	{
		bool is_err=false;
		node &nod=seq.nodes[i];
		for(j=0;j<nod.clique_num;j++)
		{
			if(!nod.cliques[j])
				continue;
			clique &cli=*(nod.cliques[j]);
			int key=0;
			for(k=0;k<cli.node_num;k++)
				key= key*c->ysize +cli.nodes[k]->key%c->ysize;
			//key = predict cli.key
			if(key!=cli.key)
			{
				is_err=true;
				for(k=0;k<cli.feature_num;k++)
				{
					res[findex(cli.fvector[k]+key)]+=cal_cost(1,cli.fvalue,k);
					res[findex(cli.fvector[k]+cli.key)]-=cal_cost(1,cli.fvalue,k);
				}
			}
		}
		if(is_err){
			obj++;
			loss++;
		}
	}
	double r=0;
	if(fmap)
		res[c->lambda_size-1]=0;
	map<int,double>::iterator it;
	for(it=res.begin();it!=res.end();it++)
	{
		loss+=c->lambda[it->first]*it->second;
		r+=it->second*it->second;
	}
	r=loss/r;
	if(r>c->sigma)
		r=c->sigma;
	for(it=res.begin();it!=res.end();it++)
	{
		c->lambda[it->first]-=r*it->second;
		c->gradient[it->first]-=r*it->second*times;
	}
	//recover
	for(i=0;i<seq.node_num;i++)
		seq.nodes[i].key=seq_key[i];
}

void crf_thread::ap_update(sequence1 &seq1)//bst_path key is node.y, seq1.key is nodes.y
{
	int i,j;
	for(i=0;i<seq1.vertex_num;i++){
		vertex &vtx=seq1.vertexes[i];
		int key=0;
		for(j=0;j<2;j++){
			if(i+j>=1)
				key=key*c->ysize+bst_path[i+j-1];
		}
		bool is_err=false;
		if(key!=vtx.key)
		{
			is_err=true;
			
			for(j=0;j<vtx.feature_num;j++)
			{
				c->lambda[findex(vtx.fvector[j]+key%c->ysize)]-=cal_cost(1,vtx.fvalue,j);
				c->gradient[findex(vtx.fvector[j]+key%c->ysize)]-=times*cal_cost(1,vtx.fvalue,j);
				c->lambda[findex(vtx.fvector[j]+vtx.key%c->ysize)]+=cal_cost(1,vtx.fvalue,j);
				c->gradient[findex(vtx.fvector[j]+vtx.key%c->ysize)]+=times*cal_cost(1,vtx.fvalue,j);
			}
			
			if(i)
			{
				if(c->chain_type==SIMPLE_CHAIN)
				{
					c->lambda[findex(c->transit+key)]--;
					c->gradient[findex(c->transit+key)]-=times;
					c->lambda[findex(c->transit+vtx.key)]++;
					c->gradient[findex(c->transit+vtx.key)]+=times;
				}else if(c->chain_type==FIRST_CHAIN){
					edge &e=seq1.edges[i-1];
					
					for(j=0;j<e.feature_num;j++)
					{
						
						c->lambda[findex(e.fvector[j]+key)]-=cal_cost(1,e.fvalue,j);
						c->gradient[findex(e.fvector[j]+key)]-=times*cal_cost(1,e.fvalue,j);
						c->lambda[findex(e.fvector[j]+vtx.key)]+=cal_cost(1,e.fvalue,j);
						c->gradient[findex(e.fvector[j]+vtx.key)]+=times*cal_cost(1,e.fvalue,j);
					}
					
				}
			}
		}
		if(is_err)
			obj++;
	}
}


void crf_thread::pa_update(sequence1 &seq1)//bst_path key is node.y, seq1.key is nodes.y
{
	int i,j;
	//calculate \delta f(x,y) = f(x,y) - f(x,y~)
	map<int,double> res;
	double loss=0;
	for(i=0;i<seq1.vertex_num;i++){
		vertex &vtx=seq1.vertexes[i];
		int key=0;
		for(j=0;j<2;j++){
			if(i+j>=1)
				key=key*c->ysize+bst_path[i+j-1];
		}
		bool is_err=false;
		if(key!=vtx.key)
		{
			is_err=true;			
			for(j=0;j<vtx.feature_num;j++)
			{
				res[findex(vtx.fvector[j]+key%c->ysize)]+=cal_cost(1,vtx.fvalue,j);
				res[findex(vtx.fvector[j]+vtx.key%c->ysize)]-=cal_cost(1,vtx.fvalue,j);
			}
			
			if(i)
			{
				if(c->chain_type==SIMPLE_CHAIN)
				{
					res[findex(c->transit+key)]++;
					res[findex(c->transit+vtx.key)]--;
				}else if(c->chain_type==FIRST_CHAIN){
					edge &e=seq1.edges[i-1];
					for(j=0;j<e.feature_num;j++)
					{
						res[findex(e.fvector[j]+key)]+=cal_cost(1,e.fvalue,j);
						res[findex(e.fvector[j]+vtx.key)]-=cal_cost(1,e.fvalue,j);
					}	
				}
			}
		}
		if(is_err){
			obj++;
			loss++;
		}
	}
	double r=0;
	if(fmap)
		res[c->lambda_size-1]=0;
	map<int,double>::iterator it;
	for(it=res.begin();it!=res.end();it++)
	{
		loss+=c->lambda[it->first]*it->second;
		r+=it->second*it->second;
	}
	r=loss/r;
	if(r>c->sigma)
		r=c->sigma;
	for(it=res.begin();it!=res.end();it++)
	{
		c->lambda[it->first]-=r*it->second;
		c->gradient[it->first]-=r*it->second*times;
	}
}

void crf_thread::node_margin(sequence &seq, vector<vector<double> >&node_p, double &z)
{
	int i,j;
	node_p.resize(seq.node_num);
	for(i=0;i<seq.node_num;i++)
		node_p[i].resize(c->ysize);
	if(c->order>0)
	{
		vector<int> first_cal(c->ysize*seq.node_num,1);
		for(i=0;i<seq.node_num;i++)
		{	
			int *cur_first=&first_cal[c->ysize*i];
			double *cur_p=&node_p[i][0];
			for(j=0;j<c->node_anum;j++)
			{
				int index = j % c->ysize;
				if(cur_first[index])
				{
					cur_first[index]=0;
					cur_p[index]=alpha[i*c->node_anum+j]+beta[i*c->node_anum+j];
				}else
					cur_p[index]=log_sum_exp(alpha[i*c->node_anum+j]+beta[i*c->node_anum+j],cur_p[index]);
			}
			for(j=0;j<c->ysize;j++)
				cur_p[j]=exp(cur_p[j]-z);
		}
	}else{
		for(i=0;i<seq.node_num;i++)
		{
			double sum=path[c->path_num*i];
			for(j=1;j<c->ysize;j++)
				sum=log_sum_exp(path[c->path_num*i+j],sum);
			for(j=0;j<c->ysize;j++)
				node_p[i][j]=exp(path[c->path_num*i+j]-sum);
		}
	}
}


void crf_thread::viterbi(sequence &seq)
{
	int i,j;
	int i1,j1,k1,i2,j2,k2;
	//last node: i1 th node, j1 th tag, k1 th best
	//current node: i2 th node, j2 th tag, k2 th best
	if(final_path.size()==0)
		final_path.resize(c->node_anum);
	fill(final_path.begin(),final_path.end(),0);
	last_best.clear();
	last_best.resize(c->node_anum);
	// i th tag, j th best: last_best[i][j]
	cur_best.clear();
	cur_best.resize(c->node_anum);
	best_prev.clear();
	best_prev.resize(seq.node_num+1);
	//best_prev of i th node j th tag, k th best: best_prev[i][j][k]
	best_link.clear();
	best_link.resize(seq.node_num+1);
	//index of the path from i th node j th tag, k th best to i-1 th node best_prev[i][j][k].first th tag, best_prev[i][j][k].second th best
	bst_path.clear();
	bst_path.resize(seq.node_num);
	best_path.clear();
	for(i=0;i<=seq.node_num;i++)
	{
		i2=i;
		//current node: i2 th node
		i1=i-1;
		//last node: i1 th node
		double *cur_path;
		if(i<seq.node_num)
		{
			best_prev[i].resize(c->node_anum);
			best_link[i].resize(c->path_num);
			cur_path=&path[c->path_num*i];
		}else{
			best_prev[i].resize(1);
			best_link[i].resize(1);
			cur_path=&final_path[0];
		}

		//search n-best path
		if(i>0)
		{
			last_best=cur_best;
			for(j=0;j<c->node_anum;j++)
				cur_best[j].clear();
			if(i==seq.node_num)
				cur_best.resize(1);
		}else{
			vector<double> init_best(1,0);
			last_best[0]=init_best;
		}
		int routine_num=i<seq.node_num?c->path_num:c->node_anum;
		for(j=0;j<routine_num;j++)
		{
			if(i<seq.node_num)
			{	
				j2 = j % c->node_anum;//current node, j2 th tag
				j1 = j / c->ysize;//last node, j1 th tag
			}else{
				j2 = 0;
				j1 = j;
			}
			for(k1=0;k1<last_best[j1].size();k1++)
			{

				double cost=last_best[j1][k1];
				cost+=cur_path[j];
				
				//last node : j1 th tag,  k1 th best
				int k_pos;
				vector_search(cur_best[j2],cost,k_pos,k2,inverse_cmp<double>());
				//current node : j2 th tag, k2 th best, now, search for k2
				if(k2<c->nbest)
				{
					vector_insert(cur_best[j2],cost,k2);
					pair<int,int> q=make_pair(j1,k1);
					vector_insert(best_prev[i2][j2], q ,k2);
					vector_insert(best_link[i2][j2], j, k2);
					//for current node, i.e. i2 th node, its j2 th tag, k2 th 
					//best previous is the j1 th tag, k1 th best of last node
					if(cur_best[j2].size()>c->nbest)//drop the nbest+1 th candidate
					{
						cur_best[j2].pop_back();
						best_prev[i2][j2].pop_back();
						best_link[i2][j2].pop_back();
					}
				}else{
					break;
				}
			}
		}
	}
	for(i=0;i<best_prev[seq.node_num][0].size();i++)
	{
		k2=i;
		j2=0;
		for(i2=seq.node_num;i2>0;i2--)
		{
			
			pair<int ,int> p=best_prev[i2][j2][k2];
			j2=p.first;
			k2=p.second;
			bst_path[i2-1]=best_link[i2-1][j2][k2] % c->ysize;
		}
		best_path.push_back(bst_path);
	}
}

void crf_thread::viterbi(sequence &seq, vector<int> &con_pos, vector<int> &con_y)
{
	int i,j;
	int i1,j1,k1,i2,j2,k2;
	//last node: i1 th node, j1 th tag, k1 th best
	//current node: i2 th node, j2 th tag, k2 th best
	if(final_path.size()==0)
		final_path.resize(c->node_anum);
	fill(final_path.begin(),final_path.end(),0);
	last_best.clear();
	last_best.resize(c->node_anum);
	// i th tag, j th best: last_best[i][j]
	cur_best.clear();
	cur_best.resize(c->node_anum);
	best_prev.clear();
	best_prev.resize(seq.node_num+1);
	//best_prev of i th node j th tag, k th best: best_prev[i][j][k]
	best_link.clear();
	best_link.resize(seq.node_num+1);
	//index of the path from i th node j th tag, k th best to i-1 th node best_prev[i][j][k].first th tag, best_prev[i][j][k].second th best
	bst_path.clear();
	bst_path.resize(seq.node_num);
	best_path.clear();
	int con_i=0;
	for(i=0;i<=seq.node_num;i++)
	{
		int cur_con_y=-1;
		if(con_i<con_pos.size() && i==con_pos[con_i])
		{
			cur_con_y=con_y[con_i];
			con_i++;
		}

		i2=i;
		//current node: i2 th node
		i1=i-1;
		//last node: i1 th node
		double *cur_path;
		if(i<seq.node_num)
		{
			best_prev[i].resize(c->node_anum);
			best_link[i].resize(c->path_num);
			cur_path=&path[c->path_num*i];
		}else{
			best_prev[i].resize(1);
			best_link[i].resize(1);
			cur_path=&final_path[0];
		}

		//search n-best path
		if(i>0)
		{
			last_best=cur_best;
			for(j=0;j<c->node_anum;j++)
				cur_best[j].clear();
			if(i==seq.node_num)
				cur_best.resize(1);
		}else{
			vector<double> init_best(1,0);
			last_best[0]=init_best;
		}
		int routine_num=i<seq.node_num?c->path_num:c->node_anum;
		for(j=0;j<routine_num;j++)
		{
			if(i<seq.node_num)
			{
				j2 = j % c->node_anum;//current node, j2 th tag
				j1 = j / c->ysize;//last node, j1 th tag
			}else{
				j2 = 0;
				j1 = j;
			}
			if(cur_con_y!=-1 && j2!=cur_con_y)
				continue;
			for(k1=0;k1<last_best[j1].size();k1++)
			{

				double cost=last_best[j1][k1];
				cost+=cur_path[j];
				
				//last node : j1 th tag,  k1 th best
				int k_pos;
				vector_search(cur_best[j2],cost,k_pos,k2,inverse_cmp<double>());
				//current node : j2 th tag, k2 th best, now, search for k2
				if(k2<c->nbest)
				{
					vector_insert(cur_best[j2],cost,k2);
					pair<int,int> q=make_pair(j1,k1);
					vector_insert(best_prev[i2][j2], q ,k2);
					vector_insert(best_link[i2][j2], j, k2);
					//for current node, i.e. i2 th node, its j2 th tag, k2 th 
					//best previous is the j1 th tag, k1 th best of last node
					if(cur_best[j2].size()>c->nbest)//drop the nbest+1 th candidate
					{
						cur_best[j2].pop_back();
						best_prev[i2][j2].pop_back();
						best_link[i2][j2].pop_back();
					}
				}else{
					break;
				}
			}
		}
	}
	for(i=0;i<best_prev[seq.node_num][0].size();i++)
	{
		k2=i;
		j2=0;
		for(i2=seq.node_num;i2>0;i2--)
		{
			
			pair<int ,int> p=best_prev[i2][j2][k2];
			j2=p.first;
			k2=p.second;
			bst_path[i2-1]=best_link[i2-1][j2][k2] % c->ysize;
		}
		best_path.push_back(bst_path);
	}
}


void crf_thread::viterbi(sequence1 &seq1)
{
	int i,j;
	int i1,j1,k1,i2,j2,k2;
	//last node: i1 th node, j1 th tag, k1 th best
	//current node: i2 th node, j2 th tag, k2 th best
	if(final_path.size()==0)//initialize
		final_path.resize(c->ysize);
	fill(final_path.begin(),final_path.end(),0);
	last_best.clear();
	cur_best.clear();
	last_best.resize(c->ysize);
	cur_best.resize(c->ysize);
	// i th tag, j th best: last_best[i][j]
	best_prev.clear();
	best_prev.resize(seq1.vertex_num+1);
	//best_prev of i th node j th tag, k th best: best_prev[i][j][k]
	best_link.clear();
	best_link.resize(seq1.vertex_num+1);
	//index of the path from i th node j th tag, k th best to i-1 th node best_prev[i][j][k].first th tag, best_prev[i][j][k].second th best
	bst_path.clear();
	bst_path.resize(seq1.vertex_num);
	best_path.clear();
	for(i=0;i<=seq1.vertex_num;i++)
	{
		int cur_path_offset=i?(c->path_num+ c->ysize)*i-c->path_num:0;
		i2=i;
		//current node: i2 th node
		i1=i-1;
		//last node: i1 th node
		double *cur_path;
		if(i<seq1.vertex_num)
		{
			best_prev[i].resize(c->ysize);
			best_link[i].resize(c->path_num);
			cur_path=&path[cur_path_offset];
		}else{
			best_prev[i].resize(1);
			best_link[i].resize(1);
			cur_path=&final_path[0];
		}

		//search n-best path
		if(i>0)
		{
			last_best=cur_best;
			for(j=0;j<c->ysize;j++)
				cur_best[j].clear();
			if(i==seq1.vertex_num)
				cur_best.resize(1);
		}else{
			vector<double> init_best(1,0);
			last_best[0]=init_best;
		}
		int routine_num=i<seq1.vertex_num?c->path_num:c->ysize;
		for(j=0;j<routine_num;j++)
		{
			if(i<seq1.vertex_num)
			{	
				j2 = j % c->ysize;//current node, j2 th tag
				j1 = j / c->ysize;//last node, j1 th tag
			}else{
				j2 = 0;
				j1 = j;
			}
			for(k1=0;k1<last_best[j1].size();k1++)
			{
				double cost=last_best[j1][k1];
				if(i && i<seq1.vertex_num)
					cost+=cur_path[j+c->ysize];//edge weight
				cost+=cur_path[j%c->ysize];//node weight
				//last node : j1 th tag,  k1 th best
				int k_pos;
				vector_search(cur_best[j2],cost,k_pos,k2,inverse_cmp<double>());
				//current node : j2 th tag, k2 th best, now, search for k2
				if(k2<c->nbest)
				{
					vector_insert(cur_best[j2],cost,k2);
					pair<int,int> q=make_pair(j1,k1);
					vector_insert(best_prev[i2][j2], q ,k2);
					vector_insert(best_link[i2][j2], j, k2);
					//for current node, i.e. i2 th node, its j2 th tag, k2 th 
					//best previous is the j1 th tag, k1 th best of last node
					if(cur_best[j2].size()>c->nbest)//drop the nbest+1 th candidate
					{
						cur_best[j2].pop_back();
						best_prev[i2][j2].pop_back();
						best_link[i2][j2].pop_back();
					}
				}else{
					break;
				}
			}
		}
	}
	for(i=0;i<best_prev[seq1.vertex_num][0].size();i++)
	{
		k2=i;
		j2=0;
		for(i2=seq1.vertex_num;i2>0;i2--)
		{
			
			pair<int ,int> p=best_prev[i2][j2][k2];
			j2=p.first;
			k2=p.second;
			bst_path[i2-1]=best_link[i2-1][j2][k2] % c->ysize;
		}
		best_path.push_back(bst_path);
	}
}

