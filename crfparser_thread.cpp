#include "crfparser_thread.h"
#include "fun.h"
#include "const.h"
#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
void crfparser_thread::run(){
	
	int i,j;
	if(_c->_algorithm==CRF_ALGORITHM||_c->_algorithm==L1_ALGORITHM){
		for(i = _start_i; i < _c->_sentences.size(); i += _c->_thread_num)
			_obj += sen_fx_gx (_c->_sentences[i]);
	}else if(_c->_algorithm==AP_ALGORITHM){
		for(i = 0; i < _c->_sentences.size(); i ++){
			build_lattice_max(_c->_sentences[i]);
			cyk(_c->_sentences[i]);
			ap_update(_c->_sentences[i]);
			_times--;
		}
	}else if(_c->_algorithm==PA_ALGORITHM){
		
		for(i = 0; i < _c->_sentences.size(); i ++){
			build_lattice_max(_c->_sentences[i]);
			cyk(_c->_sentences[i]);
			pa_update(_c->_sentences[i]);
			_times--;
		}
	}else if(_c->_algorithm==FACTOR_ALGORITHM){
		
		int ysize=base_feature::_labels.size();
		for(i = 0; i < _c->_sentences.size(); i ++){
			double s,z;
			sentence &sen=_c->_sentences[i];
			build_lattice(sen);
			//scale
			for(j=0;j<sen.l*sen.l;j++){
				_lattice[j]*=_c->_factor;
			}
			for(j=0;j<sen.l*sen.l*ysize;j++){
				_edgep[j]*=_c->_factor;
			}
			collapse_lattice(sen);
			s=tree_cost(sen);
			inside_outside(sen,z);
			calculate_margin(sen,z);
			calculate_gradient_factor(sen,z);
			_c->_factor_gradient-=s;
			_obj += z-s;
		}
	}
}

void crfparser_thread::ap_update(sentence &sen){
	int i,j,k;
	int l=sen.l;
	int grid;
	int ysize=base_feature::_labels.size();
	for(i=0;i<l;i++){
		if(sen.parent[i]!=_par[0][i] || sen.label[i]!=_label[0][i]){
			_obj++;
			if(sen.parent[i]!=l){
				grid=get_feature_grid(l,sen.parent[i],i,true,true);
				k=sen.label[i];
				for(j=0;j<sen.fnum[grid][k];j++){
					_c->_lambda[sen.fvector[grid][k][j]]++;
					_c->_gradient[sen.fvector[grid][k][j]]+=_times;
				}
				grid=get_feature_grid(l,sen.parent[i],i,false,true);
				for(j=0;j<sen.fnum[grid][k];j++){
					_c->_lambda[sen.fvector[grid][k][j]]++;
					_c->_gradient[sen.fvector[grid][k][j]]+=_times;
				}
				grid=get_feature_grid(l,sen.parent[i],i,true,false);
				for(j=0;j<sen.fnum[grid][k];j++){
					_c->_lambda[sen.fvector[grid][k][j]]++;
					_c->_gradient[sen.fvector[grid][k][j]]+=_times;
				}
			}else{
				grid=get_feature_grid(l,sen.parent[i],i,false,false);
				for(j=0;j<sen.fnum[grid][0];j++){
					_c->_lambda[sen.fvector[grid][0][j]]++;
					_c->_gradient[sen.fvector[grid][0][j]]+=_times;
				}
			}

			if(_par[0][i]!=l){
				grid=get_feature_grid(l,_par[0][i],i,true,true);
				int k=_label[0][i];
				for(j=0;j<sen.fnum[grid][k];j++){
					_c->_lambda[sen.fvector[grid][k][j]]--;
					_c->_gradient[sen.fvector[grid][k][j]]-=_times;
				}
				grid=get_feature_grid(l,_par[0][i],i,false,true);
				for(j=0;j<sen.fnum[grid][k];j++){
					_c->_lambda[sen.fvector[grid][k][j]]--;
					_c->_gradient[sen.fvector[grid][k][j]]-=_times;
				}
				grid=get_feature_grid(l,_par[0][i],i,true,false);
				for(j=0;j<sen.fnum[grid][k];j++){
					_c->_lambda[sen.fvector[grid][k][j]]--;
					_c->_gradient[sen.fvector[grid][k][j]]-=_times;
				}
			}else{
				grid=get_feature_grid(l,_par[0][i],i,false,false);
				for(j=0;j<sen.fnum[grid][0];j++){
					_c->_lambda[sen.fvector[grid][0][j]]--;
					_c->_gradient[sen.fvector[grid][0][j]]-=_times;
				}
			}
		}
	}
}

void crfparser_thread::pa_update(sentence &sen){
	int i,j,k;
	int l=sen.l;
	int grid;
	int ysize=base_feature::_labels.size();
	//calculate \delta f(x,y) = f(x,y) - f(x,y~)
	map<int,double> res;
	double loss=0;
	for(i=0;i<l;i++){
		if(sen.parent[i]!=_par[0][i] || sen.label[i]!=_label[0][i]){
			_obj++;
			loss++;
			if(sen.parent[i]!=l){
				grid=get_feature_grid(l,sen.parent[i],i,true,true);
				k=sen.label[i];
				for(j=0;j<sen.fnum[grid][k];j++){
					res[sen.fvector[grid][k][j]]--;
				}
				grid=get_feature_grid(l,sen.parent[i],i,false,true);
				for(j=0;j<sen.fnum[grid][k];j++){
					res[sen.fvector[grid][k][j]]--;
				}
				grid=get_feature_grid(l,sen.parent[i],i,true,false);
				for(j=0;j<sen.fnum[grid][k];j++){
					res[sen.fvector[grid][k][j]]--;
				}
			}else{
				grid=get_feature_grid(l,sen.parent[i],i,false,false);
				for(j=0;j<sen.fnum[grid][0];j++){
					res[sen.fvector[grid][0][j]]--;
				}
			}

			if(_par[0][i]!=l){
				grid=get_feature_grid(l,_par[0][i],i,true,true);
				int k=_label[0][i];
				for(j=0;j<sen.fnum[grid][k];j++){
					res[sen.fvector[grid][k][j]]++;
				}
				grid=get_feature_grid(l,_par[0][i],i,false,true);
				for(j=0;j<sen.fnum[grid][k];j++){
					res[sen.fvector[grid][k][j]]++;
				}
				grid=get_feature_grid(l,_par[0][i],i,true,false);
				for(j=0;j<sen.fnum[grid][k];j++){
					res[sen.fvector[grid][k][j]]++;
				}
			}else{
				grid=get_feature_grid(l,_par[0][i],i,false,false);
				for(j=0;j<sen.fnum[grid][0];j++){
					res[sen.fvector[grid][0][j]]++;
				}
			}
		}
	}
	double r=0;
	map<int,double>::iterator it;
	for(it=res.begin();it!=res.end();it++)
	{
		loss+=_c->_lambda[it->first]*it->second;
		r+=it->second*it->second;
	}
	r=loss/r;
	if(r>_c->_sigma)
		r=_c->_sigma;
	for(it=res.begin();it!=res.end();it++)
	{
		_c->_lambda[it->first]-=r*it->second;
		_c->_gradient[it->first]-=r*it->second*_times;
	}
}


double crfparser_thread::sen_fx_gx(sentence &sen){
	double s;
	double z;
	build_lattice(sen);
	collapse_lattice(sen);
	s=tree_cost(sen);
	inside_outside(sen,z);
	calculate_margin(sen,z);
	calculate_gradient(sen,z);

	return z-s;
}

void crfparser_thread::inside_outside(sentence &sen,double &z){
/*refrence:
Mark A.Paskin, Cubic-time Parsing and Learning Algorithm for Grammatical Bigram Models, Technical Report, 2001

chart cell order:
						bL bR s
0						F  F  (T, if |i-j|=1; F else)
1						F  T  F
2						T  F  F
3						F  T  T (total)
4						T  F  T (total)

chart[i*l+j]=edge(i,j), i<j
             edge(i,root), i==j

JOINT(i,k) (k,j)

if(k==i+1)
	0+2->0
	0+4->0
3+0->0
3+1->1
3+3->1
4+2->2
4+4->2

JOINT(i,j) (j,n)
4+4->2
4+2->2
3+0->0 + close
*/
	int ysize=base_feature::_labels.size();

	int l=sen.l;
	
	int cell_num=l*(l+1)/2;
	int i,j,k;
	int index,index1,index2,index3,index4,index5;
	int n;
	if(_alpha.size()<cell_num*5){
		_alpha.resize(cell_num*5);
		_has_alpha.resize(cell_num*5);
		_beta.resize(cell_num*5);
		_has_beta.resize(cell_num*5);
	}
	fill(_alpha.begin(),_alpha.end(),0);
	fill(_beta.begin(),_beta.end(),0);
	fill(_has_alpha.begin(),_has_alpha.end(),false);
	fill(_has_beta.begin(),_has_beta.end(),false);
	for(i=0;i<l;i++){
		index=cell_index(i,i+1,l)*5;
		//ADD(SEED(i))
		_has_alpha[index]=true;
		//CLOSE LEFT, i->i+1
		close_span(i,i+1,l);
	}
	
	
	for(n=2;n<=l;n++){//n=span len
		for(i=0;i<=l-n;i++){
			j=i+n;
			//(i,i+1) (i+1,j): 0+2->0, 0+4->0
			index3=cell_index(i,j,l)*5;
			if(j<l){// if j==l, F F + T F is not allowed, because of single root
				index1=cell_index(i,i+1,l)*5;
				index2=cell_index(i+1,j,l)*5;
				join_span(index1,index2+2,index3);
				join_span(index1,index2+4,index3);
			}
			for(k=i+1;k<j;k++){
				index1=cell_index(i,k,l)*5;
				index2=cell_index(k,j,l)*5;
				//3+0->0
				join_span(index1+3,index2,index3);
				//3+1->1
				if(j<l)
					join_span(index1+3,index2+1,index3+1);
				//3+3->1
				if(j<l)
					join_span(index1+3,index2+3,index3+1);
				//4+2->2
				join_span(index1+4,index2+2,index3+2);
				//4+4->2
				join_span(index1+4,index2+4,index3+2);
			}
			close_span(i,j,l);
		}
	}
	z=log_sum_exp(_alpha[2],_alpha[4]);

/*
chart cell order:
						bL bR s
0						F  F  (T, if |i-j|=1; F else)
1						F  T  F
2						T  F  F
3						F  T  T (total)
4						T  F  T (total)

chart[i*l+j]=edge(i,j), i<j
             edge(i,root), i==j

JOINT(i,k) (k,j)

if(k==i+1)
	0+2->0
	0+4->0
3+0->0
3+1->1
3+3->1
4+2->2
4+4->2

JOINT(i,j) (j,n)
4+4->2
4+2->2
3+0->0
*/
	
	index=cell_index(0,l,l)*5;
	if(l>1)
		_has_beta[index+2]=true;
	_has_beta[index+4]=true;
	for(n=l;n>=1;n--){
		
		for(i=0;i<=l-n;i++){
			j=i+n;
//			if(i==1 && j==3)
//				int xx=1;
			disclose_span(i,j,l);
			if(n==1)
				continue;
			//(i,i+1) (i+1,j): 0+2->0, 0+4->0
			index5=cell_index(i,j,l)*5;
			if(j<l){
				index1=cell_index(i,i+1,l)*5;
				index2=cell_index(i+1,j,l)*5;
				index3=cell_index(i,i+1,l)*5;
				index4=cell_index(i+1,j,l)*5;
				split_span(index1,index2+2,index3,index4+2,index5);
				split_span(index1,index2+4,index3,index4+4,index5);
			}
			for(k=i+1;k<j;k++){
				index1=cell_index(i,k,l)*5;
				index2=cell_index(k,j,l)*5;
				index3=cell_index(i,k,l)*5;
				index4=cell_index(k,j,l)*5;

				//3+0->0
				split_span(index1+3,index2,index3+3,index4,index5);
				//3+1->1
				if(j<l)
					split_span(index1+3,index2+1,index3+3,index4+1,index5+1);
				//3+3->1
				if(j<l)
					split_span(index1+3,index2+3,index3+3,index4+3,index5+1);
				//4+2->2
				split_span(index1+4,index2+2,index3+4,index4+2,index5+2);
				//4+4->2
				split_span(index1+4,index2+4,index3+4,index4+4,index5+2);
			}
		}
	}
}

void crfparser_thread::join_span(int index1,int index2,int index3){
	if(_has_alpha[index1] && _has_alpha[index2]){
		if(_has_alpha[index3]){
			_alpha[index3]=log_sum_exp(_alpha[index3],_alpha[index1]+_alpha[index2]);
		}else{
			_has_alpha[index3]=true;
			_alpha[index3]=_alpha[index1]+_alpha[index2];
		}
	}
}


bool crfparser_thread::merge_span_nbest(int pos1,int pos2,int pos3,int l){
	pair<vector< span >::const_iterator,vector< span >::const_iterator> ip;
	int ysize=base_feature::_labels.size();
	int nbest=_c->_nbest;
	int index1,index2,index3;
	index1=cell_index(pos1,pos2,l)*5*nbest;
	index2=cell_index(pos2,pos3,l)*5*nbest;
	index3=cell_index(pos1,pos3,l)*5*nbest;
	int k1,k2;
	int pos;
	int k;
	span t;
	t.y=-1;
	t.head=pos1;
	t.tail=pos3;
	if(pos2==pos1+1 && pos3<l){
		//(i,i+1) (i+1,j): 0+2->0
		assert(_spans[index1].y==-1);
		for(k2=0;k2<nbest && _spans[index2+2*nbest+k2].y!=-2;k2++){
			t.score=_spans[index1].score+_spans[index2+2*nbest+k2].score;
			ip = equal_range(_spans.begin()+index3, _spans.begin()+index3+nbest, t, span_cmp());
			pos=ip.second-(_spans.begin()+index3);
			if(pos>=nbest) break;
			//insert
			for(k=index3+nbest-1;k>index3+pos;k--)
				_spans[k]=_spans[k-1];
			_spans[k]=t;
			_spans[k].left=&_spans[index1];
			_spans[k].right=&_spans[index2+2*nbest+k2];
		}
		//0+4->0
		for(k2=0;k2<nbest && _spans[index2+4*nbest+k2].y!=-2;k2++){
			t.score=_spans[index1].score+_spans[index2+4*nbest+k2].score;
			ip = equal_range(_spans.begin()+index3, _spans.begin()+index3+nbest, t, span_cmp());
			pos=ip.second-(_spans.begin()+index3);
			if(pos>=nbest) break;
			//insert
			for(k=index3+nbest-1;k>index3+pos;k--)
				_spans[k]=_spans[k-1];
			_spans[k]=t;
			_spans[k].left=&_spans[index1];
			_spans[k].right=&_spans[index2+4*nbest+k2];
		}
	}
	//3+0->0
	for(k1=0;k1<nbest && _spans[index1+3*nbest+k1].y!=-2;k1++){
		for(k2=0;k2+k1+k1*k2<nbest && _spans[index2+k2].y!=-2;k2++){
			t.score=_spans[index1+3*nbest+k1].score+_spans[index2+k2].score;
			ip = equal_range(_spans.begin()+index3, _spans.begin()+index3+nbest, t, span_cmp());
			pos=ip.second-(_spans.begin()+index3);
			if(pos>=nbest) break;
			//insert
			for(k=index3+nbest-1;k>index3+pos;k--)
				_spans[k]=_spans[k-1];
			_spans[k]=t;
			_spans[k].left=&_spans[index1+3*nbest+k1];
			_spans[k].right=&_spans[index2+k2];
		}
	}
	
	if(pos3<l){
	//3+1->1
		for(k1=0;k1<nbest && _spans[index1+3*nbest+k1].y!=-2;k1++){
			for(k2=0;k2+k1+k1*k2<nbest && _spans[index2+nbest+k2].y!=-2;k2++){
				t.score=_spans[index1+3*nbest+k1].score+_spans[index2+nbest+k2].score;
				ip = equal_range(_spans.begin()+index3+nbest, _spans.begin()+index3+2*nbest, t, span_cmp());
				pos=ip.second-(_spans.begin()+index3+nbest);
				if(pos>=nbest) break;
				//insert
				for(k=index3+2*nbest-1;k>index3+nbest+pos;k--)
					_spans[k]=_spans[k-1];
				_spans[k]=t;
				_spans[k].left=&_spans[index1+3*nbest+k1];
				_spans[k].right=&_spans[index2+nbest+k2];
			}
		}
	}

	if(pos3<l){
	//3+3->1
		for(k1=0;k1<nbest && _spans[index1+3*nbest+k1].y!=-2;k1++){
			for(k2=0;k2+k1+k1*k2<nbest && _spans[index2+3*nbest+k2].y!=-2;k2++){
				t.score=_spans[index1+3*nbest+k1].score+_spans[index2+3*nbest+k2].score;
				ip = equal_range(_spans.begin()+index3+nbest, _spans.begin()+index3+2*nbest, t, span_cmp());
				pos=ip.second-(_spans.begin()+index3+nbest);
				if(pos>=nbest) break;
				//insert
				for(k=index3+2*nbest-1;k>index3+nbest+pos;k--)
					_spans[k]=_spans[k-1];
				_spans[k]=t;
				_spans[k].left=&_spans[index1+3*nbest+k1];
				_spans[k].right=&_spans[index2+3*nbest+k2];
			}
		}
	}
	//4+2->2
	for(k1=0;k1<nbest && _spans[index1+4*nbest+k1].y!=-2;k1++){
		for(k2=0;k2+k1+k1*k2<nbest && _spans[index2+2*nbest+k2].y!=-2;k2++){
			t.score=_spans[index1+4*nbest+k1].score+_spans[index2+2*nbest+k2].score;
			ip = equal_range(_spans.begin()+index3+2*nbest, _spans.begin()+index3+3*nbest, t, span_cmp());
			pos=ip.second-(_spans.begin()+index3+2*nbest);
			if(pos>=nbest) break;
			//insert
			for(k=index3+3*nbest-1;k>index3+2*nbest+pos;k--)
				_spans[k]=_spans[k-1];
			_spans[k]=t;
			_spans[k].left=&_spans[index1+4*nbest+k1];
			_spans[k].right=&_spans[index2+2*nbest+k2];
		}
	}
	//4+4->2
	for(k1=0;k1<nbest && _spans[index1+4*nbest+k1].y!=-2;k1++){
		for(k2=0;k2+k1+k1*k2<nbest && _spans[index2+4*nbest+k2].y!=-2;k2++){
			t.score=_spans[index1+4*nbest+k1].score+_spans[index2+4*nbest+k2].score;
			ip = equal_range(_spans.begin()+index3+2*nbest, _spans.begin()+index3+3*nbest, t, span_cmp());
			pos=ip.second-(_spans.begin()+index3+2*nbest);
			if(pos>=nbest) break;
			//insert
			for(k=index3+3*nbest-1;k>index3+2*nbest+pos;k--)
				_spans[k]=_spans[k-1];
			_spans[k]=t;
			_spans[k].left=&_spans[index1+4*nbest+k1];
			_spans[k].right=&_spans[index2+4*nbest+k2];
		}
	}
	return false;
}

void crfparser_thread::split_span(int index1,int index2,int index3, int index4, int index5){
	if(_has_beta[index5] && _has_alpha[index4]){
		if(_has_beta[index1]){
			_beta[index1]=log_sum_exp(_beta[index1],_alpha[index4]+_beta[index5]);
		}else{
			_has_beta[index1]=true;
			_beta[index1]=_alpha[index4]+_beta[index5];
		}
	}
	if(_has_beta[index5] && _has_alpha[index3]){
		if(_has_beta[index2]){
			_beta[index2]=log_sum_exp(_beta[index2],_alpha[index3]+_beta[index5]);
		}else{
			_has_beta[index2]=true;
			_beta[index2]=_alpha[index3]+_beta[index5];
		}
	}
}

void crfparser_thread::close_span(int i,int j,int l){
	//i<j
	int k;
	int ysize=base_feature::_labels.size();
	int index=cell_index(i,j,l)*5;
	if(j!=l){//CLOSE LEFT i->j, 0->3
		_has_alpha[index+3]=true;
		_alpha[index+3]=_alpha[index]+_lattice[i*l+j];
		_has_alpha[index+4]=true;
		_alpha[index+4]=_alpha[index]+_lattice[j*l+i];
	}else{//CLOSE RIGHT i<-l
		_alpha[index+4]=_alpha[index]+_lattice[i*l+i];
		_has_alpha[index+4]=true;
	}
}

void crfparser_thread::close_span_nbest(int i,int j,int l){
	//i<j
	int k1,k2,pos,k;
	int nbest=_c->_nbest;
	span t;
	t.head=i;
	t.tail=j;
	pair<vector< span >::const_iterator,vector< span >::const_iterator> ip;
	int ysize=base_feature::_labels.size();
	int index=cell_index(i,j,l)*5*nbest;
	if(j!=l){//CLOSE LEFT i->j
		for(k1=0;k1<nbest && _spans[index+k1].y!=-2;k1++){
			//0->3
			span &s=_spans[index+k1];
			for(k2=0;k1*k2+k1+k2<nbest && _best_label[nbest*(i*l+j)+k2]!=-2;k2++){
				t.score=s.score+_lattice[nbest*(i*l+j)+k2];
				t.y=_best_label[nbest*(i*l+j)+k2];
				ip=equal_range(_spans.begin()+index+3*nbest, _spans.begin()+index+4*nbest, t, span_cmp());
				pos=ip.second-(_spans.begin()+index+3*nbest);
				if(pos>=nbest) break;
				//insert
				for(k=index+4*nbest-1;k>index+3*nbest+pos;k--)
					_spans[k]=_spans[k-1];
				_spans[k]=t;
				_spans[k].left=&_spans[index+k1];
				_spans[k].right=NULL;
			}
			//0->4	CLOSE RIGHT i->j
			for(k2=0;k1*k2+k1+k2<nbest && _best_label[nbest*(j*l+i)+k2]!=-2;k2++){
				t.score=s.score+_lattice[nbest*(j*l+i)+k2];
				t.y=_best_label[nbest*(j*l+i)+k2];
				ip=equal_range(_spans.begin()+index+4*nbest, _spans.begin()+index+5*nbest, t, span_cmp());
				pos=ip.second-(_spans.begin()+index+4*nbest);
				if(pos>=nbest) break;
				//insert
				for(k=index+5*nbest-1;k>index+4*nbest+pos;k--)
					_spans[k]=_spans[k-1];
				_spans[k]=t;
				_spans[k].left=NULL;
				_spans[k].right=&_spans[index+k1];
			}
		}
	}else{//CLOSE RIGHT i<-l
		//0->4
		for(k1=0;k1<nbest && _spans[index+k1].y!=-2;k1++){
			span &s=_spans[index+k1];
			t.score=s.score+_lattice[nbest*(i*l+i)];
			t.y=-1;
			ip=equal_range(_spans.begin()+index+4*nbest, _spans.begin()+index+5*nbest, t, span_cmp());
			pos=ip.second-(_spans.begin()+index+4*nbest);
			if(pos>=nbest) break;
			//insert
			for(k=index+5*nbest-1;k>index+4*nbest+pos;k--)
				_spans[k]=_spans[k-1];
			_spans[k]=t;
			_spans[k].left=NULL;
			_spans[k].right=&_spans[index+k1];
		}
	}
}


void crfparser_thread::disclose_span(int i,int j,int l){
	//i<j
/*
chart cell order:
						bL bR s
0						F  F  (T, if |i-j|=1; F else)
1						F  T  F
2						T  F  F
3						F  T  T
4						T  F  T
*/

	int ysize=base_feature::_labels.size();
	int index=cell_index(i,j,l)*5;
	if(j!=l){//DISCLOSE LEFT
		if(_has_beta[index]){
			_beta[index]=log_sum_exp(_beta[index],_beta[index+3]+_lattice[i*l+j]);
		}else{
			_has_beta[index]=true;
			_beta[index]=_beta[index+3]+_lattice[i*l+j];
		}
		//DISCLOSE RIGHT
		if(_has_beta[index]){
			_beta[index]=log_sum_exp(_beta[index],_beta[index+4]+_lattice[j*l+i]);
		}else{
			_has_beta[index]=true;
			_beta[index]=_beta[index+4]+_lattice[j*l+i];
		}
	}else{//DISCLOSE RIGHT
		if(_has_beta[index]){
			_beta[index]=log_sum_exp(_beta[index],_beta[index+4]+_lattice[i*l+i]);
		}else{
			_beta[index]=_beta[index+4]+_lattice[i*l+i];
			_has_beta[index]=true;
		}
	}
}
void crfparser_thread::build_lattice(sentence &sen){
//_lattice[i*l+j]=edge(i->j)  i!=j
//_lattice[i*l+j]=edge(i<-root), i==j
//_edgep[(i*l+j)*ysize+y]=edge(i->j,y), i!=j
	int l=sen.l;
	int i,j,k,ii;
	int ysize=base_feature::_labels.size();
	if(_lattice.size()<l*l)
		_lattice.resize(l*l);
	if(_edgep.size()<l*l*ysize)
		_edgep.resize(l*l*ysize);//relation type probability on each edge
	fill(_lattice.begin(),_lattice.end(),0);
	fill(_edgep.begin(),_edgep.end(),0);
	for(i=0;i<l;i++){
		for(j=i+1;j<l;j++){
			//i->j
			for(k=0;k<ysize;k++)
				for(ii=0;ii<sen.fnum[i*l+j][k];ii++)
					_edgep[ysize*(i*l+j)+k]+=_c->_lambda[sen.fvector[i*l+j][k][ii]];
			//i<-j
			for(k=0;k<ysize;k++)
				for(ii=0;ii<sen.fnum[j*l+i][k];ii++)
					_edgep[ysize*(j*l+i)+k]+=_c->_lambda[sen.fvector[j*l+i][k][ii]];
		}
		//i<-root
		for(j=0;j<sen.fnum[i*l+i][0];j++){
			_lattice[i*l+i]+=_c->_lambda[sen.fvector[i*l+i][0][j]];
		}

		for(k=0;k<ysize;k++){
			//i act as parent, i->j
			for(ii=0;ii<sen.fnum[l*l+i][k];ii++){
				for(j=0;j<i;j++)
					_edgep[ysize*(i*l+j)+k]+=_c->_lambda[sen.fvector[l*l+i][k][ii]];
				for(j=i+1;j<l;j++)
					_edgep[ysize*(i*l+j)+k]+=_c->_lambda[sen.fvector[l*l+i][k][ii]];

			}
			//i act as child, i<-j
			for(ii=0;ii<sen.fnum[l*l+l+i][k];ii++){
				for(j=0;j<i;j++)
					_edgep[ysize*(j*l+i)+k]+=_c->_lambda[sen.fvector[l*l+l+i][k][ii]];
				for(j=i+1;j<l;j++)
					_edgep[ysize*(j*l+i)+k]+=_c->_lambda[sen.fvector[l*l+l+i][k][ii]];

			}
		}
	}
}

void crfparser_thread::collapse_lattice(sentence &sen){
	int l=sen.l;
	int i,j,k,ii;
	int ysize=base_feature::_labels.size();
	for(i=0;i<l;i++){
		for(j=i+1;j<l;j++){//collapse feature _edgep->_lattice
			//i->j
			_lattice[i*l+j]=_edgep[ysize*(i*l+j)];
			for(k=1;k<ysize;k++){
				_lattice[i*l+j]=log_sum_exp(_lattice[i*l+j],_edgep[ysize*(i*l+j)+k]);
			}
			//i<-j
			_lattice[j*l+i]=_edgep[ysize*(j*l+i)];
			for(k=1;k<ysize;k++){
				_lattice[j*l+i]=log_sum_exp(_lattice[j*l+i],_edgep[ysize*(j*l+i)+k]);
			}
		}
	}
}


void crfparser_thread::build_lattice_max(sentence &sen){
//_lattice[i*l+j]=edge(i->j)  i!=j
//_lattice[i*l+j]=edge(i<-root), i==j
//_edgep[ysize*(i*l+j)+y]=edge(i->j,y), i!=j
	int l=sen.l;
	int i,j,k,ii;
	int pos;
	int ysize=base_feature::_labels.size();
	int nbest=_c->_nbest;
	pair<vector< double >::const_iterator,vector< double >::const_iterator> ip;
	
	if(_lattice.size()<l*l*nbest){
		_lattice.resize(l*l*nbest);
		_best_label.resize(l*l*nbest);
	}
	fill(_lattice.begin(),_lattice.end(),-INF);
	fill(_best_label.begin(),_best_label.end(),-2);
	for(i=0;i<l;i++){
		for(j=0;j<l;j++){
			if(i==j) continue;
			//i->j
			for(k=0;k<ysize;k++){
				double logp=0;
				for(ii=0;ii<sen.fnum[i*l+j][k];ii++)
					logp+=_c->_lambda[sen.fvector[i*l+j][k][ii]];
				for(ii=0;ii<sen.fnum[l*l+i][k];ii++)
					logp+=_c->_lambda[sen.fvector[l*l+i][k][ii]];
				for(ii=0;ii<sen.fnum[l*l+l+j][k];ii++)
					logp+=_c->_lambda[sen.fvector[l*l+l+j][k][ii]];
				ip = equal_range(_lattice.begin()+(i*l+j)*nbest, _lattice.begin()+(i*l+j+1)*nbest, logp, inverse_cmp<double>());
				pos=ip.second-(_lattice.begin()+(i*l+j)*nbest);
				if(pos>=nbest) continue;
				//insert
				for(ii=(i*l+j+1)*nbest-1;ii>(i*l+j)*nbest+pos;ii--){
					_best_label[ii]=_best_label[ii-1];
					_lattice[ii]=_lattice[ii-1];
				}
				_lattice[ii]=logp;
				_best_label[ii]=k;
			}
		}
		//i<-root
		_lattice[(i*l+i)*nbest]=0;
		_best_label[(i*l+i)*nbest]=-1;
		for(j=0;j<sen.fnum[i*l+i][0];j++)
			_lattice[(i*l+i)*nbest]+=_c->_lambda[sen.fvector[i*l+i][0][j]];
	}
}



int crfparser_thread::cell_index(int i, int j,int l){//return chart cell(i,j) position in chart,  i<j
	if(i>j){
		int k=i;
		i=j;
		j=k;
	}
	if(j<l)
		return (l+l-(i-1))*i/2+j-i;
	else
		return (l+l-(i-1))*i/2;// on the diagonal 
}


void crfparser_thread::calculate_margin(sentence &sen, double &z){
//	_margin[ysize*(i*l+j)+k], 
//		i!=j
//			marginal probability of i->j with label k
//		i==j
//			root->i
//	_par_margin[k*l+i]+=\sum_j _margin[ysize*(i*l+j)+k];, marginal probability of node i being parent.
//	_child_margin[k*l+j]+=\sum_i _margin[ysize*(i*l+j)+k];, marginal probability of node j being child.
	int i,j,k;
	int l=sen.l;
	int ysize=base_feature::_labels.size();
	if(_margin.size()<ysize*l*l){
		_margin.resize(ysize*l*l);
		_par_margin.resize(ysize*l);
		_child_margin.resize(ysize*l);
	}
	fill(_margin.begin(),_margin.end(),0);
	fill(_par_margin.begin(),_par_margin.end(),0);
	fill(_child_margin.begin(),_child_margin.end(),0);
	for(i=0;i<l;i++){
		for(j=i+1;j<l;j++){
			int index=cell_index(i,j,l);
			double p3=_alpha[index*5+3]+_beta[index*5+3]-z;//log(p(i->j))
			double p4=_alpha[index*5+4]+_beta[index*5+4]-z;//log(p(i->j))
			for(k=0;k<ysize;k++){
				//i->j
				_margin[ysize*(i*l+j)+k]=exp(_edgep[ysize*(i*l+j)+k]-_lattice[i*l+j]+p3);
				_par_margin[k*l+i]+=_margin[ysize*(i*l+j)+k];
				_child_margin[k*l+j]+=_margin[ysize*(i*l+j)+k];
				//i<-j
				_margin[ysize*(j*l+i)+k]=exp(_edgep[ysize*(j*l+i)+k]-_lattice[j*l+i]+p4);
				_par_margin[k*l+j]+=_margin[ysize*(j*l+i)+k];
				_child_margin[k*l+i]+=_margin[ysize*(j*l+i)+k];
			}
		}
		//i<-root
		int index=cell_index(i,l,l);
		_margin[ysize*(i*l+i)]=exp(_alpha[index*5+4]+_beta[index*5+4]-z);
	}
}
void crfparser_thread::calculate_gradient(sentence &sen, double &z){
/*
chart cell order:
						bL bR s
0						F  F  (T, if |i-j|=1; F else)
1						F  T  F
2						T  F  F
3						F  T  T
4						T  F  T
*/
	int i,j,k,ii,y;
	int l=sen.l;
	int ysize=base_feature::_labels.size();
	for(i=0;i<l;i++){
		for(j=i+1;j<l;j++){
			for(y=0;y<ysize;y++){
				//i->j
				for(k=0;k<sen.fnum[i*l+j][y];k++)
					_gradient[sen.fvector[i*l+j][y][k]]+=_margin[ysize*(i*l+j)+y];
				//i<-j
				for(k=0;k<sen.fnum[j*l+i][y];k++)
					_gradient[sen.fvector[j*l+i][y][k]]+=_margin[ysize*(j*l+i)+y];
			}
		}
		//i<-root
		for(j=0;j<sen.fnum[i*l+i][0];j++)
			_gradient[sen.fvector[i*l+i][0][j]]+=_margin[ysize*(i*l+i)];
//i*l+j i!=j, i->j
//      i==j, i<-root
//l*l+i i act as parent
//l*l+l+i act as child
		for(y=0;y<ysize;y++){
			//i act as parent
			for(j=0;j<sen.fnum[l*l+i][y];j++)
				_gradient[sen.fvector[l*l+i][y][j]]+=_par_margin[y*l+i];
			//i act as child
			for(j=0;j<sen.fnum[l*l+l+i][y];j++)
				_gradient[sen.fvector[l*l+l+i][y][j]]+=_child_margin[y*l+i];
		}
	}
	
	for(j=0;j<l;j++){
		i=sen.parent[j];//i->j
		if(i==l){
			for(k=0;k<sen.fnum[j*l+j][0];k++)
				_gradient[sen.fvector[j*l+j][0][k]]--;
		}else{
			k=sen.label[j];
			//i->j
			for(ii=0;ii<sen.fnum[i*l+j][k];ii++)
				_gradient[sen.fvector[i*l+j][k][ii]]--;
			//i act as parent
			for(ii=0;ii<sen.fnum[l*l+i][k];ii++)
				_gradient[sen.fvector[l*l+i][k][ii]]--;
			//j act as child
			for(ii=0;ii<sen.fnum[l*l+l+j][k];ii++)
				_gradient[sen.fvector[l*l+l+j][k][ii]]--;
		}
	}
}

void crfparser_thread::calculate_gradient_factor(sentence &sen, double &z){
//_lattice[i*l+j]=edge(i->j)  i!=j
//_lattice[i*l+j]=edge(i<-root), i==j
//_edgep[(i*l+j)*ysize+y]=edge(i->j,y), i!=j
	int i,j,k,ii,y;
	int l=sen.l;
	int ysize=base_feature::_labels.size();
	for(i=0;i<l;i++){
		for(j=i+1;j<l;j++){
			for(y=0;y<ysize;y++){
				//i->j
				_c->_factor_gradient+=_margin[ysize*(i*l+j)+y]*_edgep[(i*l+j)*ysize+y];
				//i<-j
				_c->_factor_gradient+=_margin[ysize*(j*l+i)+y]*_edgep[(j*l+i)*ysize+y];
			}
		}
		//i<-root
		_c->_factor_gradient+=_margin[ysize*(i*l+i)]*_lattice[i*l+i];
	}
}

double crfparser_thread::tree_cost(sentence &sen){
	double s=0;
	int ysize=base_feature::_labels.size();
	int i,j,k;
	int l=sen.l;
	for(i=0;i<l;i++){
		j=sen.parent[i];
		if(j!=l){
			k=sen.label[i];
			s+=_edgep[ysize*(j*l+i)+k];
		}else{
			s+=_lattice[i*l+i];
		}
	}
	return s;
}

void crfparser_thread::cyk(sentence &sen){
/*
chart cell order:
						bL bR s
0						F  F  (T, if |i-j|=1; F else)
1						F  T  F
2						T  F  F
3						F  T  T (max)
4						T  F  T (max)

chart[i*l+j]=edge(i,j), i<j
             edge(i,root), i==j

JOINT(i,k) (k,j)

if(k==i+1)
	0+2->0
	0+4->0
3+0->0
3+1->1
3+3->1
4+2->2
4+4->2

JOINT(i,j) (j,n)
4+4->2
4+2->2
3+0->0
*/
	int ysize=base_feature::_labels.size();
	int l=sen.l;
	int i,j,k;
	int index;
	int n;
	int cell_num=l*(l+1)/2;
	int nbest=_c->_nbest;
	if(_spans.size()<cell_num*5*nbest){
		_spans.resize(cell_num*5*nbest);
		_cur_par.resize(l);
		_cur_label.resize(l);
	}

	for(i=0;i<cell_num*5*nbest;i++)
		_spans[i].y=-2;

	for(i=0;i<l;i++){
		index=cell_index(i,i+1,l)*5*nbest;
		//ADD(SEED(i))
		_spans[index].score=0;
		_spans[index].y=-1;
		_spans[index].left=NULL;
		_spans[index].right=NULL;
		_spans[index].head=i;
		_spans[index].tail=i+1;
		close_span_nbest(i,i+1,l);
	}
	for(n=2;n<=l;n++){//n=span len
		for(i=0;i<=l-n;i++){
			j=i+n;
			for(k=i+1;k<j;k++)
				merge_span_nbest(i,k,j,l);
			close_span_nbest(i,j,l);
		}
	}
	
	index=cell_index(0,l,l)*5*nbest;
	_par.clear();
	_label.clear();

	for(i=0,j=0;_par.size()<nbest;){
		fill(_cur_par.begin(),_cur_par.end(),-1);
		fill(_cur_label.begin(),_cur_label.end(),-2);
		if(_spans[index+4*nbest+j].y==-2){
			if(_spans[index+2*nbest+i].y==-2)
				break;
			cyk_trace(_spans[index+2*nbest+i]);
			i++;
		}else if(_spans[index+2*nbest+i].y==-2){
			_cur_label[0]=-1;
			_cur_par[0]=l;
			cyk_trace(_spans[index+4*nbest+j]);
			j++;
		}else if(_spans[index+2*nbest+i].score>=_spans[index+4*nbest+j].score){
			cyk_trace(_spans[index+2*nbest+i]);
			i++;
		}else{
			_cur_label[0]=-1;
			_cur_par[0]=l;
			cyk_trace(_spans[index+4*nbest+j]);
			j++;
		}
		_par.push_back(_cur_par);
		_label.push_back(_cur_label);
	}
}

void crfparser_thread::cyk_trace(span &s){
	if(s.y==-2)
		return;
	if(s.right==NULL && s.left==NULL)
		return; //SEED(i)
	if(s.right==NULL){
		_cur_par[s.tail]=s.head;
		_cur_label[s.tail]=s.y;
		cyk_trace(*s.left);
	}else if(s.left==NULL){
		_cur_par[s.head]=s.tail;
		_cur_label[s.head]=s.y;
		cyk_trace(*s.right);
	}else{
		cyk_trace(*s.left);
		cyk_trace(*s.right);
	}
}

void crfparser_thread::assign_label(sentence &s,int *par, int *label){
	int i;
	for(i=0;i<s.l;i++){
		s.label[i]=label[i];
		s.parent[i]=par[i];

	}
}
