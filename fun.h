#ifndef FUN_H
#define FUN_H
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <string.h>
#include <set>
using namespace std;
const int MINUS_LOG_EPSILON=50;

char* StrStr(char *s, char *c);
bool split_string(char *str, const char *cut, vector<char *> &strs);
char* catch_string(char *str, char *head, char* tail,char* catched);//catch the first substring catched between head and tail in str
char* catch_string(char *str, char* tail, char* catched);//catch the first substring before tail
void trim_line(char *line);//get rid of '\r\n'
int vector_cmp(const vector<int> & s, const vector<int> & t);

template <class T> inline T _min(T x, T y) { return(x < y) ? x : y; }
template <class T> inline T _max(T x, T y) { return(x > y) ? x : y; }

inline double log_sum_exp(double x,double y)
{
    double vmin = _min(x, y);
    double vmax = _max(x, y);
    if (vmax > vmin + MINUS_LOG_EPSILON) {
      return vmax;
    } else {
      return vmax + log(exp(vmin - vmax) + 1.0);
    }
}




template <class T, class cmp_func>
bool vector_search(vector<T> &v, const T & s, int &index, int &insert_pos, cmp_func cmp)
{
	pair<typename vector< T >::const_iterator, typename vector< T >::const_iterator> ip;
	ip = equal_range( v.begin( ), v.end( ), s , cmp);
	index=ip.first-v.begin();
	insert_pos=ip.second-v.begin();
	if ( ip.first == ip.second )//not found
		return false;
	return true;
}



template <class T>
bool vector_search(vector<T> &v, const T & s, int &index, int &insert_pos)
{
	pair<typename vector< T >::const_iterator,typename vector< T >::const_iterator> ip;
	ip = equal_range( v.begin( ), v.end( ), s );
	index=ip.first-v.begin();
	insert_pos=ip.second-v.begin();
	if ( ip.first == ip.second )//not found
	{
		return false;
	}
	return true;
}

template <class T>
bool vector_insert(vector<T> &v, const T & s, int index)
{
	if(index>v.size())return false;
	typename vector<T>::iterator t=v.begin()+index;
	v.insert(t,s);
	return true;
}
template <class T>
void combine(vector<T> &v, int left, int m, int right, vector<int> &index)
{
	vector<T> tempv(v.begin()+left,v.begin()+right+1);
	vector<int> tempindex(index.begin()+left,index.begin()+right+1);

	int left_size=m-left+1;
	int size=right-left+1;
	int middle=m-left+1;
	int i=0;
	int j=middle;
	int k=left;
	while(i<left_size && j<size)
	{
		if(tempv[i]>=tempv[j])
		{
			v[k]=tempv[i];
			index[k]=tempindex[i];
			k++;
			i++;
		}else{
			v[k]=tempv[j];
			index[k]=tempindex[j];
			k++;
			j++;
		}
	}
	while(i<left_size)
	{
		v[k]=tempv[i];
		index[k]=tempindex[i];
		k++;
		i++;
	}
}


template <class T>
void merge_sort(vector<T> &v, int left, int right, vector<int> &index)
{
    if (left<right)
    {
        int m=(left+right)/2;
        merge_sort(v,left, m,index);
        merge_sort(v,m+1, right, index);
        combine(v,left, m, right,index);
    }
}


template <class T>
void merge_sort(vector<T> v, vector<int> &index)
//index[i] is the original index of i th best element
{
	index.clear();
	index.resize(v.size());
	for(int i=0;i<v.size();i++)
		index[i]=i;
	merge_sort(v,0,v.size()-1,index);
}

class str_cmp
{
public:
	bool operator()(const char* s, const char* t) const 
	{
		return strcmp(s,t)<0;
	}
};


template <class T>
int bisearch(T *v, const T & s, int low, int high)
{
	if(high<low)
		return -1;
	int mid=(low+high)/2;
	if(v[mid]==s)
		return mid;
	else if(v[mid]>s)
		return bisearch(v,s,low,mid-1);
	else
		return bisearch(v,s,mid+1,high);
}


template <class T>
int bisearch(T *v, const T & s, int length)//return -1 if not found
{
	return bisearch(v,s,0,length-1);
}

template <class T>
class inverse_cmp
{
	public:
	bool operator()(T s, T t) const 
	{
		return s>t;
	}
};


bool build_table(ifstream &fin, vector<vector<string> > & table, int &lines);

int get_feature_grid(int len, int p, int c, bool has_p, bool has_c);
int get_lattice_grid(int i, int j, int l);

class str_length{
	public:
		size_t operator()(const char *str) const  {return strlen(str)+1;}
};


template<class T>
bool is_subset(set<T> &s1, set<T> &s2){//if s1 is empty, return true; else if s2 is empty return false;
	typename set<T>::iterator iter;
	for (iter=s1.begin();iter!=s1.end();++iter){
		if (s2.find(*iter)==s2.end())
			return false;
	}
	return true;
}

template<class T>
bool is_inter(set<T>& s1, set<T>& s2){//if s1 or s2 is empty, return false;
	typename set<T>::iterator iter;
	for (iter=s1.begin();iter!=s1.end();++iter){
		if (s2.find(*iter)!=s2.end())
			return true;
	}
	return false;
}

template<class T>
set<T> set_sub(set<T>& s1, set<T>& s2){//s1-s2
	set<T> r;
	typename set<T>::iterator iter;
	for (iter=s1.begin();iter!=s1.end();++iter){
		if (s2.find(*iter)==s2.end())
			r.insert(*iter);
	}
	return r;
}

template<class T>
set<T> set_union(set<T>& s1, set<T>& s2){//s1-s2
	set<T> r;
	typename set<T>::iterator iter;
	for (iter=s1.begin();iter!=s1.end();++iter)
			r.insert(*iter);
	for (iter=s2.begin();iter!=s2.end();++iter)
			r.insert(*iter);
	return r;
}


#endif