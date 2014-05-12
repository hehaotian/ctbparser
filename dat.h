#ifndef DAT_H
#define DAT_H

#include <vector>
#include <cstdio>
#include <fstream>
#include <string.h>
#include <stdlib.h>
class double_array_trie{
private:
	struct data_t{
		int base;
		unsigned int check;
	};
	
	struct node_t{
		unsigned int begin,end,depth,pos;
		node_t(unsigned int begin, unsigned int end, unsigned int depth, unsigned int pos){
			this->begin = begin;
			this->end = end;
			this->depth = depth;
			this->pos = pos;
		}
	};
	
	unsigned int _length,_next_pos,_size;
	
	data_t* _data;
	
	void realloc_data(unsigned int size){
		unsigned int pre = _length;
		data_t* temp = new data_t[size];
		data_t zero;
		zero.base = 0;
		zero.check = 0;
		if (!temp) exit(-1);	
		for (unsigned int i = 0; i < pre; ++i)
			temp[i] = _data[i];
		for (unsigned int i = pre; i < size; ++i)
			temp[i] = zero;
		delete []_data;
		_data = temp;
		_length = size;
	}
	
	unsigned int get_pos(unsigned char* list, unsigned int num){
		unsigned int i,pos,ans = 0, nonzeronum=0;
		pos = _next_pos;
		bool first = false;
		if (pos<list[0])
			pos = list[0];
		while (pos >= _length)
			realloc_data((unsigned int)(_length*1.1));
		while (true){
			pos++;
			while (pos >= _length)
				realloc_data((unsigned int)(_length*1.1));
			if (_data[pos].check!=0){
				if (first) nonzeronum++;
				continue;
			}else if (!first){
				_next_pos = pos;
				first = true;
			}
			while (_length <= pos + list[num-1] - list[0])
				realloc_data((unsigned int)(_length*1.1));
			for (i = 1; i<num; ++i)
				if (_data[pos+list[i]-list[0]].check!=0) 
				break;
			if (i==num){
				ans = pos - list[0];
				break;
			}
		}
		if (1.0*nonzeronum/(pos - _next_pos + 1) > 0.95)
			_next_pos = pos;
		return ans;
	}

    void clear(){
		if(_data)
			delete []_data;
		_data=NULL;
	}
	

	void insert(node_t now , unsigned char** keys, int* ids){
		std::vector<node_t> sonlist;
		unsigned char** iter;
		unsigned int left,right;
		unsigned char list[256];
		unsigned int i;
		iter = keys + now.begin;
		if (strlen((char*)(*iter))<now.depth){
			_data[now.pos].base = ids[now.begin];
			return;
		}
		list[0] = (*iter)[now.depth];
		i = 0;
		iter++;
		unsigned char** end_pos = keys + now.end;
		for (;iter!=end_pos;++iter)
			if ((*iter)[now.depth]!=list[i])
				list[++i] = (*iter)[now.depth];	
		unsigned int new_base = get_pos(list,i+1);
		_data[now.pos].base = new_base;
		left = now.begin;
		i = 0;
		for (unsigned int k = now.begin;k<now.end;++k){
			if (keys[k][now.depth]!=list[i]){
				unsigned int pos = new_base + list[i];
				_data[pos].check = now.pos;
				i++;
				right = k;
				sonlist.push_back(node_t(left,right,now.depth+1,pos));
				left = k;
			}
		}
		unsigned int pos = new_base + list[i];
		_data[pos].check = now.pos;
		right = now.end;
		sonlist.push_back(node_t(left,right,now.depth+1,pos));
		std::vector<node_t>::iterator soniter;
		for (soniter = sonlist.begin();soniter!=sonlist.end();++soniter) 
			insert(*soniter, keys, ids);
	}
public:
	
	double_array_trie(){
		_data = NULL;
	}
	
	~double_array_trie(){
		clear();
	}

	unsigned int size(){
		return _size;
	}

	int build(unsigned int size, unsigned char** keys, int* ids){
		_size=size;
		_length = 1024*1024;
		_data = new data_t[_length];
		if (!_data) exit(-1);
		for (unsigned int i=0;i<_length;++i){
			_data[i].base = 0;
			_data[i].check = 0;
		}
		_next_pos = 1;
		_data[1].base = 1;
		if(size>0)
			insert(node_t(0,size,0,1), keys, ids);
		realloc_data(_length+256);
		return 0;
	}
	
	inline bool search(const unsigned char* key, int& result){
		register unsigned int s = 1,t;
		register const unsigned char* k = key; 
		while (*k)
		{
			t = _data[s].base + *k;
			if (_data[t].check!=s) return false;
			s = t;
			++k;
		}
		t = _data[s].base;
		if (_data[t].check!=s) return false;
		result = _data[t].base;
		return true;
	}
	
	inline unsigned int prefix(const unsigned char* key, int*& result){
		register unsigned int num = 0;
		register unsigned int s = 1,t;
		register const unsigned char* k = key; 
		while (*k)
		{
			t = _data[s].base + *k;
			if (_data[t].check!=s) break;
			s = t;
			t = _data[s].base;//+'\0'
			if (_data[t].check==s) result[num++] = _data[t].base;
			++k;
		}
		return num;
	}
	
	bool load(FILE*& fp){
		if (!fread(&_length, sizeof(unsigned int), 1, fp)) return false;
		if (!fread(&_size, sizeof(unsigned int), 1, fp)) return false;
		_data = new data_t[_length];
		if (_length != fread(_data,sizeof(data_t),_length,fp)) return false;
		return true;
	}
	
	bool load(std::ifstream& in){
		in.read((char*)&_length,sizeof(unsigned int));
		in.read((char*)&_size,sizeof(unsigned int));
		_data = new data_t[_length];
		in.read((char*)_data,sizeof(data_t)*_length);
		return true;
	}
	
	bool write(FILE*& fp){
		if (!fwrite(&_length, sizeof(unsigned int), 1, fp)) return false;
		if (!fwrite(&_size, sizeof(unsigned int), 1, fp)) return false;
		if (_length!=fwrite(_data, sizeof(data_t), _length, fp)) return false;
		return true;
	}
	
	bool write(std::ofstream& out){
		out.write((char*)&_length,sizeof(unsigned int));
		out.write((char*)&_size,sizeof(unsigned int));
		out.write((char*)_data,sizeof(data_t)*_length);
		return true;
	}
	
};
#endif
