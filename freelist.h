/*
  freelist.h 2007-03-29

  Modified from CRF++ toolkit http://sourceforge.net/projects/crfpp/
  
  Author taku
*/

#ifndef FREELIST_H
#define FREELIST_H
//block array
#include <vector>
using namespace std;

template <class T>
class length//for any T[], its default length is 1
{
	public:
	size_t operator()(const T *str) const { return 1; }
};

template <class T, class length_func = length<T> >
class freelist {
private:
	vector <T *> element;
	size_t c;//col
	size_t r;//row
	size_t size;
	bool empty;
public:
	void free(){r=c=0;}
	T* alloc(size_t len = 1)
	{
		if ((c+len) >= size)
		{	
			r++;
			c=0;
		}
		if (r == element.size()) 
			element.push_back(new T[size]);
		T* ret=element[r]+c;
		c+=len;
		if(len)
			empty=false;
		return ret;
	}

	T* push_back(T *src, size_t len = 0)
	{
		if(!len)
			len=length_func()(src);
		T *ret=alloc(len);
		if(src==0)
			memset(ret,0,len*sizeof(T));
		else
			memcpy(ret,src,len*sizeof(T));
		return ret;
	}

	void set_size(size_t n) { size = n; }
	explicit freelist(size_t _size): c(0), r(0), size(_size) ,empty(true){}
	explicit freelist(): c(0), r(0), size(0) {}

	void clear()
	{
		if(empty) return;
		for (r = 0; r < element.size(); r++)
			delete [] element[r];
		r=c=0;
		element.clear();
		vector<T *>(element).swap(element);
		empty=true;
	}
	virtual ~freelist()
	{
		clear();
	}
};


#endif
