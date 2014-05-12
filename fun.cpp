#include "const.h"
#include "fun.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


char* StrStr(char *s, char *c)
{
	int len=strlen(s);
	for(int i=0;i<len;i+=2)
		if(!strncmp(s+i,c,2))
			return s+i;
	return NULL;
}



void trim_line(char *line)//get rid of '\r\n'
{
	size_t len=strlen(line);
	if(len && line[len-1]=='\n')
		line[len-1]=0;
	else
		return;
	len--;
	if(len && line[len-1]=='\r')
		line[len-1]=0;
}
bool split_string(char *str, const char *cut, vector<char *> &strs)
{
	char *p=str,*q;
	bool ret=false;
	while (q=strstr(p,cut))
	{
		*q=0;
		strs.push_back(p);
		p=q+strlen(cut);
		ret=true;
	}
	if(p) strs.push_back(p);
	return ret;
}



char* catch_string(char *str, char *head, char* tail, char* catched)
{
	char *p=str;
	char *q;
	q=strstr(str,head);
	if(!q)		return NULL;
	p=q+strlen(head);
	q=strstr(p,tail);
	if(!q)		return NULL;
	strncpy(catched,p,q-p);
	catched[q-p]=0;
	return q+strlen(tail);
}



char* catch_string(char *str, char* tail, char* catched)
// catch_string("12345","3",catched) => catched="12", return "45"
{
	char *q;
	q=strstr(str,tail);
	if(!q)
	{
		strcpy(catched,str);
		return NULL;
	}
	strncpy(catched,str,q-str);
	catched[q-str]=0;
	return q+strlen(tail);
}

bool build_table(ifstream &fin, vector<vector<string> >& table,int &lines)
{
	char line[MAXSTRLEN];
	table.clear();
	int i;
	bool not_end=true;
	while(not_end=(!fin.eof())){
		fin.getline(line,MAXSTRLEN-1);
		lines++;
		if(line[0]){
			
			vector<char *>columns;
			split_string(line,"\t",columns);
			vector<string> row(columns.size());
			for(i=0;i<columns.size();i++)
				row[i]=columns[i];	
			table.push_back(row);
		}else{
			return not_end;
		}
	}
	return not_end;
}


int vector_cmp(const vector<int> & s, const vector<int> & t) 
//-1 : s<t 1:s>t 0:s==t
{
	int size=s.size()>t.size()?t.size():s.size();
	int i;
	for(i=0;i<size;i++){
		int j1=s.size()-1-i;
		int j2=t.size()-1-i;
		if(s[j1]>t[j2])
			return 1;
		else if(s[j1]<t[j2])
			return -1;
	}
	if(s.size()<t.size())
		return -1;
	else if(s.size()>t.size())
		return 1;
	else
		return 0;
}

int get_feature_grid(int l, int p, int c, bool has_p, bool has_c){
	if(p<l){//not root
		if(!has_p){
			return c+l*l+l;
		}else if(!has_c){
			return l*l+p;
		}else{
			return p*l+c;
		}
	}else{
		return c*l+c;
	}
}

int get_lattice_grid(int i, int j, int l){
	//i->j
	if(i==l){
		return j*l+j;
	}else{
		return i*l+j;
	}
}