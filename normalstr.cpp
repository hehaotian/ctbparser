
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <fstream>
#include <string>
#include <string.h>
#include "normalstr.h"
#include "const.h"
using namespace std;
bool normalstr::certify(char * s)//word must be ended with '\0'
{
	int j;
    for(j=0;(unsigned char)s[j]>=0xA0&&(unsigned char)s[j]<=0xFE;j++);
    if(s[j]) return false;
    return true;
}

bool normalstr::load(char *sFilename)
{
	ifstream fin;
	int i,j;
	char c;
	fin.open(sFilename);
	if(!fin.is_open())
		return false;
	for(i=0;i<CODENUM;i++){
		j=0;
		while((c=fin.get())!='\n') 
			code[i][j++]=c;
		code[i][j]=0;
		assert(certify(code[i]));
	}
	fin.close();
	return true;
}

void normalstr::convert_string(char *in, char *out)
{
	int i,j,k,c1,c2;
	out[0]=0;
	j=0;
	for(i=0;in[i];)
	{
		if(in[i]<0)//word
		{
			c1=(unsigned char)in[i];
			c2=(unsigned char)in[i+1];
			if(c1>=0xB0&&c1<=0xF7&&c2>=0xA0&&c2<=0xFE)//GB2312
			{
				strncat(out,in+i,2);
			}else if(c1>=0xA1&&c1<=0xA9&&c2>=0xA0&&c2<=0xFE)//GB2312
			{	
				k=(c1-0xA1)*(0xFF-0xA0);
				k+=c2-0xA0+127;
				strcat(out,code[k]);
			}else if(c1<=0xA0&&c1>=0x81&&c2>=0x40&&c2<=0xFE&&c2!=0x7F)//GBK/3
			{
				k=(c1-0x81)*190;
				if(c2>0x7F)
					k=k+c2-0x40-1;
				else
					k=k+c2-0x40;
				k=k+(0xA9-0xA1+1)*(0xFF-0xA0)+127;//skip over GB2312
				strcat(out,code[k]);
			}else if(c1>=0xAA&&c1<=0xFE&&c2>=0x40&&c2<=0xA0&&c2!=0x7F)//GBK/4
			{
				k=(c1-0xAA)*96;
				if(c2>0x7F)
					k+=c2-0x40-1;
				else
					k+=c2-0x40;
				k+=(0xA9-0xA1+1)*(0xFF-0xA0)+(0xA0-0x81+1)*190+127;
				strcat(out,code[k]);
			}else if(c1>=0xA8&&c1<=0xA9&&c2>=0x40&&c2<=0xA0&&c2!=0x7F){//GBK/5
				k=(c1-0xA8)*96;
				if(c2>0x7F)
					k+=c2-0x40-1;
				else
					k+=c2-0x40;
				k+=(0xA9-0xA1+1)*(0xFF-0xA0)+(0xA0-0x81+1)*190+85*6*16+127;
				strcat(out,code[k]);
			}//else strcat(out "");
			i+=2;
		}
		else{
			strcat(out,code[in[i]-1]);
			i++;
		}
	}
}

void normalstr::inv_convert(char *gb, char *point, char *out)
{
	out[0]=0;
	char *p;
	char *q;
	char c[3];
	char w[5];
	int wl;
	c[2]=0;
	for(p=gb,q=point;*p;)
	{
		if(*p<0)
		{
			if(*q>0)
			{
				c[0]=*q;
				c[1]=0;
				q++;
				convert_string(c,w);
				strcat(out,c);
				wl=strlen(w);
				p+=wl;
			}else{
				c[0]=*q;
				c[1]=*(q+1);
				q+=2;
				convert_string(c,w);
				strcat(out,c);
				wl=strlen(w);
				p+=wl;
			}
		}else{
			c[0]=*p;
			c[1]=0;
			strcat(out,c);
			p++;
		}
	}
	if(*q)
	{
		strcat(out,q);
		strcat(out,"/w  ");
	}
}


normalstr::normalstr()
{

}

normalstr::~normalstr()
{
}

bool normalstr::convert_file(char *fi, char* fo)
{
	ifstream fin;
	ofstream fout;
	fin.open(fi);
	if(!fin.is_open())
		return false;
	fout.open(fo);
	if(!fout.is_open())
	{
		fin.close();
		return false;//fail while opening the file
	}
	char c[3]={0};
	char C[3];
	while(!fin.eof()){
		c[0]=fin.get();
		c[1]=0;
		if(c[0]<0){
			c[1]=fin.get();
		}
		if(!c[1]&& (c[0]=='\r' || c[0]=='\n'))//
			strcpy(C,c);
		else
			convert_string(c,C);
		fout<<C;
	}
	fin.close();
	fout.close();
	return true;
}
