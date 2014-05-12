#include <fstream>
#include <iostream>
#include <ctime>
#include "ctbparser.h"
using namespace std;
int main()
{
	ctbparser *c=new ctbparser();
	if(!c->load_config("config.txt")){
		delete c;
		c=0;
		return 1;
	}
	char s[10000];
	char t[10000];
	cout<<"请输入一句话：(输入'Q'退出)"<<endl;
	do{
		cin>>s;
		if(!strcmp(s,"Q"))
			break;
		c->decode_string(s,t);
		cout<<t<<endl;
	}while(1);
	c->decode_file("in.txt","out.txt");
	cout<<"ctbparser test success!"<<endl;
	delete c;
	c=0;
	return 0;
}
