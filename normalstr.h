#ifndef NORMALSTR_H
#define NORMALSTR_H
const int CODENUM=15414;//(0xA9-0xA1+1)*(0xFF-0xA0)+(0xA0-0x81+1)*(12*16-2)+(0xFE-0xAA+1)*6*16+2*(6*16+1)
class normalstr
{
public:
	normalstr();
	~normalstr();
	void convert_string(char *in, char *out);//convert GBK coded "in" to GB2312 coded "out", neither ASC II nor GBK is allowed in "out"
	void inv_convert(char *gb,char *point,char *out);//gb: the gb coded string, convert gb to out according to point, and move point to new position
	bool load(char *sFilename);//Load convert dict
	bool convert_file(char *fin, char* fout);
	bool certify(char * s);//certify s is GB2312 coded
private:
	char code[CODENUM][5];
};
#endif