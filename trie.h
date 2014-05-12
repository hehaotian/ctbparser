/////////////////////////////////////////////
//NOTICE:  TrieNode->wordptr&&Trie->link=NULL;
//All word in trie must be GB2312
/////////////////////////////////////////////
#ifndef TRIE_H
#define TRIE_H
const int WORD_MAXLENGTH = 1000;
const int POS_MAXLENGTH = 12;
inline int CC_INT(char c)
{
	return (unsigned char)c>0?(unsigned char)c-0xA0:0;
}

#define CC_NUM 6768
//pos freq pair
typedef struct PosFreqChain{
	char pos[POS_MAXLENGTH];
	double freq;
	PosFreqChain *next;
}POSFREQCHAIN, *PPOSFREQCHAIN;

typedef struct WordType{
	char* word;
	double freq;//=sum of posfreq;
	PPOSFREQCHAIN posfreq;
}WORDTYPE, *PWORDTYPE;



typedef struct TrieNode{
	int k;//no null link number
	TrieNode **link; //if no branch,link=NULL else *link[95]
	PWORDTYPE wordptr; //TrieNode->wordptr&&Trie->link=NULL;TrieNode->wordptr||Trie->link!=NULL
}TRIENODE, * PTRIENODE;

class Trie {
private:
	void PrintNode(PTRIENODE,char*);//print every elem in trie in depth first
	void PrintWordNode(PWORDTYPE,char*);


	void PrintNode2File(PTRIENODE,char*);//print every elem in trie in depth first
	void PrintWordNode2File(PWORDTYPE,char*);

public:

	bool DeleteWord (char * word , char *pos=0 );//pos=NULL:delete word regardless pos
	bool DeleteNode (PTRIENODE p);//free *p
	PTRIENODE root;//root of trie tree
	Trie();
	~Trie();
	void Clear();
	bool Load(char *sFilename);//construct trie
	PWORDTYPE Insert (char * word , char *pos, double freq);
	PWORDTYPE Search (char * word);
	bool IsExist(char * word, char *pos=0);
	PTRIENODE NearSearch (char * word, PTRIENODE start);
	PTRIENODE NearSearch (char * word);
	double GetFreq (char * word , char *pos=0 );//pos=NULL:get all freq of word regardless pos
	void Print(char*);//print every elem in trie in depth first
	void Print2File(char*);//print every elem in trie in depth first
};
#endif