#include "trie.h"
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <stdio.h>

bool Trie::Load(char *sFilename)
{
	FILE *fp;
	
	//读取每个词条
	if((fp=fopen(sFilename,"r"))==NULL)
		return false;//fail while opening the file
	double freq;
	char pos[POS_MAXLENGTH];
	char word[WORD_MAXLENGTH];
	while(fscanf(fp,"%s\t%s\t%lf",word,pos,&freq)!=EOF)
		Insert(word,pos,freq);
	fclose(fp);

	return true;
}
bool Trie::DeleteNode (PTRIENODE p)//free p
{
	PWORDTYPE w;
	PPOSFREQCHAIN pf,pf1;
	int i;
	if(!p)	return false;
	if(p->link)
	{
		for(i=0;i<95;i++) 
		{
			DeleteNode(p->link[i]);
			p->link[i]=NULL;
		}
		free(p->link);
		p->link=NULL;
	}else{
		w=p->wordptr;
		for(pf=w->posfreq;pf;)
		{
			pf1=pf;
			pf=pf->next;
			free(pf1);
			pf1=NULL;
		}
		free(w->word);
		free(w);
		p->wordptr=NULL;
	}
	free(p);
	return true;
}



PWORDTYPE Trie::Insert ( char * word , char *pos, double freq)
{
	PWORDTYPE w;
	PTRIENODE p=root;
	PTRIENODE pp;
	PTRIENODE t;
	PPOSFREQCHAIN pf,pf1;
	int i;
	if(!p)//initialization
	{
		pf=(PPOSFREQCHAIN)malloc(sizeof(POSFREQCHAIN));
		pf->freq=freq;
		strcpy(pf->pos,pos);
		pf->next=NULL;
		w=(PWORDTYPE)malloc(sizeof(WORDTYPE));
		w->posfreq=pf;
		w->freq=freq;
		w->word=(char *)malloc(sizeof(char)*(strlen(word)+1));
		strcpy(w->word,word);
		t=(PTRIENODE)malloc(sizeof(TRIENODE));
		t->link=NULL;
		t->wordptr=w;
		t->k=0;
		root=t;
		return w;
	}
	for(i=0;p&&p->link;i++)
	{
		pp=p;
		p=p->link[CC_INT(word[i])];
	}
	if(!p)//p=NULL
	{
		p=pp;//since pp->link!=NULL , pp->wordptr should be NULL;
		pf=(PPOSFREQCHAIN)malloc(sizeof(POSFREQCHAIN));
		pf->freq=freq;
		strcpy(pf->pos,pos);
		pf->next=NULL;
		w=(PWORDTYPE)malloc(sizeof(WORDTYPE));
		w->posfreq=pf;
		w->freq=freq;
		w->word=(char *)malloc(sizeof(char)*(strlen(word)+1));
		strcpy(w->word,word);
		t=(PTRIENODE)malloc(sizeof(TRIENODE));
		t->link=NULL;
		t->wordptr=w;
		t->k=0;
		p->link[CC_INT(word[i-1])]=t;
		p->k++;
	}else{//p->link=NULL,=>p->wordptr!=NULL
		w=p->wordptr;
		if(!strcmp(w->word,word))//word is already exists
		{
			pf1=NULL;
			for(pf=w->posfreq;pf&&strcmp(pf->pos,pos)<0;pf1=pf,pf=pf->next);
			if(pf&&!strcmp(pf->pos,pos))
			{
				pf->freq+=freq;
				w->freq+=freq;
				return w;//find
			}
			pf=(PPOSFREQCHAIN)malloc(sizeof(POSFREQCHAIN));
			pf->freq=freq;
			strcpy(pf->pos,pos);
			if(!pf1)//insert as w->posfreq
			{
				pf->next=w->posfreq;
				w->posfreq=pf;
			}else{
				pf->next=pf1->next;
				pf1->next=pf;
			}
			w->freq+=freq;
		}else{//not exists, create new nodes
			p->wordptr=NULL;//set p->wordptr=NULL, and alloc p->link
			while(w->word[i]==word[i])//insert links
			{
				p->link=(PTRIENODE *)malloc(95*sizeof(PTRIENODE));//alloc 95 points to p->link
				memset(p->link,NULL,sizeof(PTRIENODE)*95);//initial
				t=(PTRIENODE)malloc(sizeof(TRIENODE));
				t->link=NULL;
				t->wordptr=NULL;
				t->k=0;
				p->link[CC_INT(word[i])]=t;
				p->k++;
				p=t;
				i++;
			}
			//insert elements
			p->link=(PTRIENODE *)malloc(95*sizeof(PTRIENODE));//alloc 95 points to p->link
			memset(p->link,NULL,sizeof(PTRIENODE)*95);//initial

			t=(PTRIENODE)malloc(sizeof(TRIENODE));
			t->link=NULL;
			t->wordptr=w;
			t->k=0;
			p->link[CC_INT(w->word[i])]=t;//MEMORY LEAK
			p->k++;

			pf=(PPOSFREQCHAIN)malloc(sizeof(POSFREQCHAIN));
			pf->freq=freq;
			strcpy(pf->pos,pos);
			pf->next=NULL;
			w=(PWORDTYPE)malloc(sizeof(WORDTYPE));
			w->posfreq=pf;
			w->freq=freq;
			w->word=(char *)malloc(sizeof(char)*(strlen(word)+1));
			strcpy(w->word,word);
			t=(PTRIENODE)malloc(sizeof(TRIENODE));
			t->link=NULL;
			t->wordptr=w;
			t->k=0;
			p->link[CC_INT(word[i])]=t;
			p->k++;
		}
	}
	return w;
}


void Trie::Print(char* s)
{
	s[0]=0;
	PrintNode(root,s);
}
void Trie::PrintNode(PTRIENODE p,char* s)
{
	if(!p)
		return;
	if(p->link)	for(int i=0;i<95;i++)PrintNode(p->link[i],s);
	else
		PrintWordNode(p->wordptr,s);
}
void Trie::PrintWordNode(PWORDTYPE p,char* s)
{
	if(!p)
		return;
	
	char temp[WORD_MAXLENGTH+20];
	sprintf(temp,"%s/%lf\t",p->word,p->freq);
	strcat(s,temp);
}

void Trie::Print2File(char* s)
{
	PrintNode2File(root,s);
}

void Trie::PrintNode2File(PTRIENODE p,char* s)
{
	if(!p)
		return;
	if(p->link)	for(int i=0;i<95;i++)PrintNode2File(p->link[i],s);
	else
		PrintWordNode2File(p->wordptr,s);
}
void Trie::PrintWordNode2File(PWORDTYPE p,char* s)
{
	if(!p)
		return;
	
	FILE *fp;
	PPOSFREQCHAIN pf;
	if((fp=fopen(s,"a"))==NULL)
		return;
	for(pf=p->posfreq;pf;pf=pf->next)
		fprintf(fp,"%s\t%s\t%lf\n",p->word,pf->pos,(double)pf->freq);
	fclose(fp);
}

Trie::Trie()
{
	root=NULL;
}



Trie::~Trie()
{
	Clear();
}

void Trie::Clear()
{
	DeleteNode(root);
	root=NULL;
}

bool Trie::DeleteWord (char * word , char *pos )
{
	PWORDTYPE w;
	PTRIENODE p=root;
	PTRIENODE pp=NULL;
	PTRIENODE ppp=root;//delete from ppp to p
	PPOSFREQCHAIN pf,pf1;
	int i;
	for(i=0;p&&p->link;i++)
	{
		if(p->k>1) ppp=p->link[CC_INT(word[i])];
		pp=p;
		p=p->link[CC_INT(word[i])];
	}
	if(!p) return false;
	w=p->wordptr;//p->link=NULL => p->wordptr!=NULL
	if(strcmp(w->word,word)) return false;
	if(pos==NULL)//delete all w->posfreq
	{
		for(pf=w->posfreq;pf;)
		{
			pf1=pf;
			pf=pf->next;
			free(pf1);
			pf1=NULL;
		}
	}else{//delete the wordpos->pos==pos
		for(pf=w->posfreq;pf&&strcmp(pf->pos,pos)<0;pf1=pf,pf=pf->next);
		if(pf&&!strcmp(pf->pos,pos))
		{//find
			pf1->next=pf->next;
			free(pf);
			pf=NULL;
		}else
			return false;
	}
	if(!w->posfreq)//no posfreq, delete w
	{	
		DeleteNode(ppp);
		ppp=NULL;
	}
	return true;
}

PWORDTYPE Trie::Search (char * word)//word must be ended with '\0'
{
	PTRIENODE p=root;
	int i;
	for(i=0;p&&p->link;i++) p=p->link[CC_INT(word[i])];
	if(!p) return NULL;//not find
	if(!strcmp(p->wordptr->word,word)) return p->wordptr;//find
	return NULL;
}

PTRIENODE Trie::NearSearch (char * word)//word must be ended with '\0'
{
	PTRIENODE p=root;
	int len=strlen(word);
	int i;
	for(i=0;i<len&&p&&p->link;i++) p=p->link[CC_INT(word[i])];
	return p;//if !p return NULL,if !p->link return the elem, if p->link,i>len return the branch
}

PTRIENODE Trie::NearSearch (char * word, PTRIENODE start)//word must be ended with '\0'
{
	PTRIENODE p=start;
	int len=strlen(word);
	int i;
	for(i=0;i<len&&p&&p->link;i++) 		p=p->link[CC_INT(word[i])];
	return p;//if !p return NULL,if !p->link return the elem, if p->link,i>len return the branch
}
double Trie::GetFreq (char * word , char *pos )
{
	PWORDTYPE w=Search (word);
	PPOSFREQCHAIN pf;
	if(!w) return 0;
	if(pos==NULL)//get the sum of freq
	{
		return w->freq;
	}else{
		for(pf=w->posfreq;pf&&strcmp(pf->pos,pos)<0;pf=pf->next);
		if(pf&&!strcmp(pf->pos,pos))
		{//find
			return pf->freq;
		}
	}
	return 0;
}
bool Trie::IsExist(char * word , char *pos)
{
	PWORDTYPE w=Search (word);
	PPOSFREQCHAIN pf;
	if(!w) return false;
	if(pos==NULL)return true;
	for(pf=w->posfreq;pf&&strcmp(pf->pos,pos)<0;pf=pf->next);
	if(pf&&!strcmp(pf->pos,pos))
		return true;
	return false;
}

