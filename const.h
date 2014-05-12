#ifndef CONST_H
#define CONST_H
const int PAGESIZE = 8192;
const int MAXSTRLEN = 131072;//16 pages
const int MAXSENLEN = 1024;//1024 characters
const int SIMPLE_CHAIN=0;
const int FIRST_CHAIN=1;//first order markov chain
const int GENERAL_CHAIN=2;

const int CRF_ALGORITHM=0;
const int L1_ALGORITHM=3;
const int AP_ALGORITHM=1;
const int PA_ALGORITHM=2;
const int FACTOR_ALGORITHM=4;

const int SPEED_PRIOR=0;
const int MEM_PRIOR=1;
const int DIRECT_LEARN=0;
const int REFINED_LEARN=1;
const double EPS=1e-15;

const int TASK_SEG=0;
const int TASK_POS=1;
const int TASK_PARSE=2;

#define SENSEP "¡££¿£¡"

const int INF_INT=1024*1024*1024;
const int MAXWORDLEN=1024;//1024 characters
const double INF=1e+37;
const double LOG_INF=37;
#endif
