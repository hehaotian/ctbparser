#include <math.h>
#include <string.h>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#include "lbfgs.h"

const double ETA=1e-7;
const double XTOL=1e-7;
const int MP = 6;
const int LP = 6;
const double GTOL = .9;
const double STPMIN = 1e-20;
const double STPMAX = 1e20;
const double FTOL=1e-4;
const int MAXFEV = 20;
const double P5=.5;
const double P66=.66;
const double XTRAPF=4;


inline double inner_product(int n,double *a, double *b){
	int i;
	double s=0;
	for(i=0;i<n;i++)
		s+=a[i]*b[i];
	return s;
}
inline void daxpy(int n,double *a, double b,double *c){//a+=b*c
	int i;
	for(i=0;i<n;i++)
		a[i]+=b*c[i];
}
LBFGS::~LBFGS()
{
	if(diag)
		delete [] diag;
	if(w)
		delete [] w;
	char fn[1000];
	int i;
	for(i=0;i<m;i++)
	{
		sprintf(fn,"__s_%d",i);
		if(!access(fn,0))//file exists
			unlink(fn);
		sprintf(fn,"__y_%d",i);
		if(!access(fn,0))//file exists
			unlink(fn);
	}
	if(!access("__x",0))
		unlink("__x");
	if(!access("__q",0))
		unlink("__q");
	if(!access("__w0",0))
		unlink("__w0");
	if(!access("__w1",0))
		unlink("__w1");

}
void LBFGS::save(double *buf, char *file_name, int index)
{
	char fn[1000];
	if(index>=0)
		sprintf(fn,"%s_%d",file_name,index);
	else
		strcpy(fn,file_name);
	FILE *fp=fopen(fn,"wb");
	fwrite(buf,sizeof(double),n,fp);

	fclose(fp);
}
void LBFGS::load(double *buf, char *file_name, int index)
{
	char fn[1000];
	if(index>=0)
		sprintf(fn,"%s_%d",file_name,index);
	else
		strcpy(fn,file_name);
	FILE *fp=fopen(fn,"rb");
	fread(buf,sizeof(double),n,fp);
	fclose(fp);
}
bool LBFGS::init(int parameter_num, int depth, int pri){
	int i;
	prior=pri;
	n = parameter_num;
	m = depth;
	s.resize(m);
	y.resize(m);
	if(prior==0)
	{
		try
		{
			w=new double[n + 2*m*n];
			if(!w)
				return false;
			for(i=0;i<m;i++)
			{
				s[i]=&w[n+i*n];
				y[i]=&w[n+m*n+i*n];
			}
			q=&w[0];
			diag=new double[n];
			if(!diag)
				return false;
		}catch(bad_alloc&){
			return false;
		}
	}
	alpha.resize(m);
	rho.resize(m);

	if(!s.size() || !y.size() || !alpha.size() || !rho.size())
		return false;
	
	return true;
}




/* Subroutine */ static int mcstep_(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy, double *stp, 
									double *fp, double *dp, int *brackt,
     double *stpmin, double *stpmax, int *info)
{
  /* System generated locals */
  double d__1, d__2, d__3;



  /* Local variables */
  static double sgnd, stpc, stpf, stpq, p, q, gamma, r__, s, theta;
  static int bound;

  *info = 0;

  /*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

  if (*brackt && ((*stp <= min(*stx,*sty) || *stp >= max(*stx,*sty)) || *dx *
                  (*stp - *stx) >=(float)0. || *stpmax < *stpmin)) {
    return 0;
  }

  /*     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN. */

  sgnd = *dp *(*dx / fabs(*dx));

  /*     FIRST CASE. A HIGHER FUNCTION VALUE. */
  /*     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER */
  /*     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN, */
  /*     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN. */

  if (*fp > *fx) {
    *info = 1;
    bound = 1;
    theta =(*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
    /* Computing MAX */
    d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = max(d__1,d__2), d__2 = fabs(
                                                                          *dp);
    s = max(d__1,d__2);
    /* Computing 2nd power */
    d__1 = theta / s;
    gamma = s * sqrt(d__1 * d__1 - *dx / s *(*dp / s));
    if (*stp < *stx) {
      gamma = -gamma;
    }
    p = gamma - *dx + theta;
    q = gamma - *dx + gamma + *dp;
    r__ = p / q;
    stpc = *stx + r__ *(*stp - *stx);
    stpq = *stx + *dx /((*fx - *fp) / (*stp - *stx) + *dx) / 2 * (*stp -
                                                                  *stx);
    if ((d__1 = stpc - *stx, fabs(d__1)) < (d__2 = stpq - *stx, fabs(d__2)))
      {
        stpf = stpc;
      } else {
        stpf = stpc + (stpq - stpc) / 2;
      }
    *brackt = 1;

    /*     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF */
    /*     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC */
    /*     STEP IS CLOSER TO STX THAN THE QUADRATIC(SECANT) STEP, */
    /*     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN. */

  } else if (sgnd < (float)0.) {
    *info = 2;
    bound = 0;
    theta =(*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
    /* Computing MAX */
    d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = max(d__1,d__2), d__2 = fabs(
                                                                          *dp);
    s = max(d__1,d__2);
    /* Computing 2nd power */
    d__1 = theta / s;
    gamma = s * sqrt(d__1 * d__1 - *dx / s *(*dp / s));
    if (*stp > *stx) {
      gamma = -gamma;
    }
    p = gamma - *dp + theta;
    q = gamma - *dp + gamma + *dx;
    r__ = p / q;
    stpc = *stp + r__ *(*stx - *stp);
    stpq = *stp + *dp /(*dp - *dx) * (*stx - *stp);
    if ((d__1 = stpc - *stp, fabs(d__1)) > (d__2 = stpq - *stp, fabs(d__2)))
      {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
    *brackt = 1;

    /*     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE */
    /*     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES. */
    /*     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY */
    /*     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC */
    /*     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE */
    /*     EITHER STPMIN OR STPMAX. THE QUADRATIC(SECANT) STEP IS ALSO */
    /*     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP */
    /*     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN. */

  } else if (fabs(*dp) < fabs(*dx)) {
    *info = 3;
    bound = 1;
    theta =(*fx - *fp) * 3 / (*stp - *stx) + *dx + *dp;
    /* Computing MAX */
    d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = max(d__1,d__2), d__2 = fabs(
                                                                          *dp);
    s = max(d__1,d__2);

    /*        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND */
    /*        TO INFINITY IN THE DIRECTION OF THE STEP. */

    /* Computing MAX */
    /* Computing 2nd power */
    d__3 = theta / s;
    d__1 = 0., d__2 = d__3 * d__3 - *dx / s *(*dp / s);
    gamma = s * sqrt((max(d__1,d__2)));
    if (*stp > *stx) {
      gamma = -gamma;
    }
    p = gamma - *dp + theta;
    q = gamma + (*dx - *dp) + gamma;
    r__ = p / q;
    if (r__ < (float)0. && gamma != (float)0.) {
      stpc = *stp + r__ *(*stx - *stp);
    } else if (*stp > *stx) {
      stpc = *stpmax;
    } else {
      stpc = *stpmin;
    }
    stpq = *stp + *dp /(*dp - *dx) * (*stx - *stp);
    if (*brackt) {
      if ((d__1 = *stp - stpc, fabs(d__1)) < (d__2 = *stp - stpq, fabs(
                                                                     d__2))) {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
    } else {
      if ((d__1 = *stp - stpc, fabs(d__1)) > (d__2 = *stp - stpq, fabs(
                                                                     d__2))) {
        stpf = stpc;
      } else {
        stpf = stpq;
      }
    }

    /*     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE */
    /*     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES */
    /*     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP */
    /*     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN. */

  } else {
    *info = 4;
    bound = 0;
    if (*brackt) {
      theta =(*fp - *fy) * 3 / (*sty - *stp) + *dy + *dp;
      /* Computing MAX */
      d__1 = fabs(theta), d__2 = fabs(*dy), d__1 = max(d__1,d__2), d__2 =
        fabs(*dp);
      s = max(d__1,d__2);
      /* Computing 2nd power */
      d__1 = theta / s;
      gamma = s * sqrt(d__1 * d__1 - *dy / s *(*dp / s));
      if (*stp > *sty) {
        gamma = -gamma;
      }
      p = gamma - *dp + theta;
      q = gamma - *dp + gamma + *dy;
      r__ = p / q;
      stpc = *stp + r__ *(*sty - *stp);
      stpf = stpc;
    } else if (*stp > *stx) {
      stpf = *stpmax;
    } else {
      stpf = *stpmin;
    }
  }

  /*     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT */
  /*     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE. */

  if (*fp > *fx) {
    *sty = *stp;
    *fy = *fp;
    *dy = *dp;
  } else {
    if (sgnd < (float)0.) {
      *sty = *stx;
      *fy = *fx;
      *dy = *dx;
    }
    *stx = *stp;
    *fx = *fp;
    *dx = *dp;
  }

  /*     COMPUTE THE NEW STEP AND SAFEGUARD IT. */

  stpf = min(*stpmax,stpf);
  stpf = max(*stpmin,stpf);
  *stp = stpf;
  if (*brackt && bound) {
    if (*sty > *stx) {
      /* Computing MIN */
      d__1 = *stx + (*sty - *stx) * (float).66;
      *stp = min(d__1,*stp);
    } else {
      /* Computing MAX */
      d__1 = *stx + (*sty - *stx) * (float).66;
      *stp = max(d__1,*stp);
    }
  }
  return 0;

  /*     LAST LINE OF SUBROUTINE MCSTEP. */

} /* mcstep_ */



static int mcsrch_(int *n, double *x, double *f, double *g, double *s, double* stp, int *info, int *nfev, double *wa, int orthant)
{

  /* System generated locals */
  int i__1;
  double d__1;

  /* Local variables */
  static double dgxm, dgym;
  static int j, infoc;
  static double finit, width, stmin, stmax;
  static int stage1;
  static double width1, ftest1, dg, fm, fx, fy;
  static int brackt;
  static double dginit, dgtest;
  static double dgm, dgx, dgy, fxm, fym, stx, sty;

  /* Parameter adjustments */
  --wa;
  --s;
  --g;
  --x;


  /* Function Body */
  if (*info == -1) {
    goto L45;
  }
  infoc = 1;

  /*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

  if (*n <= 0 || *stp <= 0) {
    return 0;
  }

  /*     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION */
  /*     AND CHECK THAT S IS A DESCENT DIRECTION. */

  dginit = 0;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    dginit += g[j] * s[j];
    /* L10: */
  }
  if (dginit >= 0) {
    return 0;
  }

  /*     INITIALIZE LOCAL VARIABLES. */

  brackt = 0;
  stage1 = 1;
  *nfev = 0;
  finit = *f;
  dgtest = FTOL * dginit;
  width = STPMAX - STPMIN;
  width1 = width / P5;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    wa[j] = x[j];
    /* L20: */
  }

  /*     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP, */
  /*     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP. */
  /*     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP, */
  /*     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF */
  /*     THE INTERVAL OF UNCERTAINTY. */
  /*     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP, */
  /*     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP. */

  stx = 0;
  fx = finit;
  dgx = dginit;
  sty = 0;
  fy = finit;
  dgy = dginit;

  /*     START OF ITERATION. */

 L30:

  /*        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND */
  /*        TO THE PRESENT INTERVAL OF UNCERTAINTY. */

  if (brackt) {
    stmin = min(stx,sty);
    stmax = max(stx,sty);
  } else {
    stmin = stx;
    stmax = *stp + XTRAPF *(*stp - stx);
  }

  /*        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN. */

  *stp = max(*stp,STPMIN);
  *stp = min(*stp,STPMAX);

  /*        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET */
  /*        STP BE THE LOWEST POINT OBTAINED SO FAR. */

  if ((brackt && ((*stp <= stmin || *stp >= stmax) || *nfev >= MAXFEV - 1 ||
                  infoc == 0)) ||(brackt && (stmax - stmin <= XTOL * stmax))) {
    *stp = stx;
  }

  /*        EVALUATE THE FUNCTION AND GRADIENT AT STP */
  /*        AND COMPUTE THE DIRECTIONAL DERIVATIVE. */
  /*        We return to main program to obtain F and G. */

  i__1 = *n;
  if(orthant)
  {
	for (j = 1; j <= i__1; ++j)
	{
		if(g[j]==0 && x[j]==0)
			continue;
		//stop at 0 if x[j] cross 0
		x[j] = wa[j] + *stp * s[j];
		if(wa[j]*x[j]<0)
			x[j]=0;
	}
  }else{
	  for (j = 1; j <= i__1; ++j) {
		x[j] = wa[j] + *stp * s[j];
		/* L40: */
	  }
  }
  *info = -1;
  return 0;

 L45:
  *info = 0;
  ++(*nfev);
  dg = 0;
  i__1 = *n;
  for (j = 1; j <= i__1; ++j) {
    dg += g[j] * s[j];
    /* L50: */
  }
  ftest1 = finit + *stp * dgtest;

  /*        TEST FOR CONVERGENCE. */

  if (brackt && ((*stp <= stmin || *stp >= stmax) || infoc == 0)) {
    *info = 6;
  }
  if (*stp == STPMAX && *f <= ftest1 && dg <= dgtest) {
    *info = 5;
  }
  if (*stp == STPMIN && (*f > ftest1 || dg >= dgtest)) {
    *info = 4;
  }
  if (*nfev >= MAXFEV) {
    *info = 3;
  }
  if (brackt && stmax - stmin <= XTOL * stmax) {
    *info = 2;
  }
  if (*f <= ftest1 && fabs(dg) <= GTOL * (-dginit)) {
    *info = 1;
  }

  /*        CHECK FOR TERMINATION. */

  if (*info != 0) {
    return 0;
  }

  /*        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED */
  /*        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE. */

  if (stage1 && *f <= ftest1 && dg >= min(FTOL,GTOL) * dginit) {
    stage1 = 0;
  }

  /*        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF */
  /*        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED */
  /*        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE */
  /*        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN */
  /*        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT. */

  if (stage1 && *f <= fx && *f > ftest1) {

    /*           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES. */

    fm = *f - *stp * dgtest;
    fxm = fx - stx * dgtest;
    fym = fy - sty * dgtest;
    dgm = dg - dgtest;
    dgxm = dgx - dgtest;
    dgym = dgy - dgtest;

    /*           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY */
    /*           AND TO COMPUTE THE NEW STEP. */

    mcstep_(&stx, &fxm, &dgxm, &sty, &fym, &dgym, stp, &fm, &dgm, &brackt,
            &stmin, &stmax, &infoc);

    /*           RESET THE FUNCTION AND GRADIENT VALUES FOR F. */

    fx = fxm + stx * dgtest;
    fy = fym + sty * dgtest;
    dgx = dgxm + dgtest;
    dgy = dgym + dgtest;
  } else {

    /*           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY */
    /*           AND TO COMPUTE THE NEW STEP. */

    mcstep_(&stx, &fx, &dgx, &sty, &fy, &dgy, stp, f, &dg, &brackt, &
            stmin, &stmax, &infoc);
  }

  /*        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE */
  /*        INTERVAL OF UNCERTAINTY. */

  if (brackt) {
    if ((d__1 = sty - stx, fabs(d__1)) >= P66 * width1) {
      *stp = stx + P5 *(sty - stx);
    }
    width1 = width;
    width =(d__1 = sty - stx, fabs(d__1));
  }

  /*        END OF ITERATION. */

  goto L30;

  /*     LAST LINE OF SUBROUTINE MCSRCH. */

} /* mcsrch_ */


int LBFGS::optimize(double *x, double *f, double *g, int orthant, double *w0, double *w1)
{
	int i,j,k;
	int index;
	double ys,yy,h_0;
	double *w2,*w3;//4 work spaces for memory saving prior
	if(prior==1)
	{
		w2=x;
		w3=g;
	}else if(prior==0){
		w1=&diag[0];
	}
	double gnorm;
	if(iflag==0)//first iteration
	{
		if(prior==1)
		{
			s[0]=w0;
			q=w1;
		}
		for(i=0;i<n;i++)
			s[0][i]=-g[i];
		gnorm=sqrt(inner_product(n,g,g));
		stp1=1/gnorm;
	}else{
		if(prior==1)
		{
			load(w0,"__w0");
			load(w1,"__w1");
			s[iter%m]=w0;
		}
		goto L172;
	}
	while(1)
	{//main loop
		//w0:null, w1:q, w2:null, w3:g
		//for first iteration
		//w0:s,w1:q,w2:x,w3:g
		info=0;
		bound=iter>m?m:iter;
		if(iter==0)
			goto L165;
		index=(iter+m-1)%m;
		if(prior==1)
		{
			s[index]=w0;
			y[index]=w2;
			load(s[index],"__s",index);
			load(y[index],"__y",index);
		}

		//w0:s, w1:q, w2:y, w3:g
		ys=inner_product(n,y[index],s[index]);
		yy=inner_product(n,y[index],y[index]);
		h_0=ys/yy;
		rho[index]=1/ys;
		for(i=0;i<n;i++)
			q[i]=-g[i];
		index=iter%m;
		for(i=0;i<bound;i++)
		{
			index=(index+m-1)%m;
			if(prior==1 && i>0)
			{
				s[index]=w0;
				y[index]=w2;
				load(s[index],"__s",index);
				load(y[index],"__y",index);
			}
			double sq=inner_product(n,q,s[index]);
			alpha[index]=sq*rho[index];
			daxpy(n,q,-alpha[index],y[index]);
		}
		for(i=0;i<n;i++)
			q[i]*=h_0;
		for(i=0;i<bound;i++)
		{
			if(prior==1 && i>0)
			{
				s[index]=w0;
				y[index]=w2;
				load(s[index],"__s",index);
				load(y[index],"__y",index);
			}
			double yr=inner_product(n,y[index],q);
			double beta=yr*rho[index];
			beta=alpha[index]-beta;
			daxpy(n,q,beta,s[index]);
			index=(index+1)%m;
		}
		if(prior==1)
			s[iter%m]=w0;
		//w0:s, w1:q, w2:y, w3:g
		for(i=0;i<n;i++)
			s[iter%m][i]=q[i];
		if(prior==1)
		{
			x=w2;
			load(x,"__x");
		}
L165://w0:s, w1:q, w2:x, w3:g
		nfev = 0;
		stp = 1;
		if (iter == 0)
			stp = stp1;
		for(i=0;i<n;i++)
			q[i] = g[i];
		if(prior==1)
			save(q,"__q");
L172://w0:s, w1:temp, w2:x, w3:g
		mcsrch_(&n, x, f, g, s[iter%m], &stp, &info, &nfev, w1, orthant);
		if (info == -1)
		{
			if(prior==1)
			{
				save(w0,"__w0");
				save(w1,"__w1");
			}
			iflag = 1;
			return 1;
		}
		if (info != 1)
			return -1;//error, see iflag
		gnorm=sqrt(inner_product(n,g,g));
		double xnorm=sqrt(inner_product(n,x,x));
		if(prior==1){
			q=w1;
			load(q,"__q");
			save(x,"__x");
			y[iter%m]=w2;
		}
		//w0:s, w1:q, w2:y, w3:g
		for(i=0;i<n;i++)
		{
			s[iter%m][i]*=stp;
			y[iter%m][i]=g[i]-q[i];
		}
		if(xnorm<1) xnorm=1;
		if (gnorm / xnorm <= ETA) {
			finish = true;
		}
		if(finish)
		{
			iflag=0;
			return 0;
		}
		if(prior==1)
		{
			save(y[iter%m],"__y",iter%m);
			save(s[iter%m],"__s",iter%m);
		}
		//w0:null, w1:q, w2:null, w3:g
		iter++;
	}
}
