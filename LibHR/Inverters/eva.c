/*******************************************************************************
*
* File eva.c
*
* Generic program that serves to compute the low-lying eigenvalues and the
* associated eigenvectors of an arbitrary hermitian operator Op acting on
* single-precision spinor fields, using a subspace iteration algorithm
*
* The program for the operator must be externally defined and is assumed to
* have the following properties:
*
*   spinor_operator Op(spinor_field *pk,spinor_field *pl)
*     Application of the operator Op to the field pk[] and assignement of the
*     result to pl[]. On exit the field pk[] is unchanged
*
* The externally accessible function is
*
*   int eva(int nev,int nevt,int init,int kmax,
*           int imax,float ubnd,float omega1,float omega2,
*           spinor_operator Op,
*           spinor_field *ev,float d[],int *status)
*     Computes the lowest eigenvalues d[0],...,d[nevt-1] of the operator Op
*     and the corresponding eigenvectors ev[0],..,ev[nevt-1]. The first
*     nev of the eigenvalues are obtained to an absolute precision omega1 or
*     a relative precision omega2 (whichever is reached first), while the
*     other eigenvalues may be less accurate. On exit the program returns 0
*     if the eigenvalues are obtained to the desired precision or a negative
*     value if this was not possible. The other parameters are
*
*     init     Specifies whether all eigenvectors should be initialized
*              (init=0), or only the last nevt-nev eigenvectors (init=1)
*              or none of them (init!=0 or 1)
*
*     kmax     Maximal degree of the Chebyshev polynomials used for the 
*              acceleration of the algorithm
*
*     imax     Maximal number of subspace iterations
*
*     ubnd     Upper bound on the eigenvalues of the operator
*
*     status   On exit this variable reports the number of times the
*              operator was applied
*
* Notes:
*
* The algorithm implemented in this program is a subspace iteration with
* Chebyshev acceleration and eigenvector locking. See
*
*   Y. Saad: Numerical methods for large eigenvalue problems
*
* for example (the book was printed in 1991 by Manchester University Press
* and can be downloaded from http://www.cs.umn.edu/saad/books.html)
*
* Error bounds on the first nev eigenvalues are deduced from the residues
* (or the residue matrices in the case of closely spaced eigenvalues) of
* the associated approximate eigenvectors. This method is not optimal but
* completely safe even if there are (approximate) degeneracies in the
* spectrum. The accuracy of the other eigenvalues is not monitored and
* may be rather poor on exit, particularly so in the case of the last
* (nevt-nev)/2 eigenvalues
*
* The spectrum of the operator must be guaranteed to be below the specified
* upper bound. Otherwise the algorithm diverges and the program terminates
* abnormally. However, setting the parameter ubnd to a value far above the
* upper end of the spectrum will unduly slow down the algorithm
*
* This program is adapted from Martin Luscher's original version, as used
* for the study of the eigenvalues of the Wilson-Dirac operator presented
* in hep-lat/0512021
*
* Author: Luigi Del Debbio <luigi.del.debbio@ed.ac.uk>
*	
*	Modified: Claudio Pica 
*	Modified: Agostino Patella 
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "error.h"
#include "update.h"
#include "complex.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "logger.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "dirac.h"

#define GAMMA 3.0
#define MAX_ROTATE 1000 /*50*/

static int nop,nvc=0;
static double *dd,*ee;
static complex *aa,*bb,*cc,*vv;

static double EPSILON=1.e-12;

static suNf_spinor *psi=NULL;

static void alloc_ws_rotate(void) {
  psi=calloc(MAX_ROTATE,sizeof(suNf_spinor));
  
  error((psi==NULL),1,"alloc_ws_rotate [linalg.c]",
        "Unable to allocate workspace");
  
}

static void rotate(int n,spinor_field *pkk,complex v[])
{

  error((n<1)||(n>MAX_ROTATE),1,"rotate [eva.c]",
        "Parameter n is out of range");

  if (psi==NULL) alloc_ws_rotate();
  
  //Avoid OMP parallel region in MASTER_FOR
#undef _OMP_PRAGMA
#define _OMP_PRAGMA(s)
  
  _MASTER_FOR(pkk->type,ix) {
    for (int k=0;k<n;k++) {
      suNf_spinor *pk=&(psi[k]);
      suNf_spinor *pj=_FIELD_AT(&pkk[0],ix);
      complex *z=&v[k];
      
      _vector_mulc_f((*pk).c[0],*z,(*pj).c[0]);
      _vector_mulc_f((*pk).c[1],*z,(*pj).c[1]);
      _vector_mulc_f((*pk).c[2],*z,(*pj).c[2]);
      _vector_mulc_f((*pk).c[3],*z,(*pj).c[3]);
      
      for (int j=1;j<n;j++) {
        pj=_FIELD_AT(&pkk[j],ix);
        z+=n;
        
        _vector_mulc_add_assign_f((*pk).c[0],*z,(*pj).c[0]);
        _vector_mulc_add_assign_f((*pk).c[1],*z,(*pj).c[1]);
        _vector_mulc_add_assign_f((*pk).c[2],*z,(*pj).c[2]);
        _vector_mulc_add_assign_f((*pk).c[3],*z,(*pj).c[3]);
      }
    }
    
    for (int k=0;k<n;k++) {
      *_FIELD_AT(&pkk[k],ix)=psi[k];
    }
  }
}

static int alloc_aux(int nevt)
{
  if (nevt>nvc)
    {
      if (nvc>0)
	{
	  free(aa);
	  free(dd);
	}
      
      aa=malloc(4*nevt*nevt*sizeof(complex));
      dd=malloc(2*nevt*sizeof(double));

      bb=aa+nevt*nevt;
      cc=bb+nevt*nevt;
      vv=cc+nevt*nevt;
      ee=dd+nevt;
      nvc=nevt;
    }
   
  error((aa==NULL)||(dd==NULL),1,"alloc_aux [eva.c]",
	"Unable to allocate auxiliary arrays");
  return 0;
}


static void project(spinor_field *pk,spinor_field *pl)
{
  complex sp;

  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(pk,pl);
#endif

  sp.re=-spinor_field_prod_re_f(pl,pk);
  sp.im=-spinor_field_prod_im_f(pl,pk);

  spinor_field_mulc_add_assign_f(pk,sp,pl);
}   


static double normalize(spinor_field *ps)
{
  double r;

  r=sqrt(spinor_field_sqnorm_f(ps));
  error(r<EPSILON,1,"normalize [eva.c]","vector has vanishing norm");

  spinor_field_mul_f(ps,1./r,ps);

  apply_BCs_on_spinor_field(ps);

  return r;
}


static void mgm_subsp(int n, spinor_field *ev)
{
  int k;
   
  for (k=0;k<n;k++)
    project(&ev[n],&ev[k]);

  normalize(&ev[n]);
}


static void init_subsp(int nev,int nevt,int init,spinor_field *ev)
{
  int n;

  for (n=0;n<nevt;n++)
    {

      if ((init==0)||((init==1)&&(n>=nev)))
	gaussian_spinor_field(&ev[n]);

      mgm_subsp(n,ev);
    }
}


static void ritz_subsp(int nlock,int nevt,spinor_operator Op,
                       spinor_field *ws,spinor_field *ev,double d[])
{
  int neff,i,j;
  complex z;

  neff=nevt-nlock;
   
  for (i=0;i<neff;i++) 
    {
      Op(&ws[0],&ev[nlock+i]);
      nop+=1;      
      
      aa[neff*i+i].re=spinor_field_prod_re_f(&ev[nlock+i],&ws[0]);
      aa[neff*i+i].im=0.0f;

      for (j=0;j<i;j++) 
	{
	  z=spinor_field_prod_f(&ws[0],&ev[nlock+j]);

	  aa[neff*i+j].re= z.re;
	  aa[neff*i+j].im= z.im;
	  aa[neff*j+i].re= z.re;
	  aa[neff*j+i].im=-z.im;
	}
    }

  jacobi2(neff,aa,d+nlock,vv);
  rotate(neff,ev+nlock,vv);
}


static void submat(int nev,int ia,int ib)
{
  int i,j,n;

  n=0;
   
  for (i=ia;i<ib;i++)
    {
      for (j=ia;j<ib;j++)
	{
	  cc[n]=bb[i*nev+j];
	  n+=1;
	}
    }
}


static double min_eva(int ia,int ib,double d[])
{
  int i;
  double r0,r1;

  r0=fabs(d[ia]);

  for (i=ia+1;i<ib;i++)
    {
      r1=fabs(d[i]);
      if (r1<r0)
	r0=r1;
    }

  return r0;
}


static int res_subsp(int nlock,int nev,double omega1,double omega2,
                     spinor_operator Op,
                     spinor_field *ws,spinor_field *ev,double d[])
{
  int i,ia,ib;
  double eps1,eps2,absd1,absd2;
  complex z;

  eps1=0.0f;
  ia=nlock;
   
  for (ib=nlock;ib<nev;ib++) 
    {
      if (ib>ia)
	{
	  submat(nev,ia,ib);
	  jacobi2(ib-ia,cc,dd,vv);
	  eps1=sqrt(dd[ib-ia-1]);
	  absd1=min_eva(ia,ib,d);

	  for (i=ia;i<ib;i++)
            ee[i]=eps1;

	  if ((eps1>omega1)&&(eps1>(omega2*absd1)))
            return ia;
	}

      Op(&ws[0],&ev[ib]);
      nop+=1;
      spinor_field_lc1_f(-d[ib],&ws[0],&ev[ib]);

      bb[nev*ib+ib].re=spinor_field_sqnorm_f(&ws[0]);
      bb[nev*ib+ib].im=0.0f;

      eps2=sqrt(bb[nev*ib+ib].re);
      absd2=fabs(d[ib]);
      ee[ib]=eps2;
      
      if (ib>ia)
	{
	  if ((d[ib]-d[ib-1])>(eps1+eps2))
            ia=ib;
	}
          
      if ((eps2>omega1)&&(eps2>(omega2*absd2)))
	return ia;

      if (ib>ia)
	{
	  Op(&ws[1],&ws[0]);
	  nop+=1;      

	  for (i=ia;i<ib;i++) 
	    {
	      z=spinor_field_prod_f(&ws[1],&ev[i]);

	      bb[nev*ib+i].re= z.re;
	      bb[nev*ib+i].im= z.im;
	      bb[nev*i+ib].re= z.re;
	      bb[nev*i+ib].im=-z.im;
	    }
	}
    }

  submat(nev,ia,ib);
  jacobi2(ib-ia,cc,dd,vv);
  eps1=sqrt(dd[ib-ia-1]);
  absd1=min_eva(ia,ib,d);

  for (i=ia;i<ib;i++)
    ee[i]=eps1;   
      
  if ((eps1>omega1)&&(eps1>(omega2*absd1)))
    return ia;
   
  return nev;
}


static double set_lbnd(int nevt,int kmax,double ubnd,double d[],int *k)
{
  double mu1,mu2,nu1,nu2,t;

  mu1=d[0];
  mu2=d[nevt-1];
   
  nu1=GAMMA*GAMMA;
  nu1=log(nu1+sqrt(nu1*nu1-1.0));

  nu2=GAMMA;
  nu2=log(nu2+sqrt(nu2*nu2-1.0));
   
  (*k)=(int)(0.5*sqrt(((ubnd-mu1)*nu1*nu1-(ubnd-mu2)*nu2*nu2)/(mu2-mu1)));
  (*k)&=~0x1;

  if ((*k)>kmax)
    (*k)=(kmax&~0x1);
  if ((*k)<2)
    (*k)=2;

  t=tanh(nu2/(double)(2*(*k)));

  return mu2+(ubnd-mu2)*t*t;
}


static void apply_cheby(int k,double lbnd,double ubnd,
                        spinor_operator Op,
                        spinor_field *ws,spinor_field *ev)
{
  int j;
  double c1,c2;
  spinor_field *psi0,*psi1,*psi2,*psi3;

  c1=2.0f/(ubnd-lbnd);
  c2=-(ubnd+lbnd)/(ubnd-lbnd);

  psi0=ev;
  psi1=&ws[0];
  psi2=&ws[1];

  Op(psi1,psi0);
  spinor_field_lc2_f(c1,c2,psi1,psi0);

  c1*=2.0f;
  c2*=2.0f;

  for (j=1;j<k;j++)
    {
      Op(psi2,psi1);
      spinor_field_lc3_f(c1,c2,psi2,psi1,psi0);

      psi3=psi0;
      psi0=psi1;
      psi1=psi3;
    }

  nop+=k;
}


int eva(int nev,int nevt,int init,int kmax,
        int imax,double ubnd,double omega1,double omega2,
        spinor_operator Op,
        spinor_field *ev,double d[],int *status)   
{
  int i,k,n;
  int nlock,nupd,nlst;
  double lbnd;
  spinor_field *ws;
  
  *status=0;
   
  error((nev<=0)||(nevt<nev)||(kmax<2)||(imax<1),1,
	"eva [eva.c]",
	"Improper parameters nev,nevt,kmax or imax");
  error((omega1<(ubnd*EPSILON))&&(omega2<EPSILON),1,
	"eva [eva.c]",
	"Improper parameters omega1 or omega2");

  nop=0;
  nlock=0;
  nupd=(nevt+nev)/2;
   
  if (alloc_aux(nevt)!=0)
    return -3;
  
  ws=alloc_spinor_field_f(2,ev->type);
  
  init_subsp(nev,nupd,init,ev);
  ritz_subsp(nlock,nupd,Op,ws,ev,d);
   
  for (i=0;i<imax;i++){
    if (i>0){
      lbnd=set_lbnd(nupd,kmax,ubnd,d,&k);
      for (n=nupd;n<nevt;n++)
        spinor_field_copy_f(&ev[n],&ev[n+nupd-nevt]);
    } else
      lbnd=set_lbnd(nupd,10,ubnd,d,&k);
    for (n=nlock;n<nevt;n++) {
      if (n<nupd){
        apply_cheby(k,lbnd,ubnd,Op,ws,&ev[n]);
        mgm_subsp(n,ev);
      }else if (i>0)
        mgm_subsp(n,ev);
    }
    
    if (i>0)
      nlst=nevt;
    else
      nlst=nupd;
    
    ritz_subsp(nlock,nlst,Op,ws,ev,d);
    n=nlock;      
    nlock=res_subsp(n,nev,omega1,omega2,Op,ws,ev,d);
    *status=nop;
    
    lprintf("EVA",10,"i=%3d, k=%3d, d[%2d]=% .6e, d[%2d]=% .6e, lbnd=% .6e, eps[%2d]=% .1e\n",
      i,k,n,d[n],nlst-1,d[nlst-1],lbnd,n,ee[n]);
    
    error(d[nlst-1]>ubnd,1,"eva [eva.c]",
          "Parameter ubnd is too low");

    if (nlock==nev){
      lprintf("EVA",10,"Computation succeded. MVM = %d\n",*status);
      free_spinor_field_f(ws);
      return 0;
    }
  }

  lprintf("EVA",10,"Unable to reach required precision. MVM = %d\n",*status);
  
  free_spinor_field_f(ws);
  
  return -1;
}

int eva_tuned(int nev,int nevt,int init,int kmax,
        int imax,double lbnd, double ubnd,double omega1,double omega2,
        spinor_operator Op,
        spinor_field *ev,double d[],int *status)   
{
  int i,k,n;
  int nlock,nupd,nlst;
  //double lbnd;
  spinor_field *ws;
  

  *status=0;
   
  error((nev<=0)||(nevt<nev)||(kmax<2)||(imax<1),1,
	"eva [eva.c]",
	"Improper parameters nev,nevt,kmax or imax");
  error((omega1<(ubnd*EPSILON))&&(omega2<EPSILON),1,
	"eva [eva.c]",
	"Improper parameters omega1 or omega2");

  nop=0;
  nlock=0;
  nupd=(nevt+nev)/2;
   
  if (alloc_aux(nevt)!=0)
    return -3;
  
  ws=alloc_spinor_field_f(2,ev->type);
    
  init_subsp(nev,nupd,init,ev);
  ritz_subsp(nlock,nupd,Op,ws,ev,d);
   
  for (i=0;i<imax;i++){
    if (i>0){
//      lbnd=set_lbnd(nupd,kmax,ubnd,d,&k);
k = kmax;
      for (n=nupd;n<nevt;n++)
        spinor_field_copy_f(&ev[n],&ev[n+nupd-nevt]);
    } else
//      lbnd=set_lbnd(nupd,10,ubnd,d,&k);
k = kmax;
    for (n=nlock;n<nevt;n++) {
      if (n<nupd){
        apply_cheby(k,lbnd,ubnd,Op,ws,&ev[n]);
        mgm_subsp(n,ev);
      }else if (i>0)
        mgm_subsp(n,ev);
    }
    
    if (i>0)
      nlst=nevt;
    else
      nlst=nupd;
    
    ritz_subsp(nlock,nlst,Op,ws,ev,d);
    n=nlock;      
    nlock=res_subsp(n,nev,omega1,omega2,Op,ws,ev,d);
    *status=nop;
    
    lprintf("EVA",10,"i=%3d, k=%3d, d[%2d]=% .6e, d[%2d]=% .6e, lbnd=% .6e, eps[%2d]=% .1e\n",
      i,k,n,d[n],nlst-1,d[nlst-1],lbnd,n,ee[n]);
    
    error(d[nlst-1]>ubnd,1,"eva [eva.c]",
          "Parameter ubnd is too low");

    if (nlock==nev){
      lprintf("EVA",10,"Computation succeded. MVM = %d\n",*status);
      free_spinor_field_f(ws);
      return 0;
    }
  }

  lprintf("EVA",10,"Unable to reach required precision. MVM = %d\n",*status);
  
  free_spinor_field_f(ws);
  
  return -1;
}
