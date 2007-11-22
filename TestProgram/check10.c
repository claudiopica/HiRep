
/******************************************************************************
*
* File check10.c
*
* Check of the program eva (free case)
*
* Author: Luigi Del Debbio <luigi.del.debbio@ed.ac.uk>
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "linear_algebra.h"
#include "update.h"
#include "inverters.h"
#include "dirac.h"
#include "representation.h"
#include "global.h"
#include "logger.h"

static int iw;
static double ms;
static spinor_field *ws;
static double EPSILON=1.e-12;


static double low_ev(void)
{
   int lmax;
   double p,sp,cp,pi;

   lmax=T;   
   if (X>lmax)
      lmax=X;
   if (Y>lmax)
      lmax=Y;
   if (Z>lmax)
      lmax=Z;

   pi=4.0*atan(1.0);
   p=(2.0*pi)/(double)(lmax);
   sp=sin(p);
   cp=ms+1.0-cos(p);

   return (double)(sp*sp+cp*cp);
}


static double normalize(spinor_field *ps)
{
   double r;

   r=spinor_field_sqnorm_f(ps);
   r=sqrt(r);
   error(r<EPSILON,1,"normalize [check9.c]","vector has vanishing norm");

   r=1.0/r;
   spinor_field_mul_f(ps,r,ps);

   return (double)(1.0/r);
}


static void Op1(spinor_field *pl,spinor_field *pk)
{
   g5Dphi(ms,&ws[iw],pk);
   g5Dphi(ms,pl,&ws[iw]);
}


static double power(int nit,spinor_operator Op,spinor_field *ws)
{
   int i;
   double ubnd;

   gaussian_spinor_field(&ws[0]);
   normalize(&ws[0]);
   Op(&ws[1],&ws[0]);
   Op(&ws[0],&ws[1]);
   ubnd=normalize(&ws[0]);   

   for (i=1;i<nit;i++)
   {
      Op(&ws[1],&ws[0]);
      Op(&ws[0],&ws[1]);
      ubnd=normalize(&ws[0]);
   }

   return (double)sqrt((double)(ubnd));
}


int main(int argc,char *argv[])
{
   int i;
   int nev,nevt,ie,status;
   double ubnd,omega1,omega2,d1[16],res,lev;
   spinor_field *ev;
   FILE *log=NULL;   

   log=freopen("check10.log","w",stdout);
   printf("\n");
   printf("Diagonalization of Qnohat^2 (free case)\n");
   printf("---------------------------------------\n\n");
   geometry_eo_lexi();
   printf("The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
   printf("size of the gluon rep: %d, size of the fermion rep: %d\n",NG,NF);
   
   rlxd_init(1,12345);

	 logger_setlevel(0,1000);

   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();

   set_spinor_len(VOLUME);
   ws=alloc_spinor_field_f(3);
   ev=alloc_spinor_field_f(4*NF+8);

   
   iw=2;

   ms=0.25;
   
   nev=4*NF+4;
   nevt=nev+4;
   ubnd=1.05f*power(30,Op1,ws);   
   printf("ubnd-test: %f\n",ubnd);
   omega1=1.0e-6f;
   omega2=1.0e-3f;
   lev=low_ev();

   printf("Accuracy parameters: omega1=%.1e, omega2=%.1e\n\n",
          omega1,omega2);

   ie=eva(VOLUME,nev,nevt,0,100,20,ubnd,omega1,omega2,Op1,ws,ev,d1,&status);
   
   printf("\nEigenvalues of Qnohat^2 (status = %d, ie = %d):\n\n",status,ie);

   for (i=0;i<(4*NF+8);i++)
   {
      if (i<4*NF)
         res=(double)fabs((double)(d1[i]-ms*ms));
      else
         res=(double)fabs((double)(d1[i]-lev));

      if (i==nev)
         printf("\n");
      printf("d[%2d] = % .3e acc = %1.5e\n",i,d1[i],res);
   }

   printf("\n");
   fclose(log);
   
   exit(0);
}

