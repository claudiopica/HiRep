
/*******************************************************************************
*
* File check9.c
*
* Check of the program eva (random field)
*
* Author: Luigi Del Debbio
*
*******************************************************************************/

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
#include "eva.h"
#include "representation.h"
#include "global.h"

static int iw;
static float hmass=-0.15f;
static suNf_spinor **ws,**ev;
static double EPSILON=1.e-12;


#define MAX_ROTATE 50

static int init=0;
static const suNf_spinor s0={{{0.0f}}};
static suNf_spinor *psi,*chi;

static void alloc_ws_rotate(void)
{
   int n;

   psi=malloc(MAX_ROTATE*sizeof(suNf_spinor));
   chi=malloc(MAX_ROTATE*sizeof(suNf_spinor));

   error((psi==NULL)||(chi==NULL),1,"alloc_ws_rotate [linalg.c]",
         "Unable to allocate workspace");

   for (n=0;n<MAX_ROTATE;n++)
   {
      psi[n]=s0;
      chi[n]=s0;
   }

   init=1;
}

static void rotate(int vol,int n,suNf_spinor **ppk,complex v[])
{
   int k,j,ix;
   complex *z;
   suNf_spinor *pk,*pj;

   if (init==0)
      alloc_ws_rotate();

   error((n<1)||(n>MAX_ROTATE),1,"rotate [eva.c]",
         "Parameter n is out of range");

   for (ix=0;ix<vol;ix++)
   {
      for (k=0;k<n;k++)
      {
         pk=&(psi[k]);
         pj=ppk[0]+ix;
         z=&v[k];

         _vector_mulc_f((*pk).c1,*z,(*pj).c1);
         _vector_mulc_f((*pk).c2,*z,(*pj).c2);
         _vector_mulc_f((*pk).c3,*z,(*pj).c3);
         _vector_mulc_f((*pk).c4,*z,(*pj).c4);

         for (j=1;j<n;j++)
         {
            pj=ppk[j]+ix;
            z+=n;

            _vector_mulc_add_assign_f((*pk).c1,*z,(*pj).c1);
            _vector_mulc_add_assign_f((*pk).c2,*z,(*pj).c2);
            _vector_mulc_add_assign_f((*pk).c3,*z,(*pj).c3);
            _vector_mulc_add_assign_f((*pk).c4,*z,(*pj).c4);
         }
      }

      for (k=0;k<n;k++)
         *(ppk[k]+ix)=psi[k];
   }
}

static float normalize(suNf_spinor *ps)
{
   double r;

   r=spinor_field_sqnorm_f(ps);
   r=sqrt(r);
   error(r<EPSILON,1,"normalize [check9.c]","vector has vanishing norm");

   r=1.0/r;
   spinor_field_mul_f(ps,r,ps);

   return (float)(1.0/r);
}

static void eva_sort(int n,float d[],complex v[])
{
   int i,j,k;
   float p;
   complex q;
    
   for (i=0;i<n-1;i++)
   {
      k=i;
      p=d[i];
      
      for (j=i+1;j<n;j++)
      {
	 if (fabs((double)(d[j]))<fabs((double)(p)))
	 {
	    k=j;
	    p=d[j];
	 }
      }
      
      if (k!=i)
      {
	 d[k]=d[i];
	 d[i]=p;

	 for (j=0;j<n;j++)
	 {
	    q=v[n*j+i];
	    v[n*j+i]=v[n*j+k];
	    v[n*j+k]=q;
	 }         
      }
   }
}


static void eva_g5(int nev,float d[],suNf_spinor *ev[])
{
   int i,j;
   complex_dble z;
   complex *aa,*vv;

   aa=malloc(2*nev*nev*sizeof(complex));
   vv=aa+nev*nev;
   error(aa==NULL,1,"eva_Qnohat [check1.c]",
         "Unable to allocate auxiliary arrays");
   
   for (i=0;i<nev;i++)
   {
      g5Dphi(hmass,ev[iw],ev[i]);

      aa[nev*i+i].re=spinor_field_prod_re_f(ev[iw],ev[i]);
      aa[nev*i+i].im=0.0f;

      for (j=(i+1);j<nev;j++)
      {
         z=spinor_field_prod_f(ev[iw],ev[j]);

         aa[nev*i+j].re= (float)z.re;
         aa[nev*i+j].im= (float)z.im;
         aa[nev*j+i].re= (float)z.re;
         aa[nev*j+i].im=-(float)z.im;
      }
   }

   jacobi2(nev,aa,d,vv);
   eva_sort(nev,d,vv);
   rotate(VOLUME,nev,ev,vv);   
   free(aa);
}


static void Op1(suNf_spinor *out,suNf_spinor *in)
{
   g5Dphi(hmass,ev[iw],in);
   g5Dphi(hmass,out,ev[iw]);
}


static float power(int nit,spinor_operator Op,suNf_spinor *ws[])
{
   int i;
   float ubnd;

   gaussian_spinor_field(ws[0]);
   normalize(ws[0]);
   Op(ws[1],ws[0]);
   Op(ws[0],ws[1]);
   ubnd=normalize(ws[0]);   

   for (i=1;i<nit;i++)
   {
      Op(ws[1],ws[0]);
      Op(ws[0],ws[1]);
      ubnd=normalize(ws[0]);
   }
   return (float)sqrt((double)(ubnd));
}


int main(int argc,char *argv[])
{
   int i;
   int nev,nevt,ie,status;
   float omega1,omega2,d1[6],res,ubnd;
   complex_dble z;
   FILE *log=NULL;   

   log=freopen("check9.log","w",stdout);
   printf("\n");
   printf("Diagonalization of Q^2 (random fields)\n");
   printf("--------------------------------------\n\n");

   printf("The lattice size is %dx%d^3\n",T,L);

   rlxs_init(0,12345);

   geometry_eo_lexi();
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();
   
   printf("Generating a random gauge field... ");fflush(stdout);
   random_u();
   printf("done.\n");
   represent_gauge_field();
   
   set_spinor_len(VOLUME);

   ws=malloc(2*sizeof(suNf_spinor*));
   for (i=0;i<2;i++)
      ws[i]=alloc_spinor_field_f();
   ev=malloc(7*sizeof(suNf_spinor*));
   for (i=0;i<7;i++)
      ev[i]=alloc_spinor_field_f();
   
   iw=6;
   nev=4;
   nevt=6;
   ubnd=1.05f*power(20,Op1,ws);
   printf("test-ubnd: %f\n",ubnd);
   omega1=1.0e-6f;
   omega2=1.0e-3f;

   printf("Accuracy parameters: omega1=%.1e, omega2=%.1e\n\n",
          omega1,omega2);

   ie=eva(VOLUME,nev,nevt,0,100,20,1,ubnd,omega1,omega2,Op1,ws,ev,d1,&status);

   printf("\nEigenvalues of Q^2 (status = %d, ie = %d):\n\n",
          status,ie);

   for (i=0;i<6;i++)
   {
      Op1(ws[0],ev[i]);
      z.re=-(double)d1[i];
      z.im=(double)0.0f;
      spinor_field_mulc_add_assign_f(ws[0],z,ev[i]);
      res=spinor_field_sqnorm_f(ws[0]);
      res=(float)(sqrt((double)(res)));

      if (i==4)
         printf("\n");
      printf("d[%d] = % .3e, acc = %.1e\n",i,d1[i],res);
   }

   eva_g5(4,d1,ev);
      
   printf("\n");
   printf("Eigenvalues of Q:\n\n");

   for (i=0;i<4;i++)
   {
      g5Dphi(hmass,ws[0],ev[i]);      
      z.re=-d1[i];
      z.im=0.0f;
      spinor_field_mulc_add_assign_f(ws[0],z,ev[i]);
      res=spinor_field_sqnorm_f(ws[0]);
      res=(float)(sqrt((double)(res)));
      
      printf("d[%d] = % .3e, acc = %.1e\n",i,d1[i],res);
   }

   printf("\n");
   fclose(log);

   exit(0);
}

