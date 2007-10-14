
/******************************************************************************
*
* File check11.c
*
* Consistency checks on the programs in the module linalg
*
* Author: luigi del debbio <luigi.del.debbio@ed.ac.uk>
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
#include "update.h"
#include "linear_algebra.h"
#include "representation.h"
#include "global.h"

#define MAX_ROTATE 50

static complex v[25];
static double EPSILON=1.e-12;
static suNf_spinor *ppk[5];

static int nsp,initr=0;
static const suNf_spinor s0={{{{{0.0f}}}}};
static suNf_spinor *psi,*chi;

static void lc1(double c1,suNf_spinor *ps1,suNf_spinor *ps2)
{
   suNf_spinor *psm;

   psm=ps1+nsp;

   for (;ps1<psm;ps1++)
   {
      _spinor_mul_add_assign_f(*ps1,c1,*ps2);
      ps2+=1;
   }
}


static void lc2(double c1,double c2,suNf_spinor *ps1,suNf_spinor *ps2)
{
   suNf_spinor *psm;

   psm=ps1+nsp;

   for (;ps1<psm;ps1++)
   {
      _spinor_lc_f(*ps1,c1,*ps1,c2,*ps2);
      ps2+=1;
   }
}

static void lc3(double c1,double c2,suNf_spinor *ps1,suNf_spinor *ps2,suNf_spinor *ps3)
{
   suNf_spinor *psm;

   psm=ps1+nsp;

   c1=-c1; c2=-c2;
   for (;ps1<psm;ps1++)
   {
      _spinor_lc_add_assign_f(*ps3,c1,*ps1,c2,*ps2);
      _spinor_minus_f(*ps3,*ps3);
      ps2+=1;
      ps3+=1;
   }
}

static void alloc_ws_rotate(void)
{
   psi=calloc(MAX_ROTATE,sizeof(suNf_spinor));
   chi=calloc(MAX_ROTATE,sizeof(suNf_spinor));

   error((psi==NULL)||(chi==NULL),1,"alloc_ws_rotate [linalg.c]",
         "Unable to allocate workspace");

   initr=1;
}

static void rotate(int vol,int n,suNf_spinor **pkk,complex v[])
{
   int k,j,ix;
   complex *z;
   suNf_spinor *pk,*pj;

   if (initr==0)
      alloc_ws_rotate();

   error((n<1)||(n>MAX_ROTATE),1,"rotate [eva.c]",
         "Parameter n is out of range");

   for (ix=0;ix<vol;ix++)
   {
      for (k=0;k<n;k++)
      {
         pk=&(psi[k]);
         pj=pkk[0]+ix;
         z=&v[k];

				 _spinor_mulc_f(*pk,*z,*pj);

         for (j=1;j<n;j++)
         {
            pj=pkk[j]+ix;
            z+=n;

						_spinor_mulc_add_assign_f(*pk,*z,*pj);
         }
      }

      for (k=0;k<n;k++)
         *(pkk[k]+ix)=psi[k];
   }
}

static void project(suNf_spinor *pk,suNf_spinor *pl)
{
   complex sp;

   sp.re=-spinor_field_prod_re_f(pl,pk);
   sp.im=-spinor_field_prod_im_f(pl,pk);

   spinor_field_mulc_add_assign_f(pk,sp,pl);
}   

static double normalize(suNf_spinor *ps)
{
   double r,ri;

   r=spinor_field_sqnorm_f(ps);
   r=sqrt(r);
   error(r<EPSILON,1,"normalize [eva.c]","vector has vanishing norm");

   ri=1.0/r;
   spinor_field_mul_f(ps,ri,ps);

   return (double)(r);
}

static complex sp(int vol,suNf_spinor *pk,suNf_spinor *pl)
{
   int ix;
   double x,y;
   complex *rpk,*rpl,z;

   x=0.0;
   y=0.0;

   rpk=(complex*)(pk);
   rpl=(complex*)(pl);
   
   for (ix=0;ix<(4*NF*vol);ix++)
   {
      x+=(double)((*rpk).re*(*rpl).re+(*rpk).im*(*rpl).im);
      y+=(double)((*rpk).re*(*rpl).im-(*rpk).im*(*rpl).re);
      rpk+=1;
      rpl+=1;
   }
   
   z.re=(double)(x);
   z.im=(double)(y);
   
   return z;
}


int main(int argc,char *argv[])
{
   int i,j,vol=VOLUME,off=0;
   double r;
   double rd,zsqd;
   double d,dmax;
   complex w;
   complex zd,wd;
   suNf_spinor **ws;
   suNf_spinor *pk,*pl;
   suNf_spinor tmp[VOLUME];
   FILE *log=NULL;   

   log=freopen("check11.log","w",stdout);
      
   printf("\n");
   printf("Consistency of the programs in the module linalg\n");
   printf("------------------------------------------------\n\n");   
   printf("The lattice size is %dx%d^3\n",T,L);

   nsp=VOLUME;
   
   rlxs_init(0,12345);

   geometry_eo_lexi();
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();

   set_spinor_len(VOLUME);

   ws=malloc(10*sizeof(suNf_spinor*));
   for (i=0;i<10;i++)
      ws[i]=alloc_spinor_field_f(1);


   for (i=0;i<10;i++)
      gaussian_spinor_field(&(ws[i][0]));

   dmax=0.0;
   
   for (i=0;i<10;i++)
   {
      pk=&ws[i][off];
      pl=&ws[9-i][off];
      w=sp(vol,pk,pl);
      
      zd=spinor_field_prod_f(pk,pl);
      rd=spinor_field_sqnorm_f(pk)*spinor_field_sqnorm_f(pl);
      d=((zd.re-(double)w.re)*(zd.re-(double)w.re)+
         (zd.im-(double)w.im)*(zd.im-(double)w.im));
      d=sqrt(d/rd);
      if (d>dmax)
         dmax=d;

      rd=spinor_field_prod_re_f(pk,pl);
      d=fabs(zd.re/rd-1.0);
      if (d>dmax)
         dmax=d;

      zd=spinor_field_prod_f(pk,pk);
      rd=spinor_field_sqnorm_f(pk);
      
      d=fabs(zd.im/rd);
      if (d>dmax)
         dmax=d;

      d=fabs(zd.re/rd-1.0f);
      if (d>dmax)
         dmax=d;
   }

   printf("Check of spinor_field_prod, spinor_field_prod_re\n");
   printf("and spinor_field_sqnorm: %.2e\n\n",dmax);
   
   dmax=0.0;
   zd.re= 0.345;
   zd.im=-0.876;
   zsqd=zd.re*zd.re+zd.im*zd.im;
   
   for (i=0;i<9;i++)
   {
      pk=&ws[i][off];
      pl=&ws[i+1][off];
      
      wd=spinor_field_prod_f(pk,pl);
      rd=spinor_field_sqnorm_f(pk)+zsqd*spinor_field_sqnorm_f(pl)
         +2.0*(zd.re*wd.re-zd.im*wd.im);
      
      spinor_field_mulc_add_assign_f(pk,zd,pl);

      d=fabs(rd/spinor_field_sqnorm_f(pk)-1.0);
      if (d>dmax)
         dmax=d;
   }

   printf("Consistency of spinor_prod, norm_square\n");
   printf("and mulc_spinor_add: %.2e\n\n",dmax);
   
   for (i=0;i<10;i++)
      gaussian_spinor_field(&(ws[i][0]));

   dmax=0.0;
   
   for (i=0;i<10;i++)
   {
      pk=&ws[i][off];
      
      if (i>0)
      {
         pl=&ws[i-1][off];
         project(pk,pl);
         zd=spinor_field_prod_f(pk,pl);
         
         d=(fabs(zd.re)+
            fabs(zd.im))/
            sqrt(spinor_field_sqnorm_f(pk));
         
         if (d>dmax)
            dmax=d;
      }
      
      normalize(pk);
      rd=spinor_field_sqnorm_f(pk);
      
      d=fabs(rd-1.0f);
      if (d>dmax)
         dmax=d;
   }
   
   printf("Consistency of spinor_prod, norm_square,\n");
   printf("normalize and project: %.2e\n\n",dmax);
   
   for (i=0;i<5;i++)
   {
      pk=&ws[i][off];
      pl=&ws[i+5][off];
      
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      
      for (j=0;j<5;j++)
      {
         v[5*i+j].re=0.1234f*(double)(i^2)-0.8976f*(double)(j);
         v[5*i+j].im=0.2231f*(double)(i)+0.9922f*(double)(j^2);
      }

      ppk[i]=pl;
   }
   
   rotate(vol,5,ppk,v);
   dmax=0.0;
   
   for (i=5;i<10;i++)
   {
      pk=&ws[i][off];
      
      for (j=0;j<5;j++)
      {
         zd.re=-(double)v[5*j+(i-5)].re;
         zd.im=-(double)v[5*j+(i-5)].im;
         
         pl=&ws[j][off];
         spinor_field_mulc_add_assign_f(pk,zd,pl);
      }
      
      rd=spinor_field_sqnorm_f(pk);
      
      d=fabs(rd);
      if (d>dmax)
         dmax=d;
   }
   
   dmax/=spinor_field_sqnorm_f(&ws[0][off]);
   dmax=sqrt(dmax);
   
   printf("Consistency of mulc_spinor_add\n");
   printf("and rotate: %.2e\n\n",dmax);
   
   dmax=0.0;
   
   for (i=0;i<5;i++)
   {
      pk=&ws[i][off];
      pl=&ws[9-i][off];
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      spinor_field_g5_f(tmp,pk);
      spinor_field_g5_f(pk,tmp);

      zd.re=-1.0;
      zd.im=0.0;

      spinor_field_mulc_add_assign_f(pl,zd,pk);
      r=spinor_field_sqnorm_f(pl)/spinor_field_sqnorm_f(pk);
      d=sqrt(r);
      if (d>dmax)
         dmax=d;

      gaussian_spinor_field(pl);
      zd=spinor_field_prod_f(pk,pl);
      spinor_field_g5_f(pk,pk);
      spinor_field_g5_f(pl,pl);
      wd=spinor_field_prod_f(pk,pl);
      
      d=(fabs(zd.re-wd.re)+fabs(zd.im-wd.im))/
                  (fabs(zd.re)+fabs(zd.im));
      if (d>dmax)
         dmax=d;
   }

   printf("Check of spinor_field_g5_f: %.2e\n\n",dmax);
   
   dmax=0.0;
   
   for (i=0;i<5;i++)
   {
      pk=&ws[i][off];
      pl=&ws[9-i][off];
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      d=-2.5;
      lc1(d,pk,pl);

      zd.re=1.5;
      zd.im=0.0;
      spinor_field_mulc_add_assign_f(pk,zd,pl);
      d=spinor_field_sqnorm_f(pk)/spinor_field_sqnorm_f(pl);
      
      if (d>dmax)
         dmax=d;
   }

   printf("Check of lc1: %.2e\n\n",dmax);
   
   dmax=0.0;
   
   for (i=0;i<5;i++)
   {
      pk=&ws[i][off];
      pl=&ws[9-i][off];
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      d=1.0;
      r=2.5;
      lc2(d,r,pk,pl);

      zd.re=-3.5;
      zd.im=0.0;
      spinor_field_mulc_add_assign_f(pk,zd,pl);
      d=spinor_field_sqnorm_f(pk)/spinor_field_sqnorm_f(pl);
      
      if (d>dmax)
         dmax=d;
   }

   printf("Check of lc2: %.2e\n\n",dmax);
   
   dmax=0.0;
   
   for (i=0;i<5;i++)
   {
      pk=&ws[i][off];
      pl=&ws[9-i][off];
      gaussian_spinor_field(pk);
      spinor_field_copy_f(pl,pk);
      d=3.5;
      r=-1.5;
      lc3(d,r,pk,pl,pk);

      zd.re=-1.0;
      zd.im=0.0;
      spinor_field_mulc_add_assign_f(pk,zd,pl);
      d=spinor_field_sqnorm_f(pk)/spinor_field_sqnorm_f(pl);
      
      if (d>dmax)
         dmax=d;
   }

   printf("Check of lc3: %.2e\n\n",dmax);
   
   fclose(log);
   exit(0);
}
