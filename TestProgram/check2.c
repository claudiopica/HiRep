/*******************************************************************************
*
* Action of the Dirac operator on plane waves
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"

double hmass=0.1;

static suNf_spinor mul_gamma(int mu,suNf_spinor s)
{
   suNf_spinor r;
   complex i,m_i,m_1;

   i.re=0.0f;
   i.im=1.0f;

   m_i.re=0.0f;
   m_i.im=-1.0f;

   m_1.re=-1.0f;
   m_1.im=0.0f;

   if (mu==0)
   {
      _vector_mulc_f(r.c1,m_1,s.c3);
      _vector_mulc_f(r.c2,m_1,s.c4);
      _vector_mulc_f(r.c3,m_1,s.c1);
      _vector_mulc_f(r.c4,m_1,s.c2);
   }
   else if (mu==1)
   {
      _vector_mulc_f(r.c1,m_i,s.c4);
      _vector_mulc_f(r.c2,m_i,s.c3);
      _vector_mulc_f(r.c3,i,s.c2);
      _vector_mulc_f(r.c4,i,s.c1);
   }
   else if (mu==2)
   {
      _vector_mulc_f(r.c1,m_1,s.c4);
      r.c2=s.c3;
      r.c3=s.c2;
      _vector_mulc_f(r.c4,m_1,s.c1);
   }
   else if (mu==3)
   {
      _vector_mulc_f(r.c1,m_i,s.c3);
      _vector_mulc_f(r.c2,i,s.c4);
      _vector_mulc_f(r.c3,i,s.c1);
      _vector_mulc_f(r.c4,m_i,s.c2);
   }
   else
   {
      r.c1=s.c1;
      r.c2=s.c2;
      _vector_mulc_f(r.c3,m_1,s.c3);
      _vector_mulc_f(r.c4,m_1,s.c4);
   }

   return r;
}


int main(int argc,char *argv[])
{
   int i,j,n,ix,bc,mu,level,seed;
   int x0,x1,x2,x3;
   int np[4];
   float ran[4],sp[4];
   double pi,p[4];
   float *rs,r,mp,sig,px;
   complex z;
   suNf_spinor s,s0,s1,s2,s3,ps0[VOLUME],ps1[VOLUME],ps2[VOLUME];

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   printf("The lattice size is %dx%d^3\n",T,L);
   printf("\n");
   printf("\n");
   printf("Action of Qhat on plane waves\n");
   printf("-----------------------------\n\n");

   level=0;
   seed=123;
   rlxs_init(level,seed);
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
 
   geometry_eo_lexi();
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   printf("allocating gfield_f\n");
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();

   pi=4.0*atan(1.0);
   n=10;
   bc=1; /* 0=> periodic ; 1=> antiperiodic in time */
   printf("bc=%d\n",bc);
   
   for (i=0;i<n;i++)
   {
      ranlxs(ran,4);
      
      np[0]=(int)(ran[0]*(float)(T));
      np[1]=(int)(ran[1]*(float)(L));
      np[2]=(int)(ran[2]*(float)(L));
      np[3]=(int)(ran[3]*(float)(L));

      p[0]=((double)(np[0])*2.0+bc)*pi/(double)(T);
      p[1]=(double)(np[1])*2.0*pi/(double)(L);
      p[2]=(double)(np[2])*2.0*pi/(double)(L);
      p[3]=(double)(np[3])*2.0*pi/(double)(L);

      mp=(float)(hmass);
      mp+=(float)(1.0-cos(p[0]));
      mp+=(float)(1.0-cos(p[1]));
      mp+=(float)(1.0-cos(p[2]));
      mp+=(float)(1.0-cos(p[3]));      
      
      sp[0]=(float)(sin(p[0]));
      sp[1]=(float)(sin(p[1]));
      sp[2]=(float)(sin(p[2]));
      sp[3]=(float)(sin(p[3]));
      
      rs=(float*)(&s);
      r=0.0f;
      while ((1.0f+r)==1.0f)
      {
         gauss(rs,8*NF);
         r=0.0f;
         
         for (j=0;j<8*NF;j++)
            r+=rs[j]*rs[j];
         r=(float)(sqrt((double)(r)));
      }
      
      for (x0=0;x0<T;x0++)
      {
         for (x1=0;x1<L;x1++)
         {
            for (x2=0;x2<L;x2++)
            {
               for (x3=0;x3<L;x3++)
               {
                  ix=ipt[x0][x1][x2][x3];
                  
                  px=p[0]*(double)(x0)+p[1]*(double)(x1)
                     +p[2]*(double)(x2)+p[3]*(double)(x3);
                  
                  z.re=(float)(cos(px));
                  z.im=(float)(sin(px));
                  
                  _vector_mulc_f(s0.c1,z,s.c1);
                  _vector_mulc_f(s0.c2,z,s.c2);
                  _vector_mulc_f(s0.c3,z,s.c3);
                  _vector_mulc_f(s0.c4,z,s.c4);
                  
                  ps0[ix]=s0;
                  
                  z.re=mp;
                  z.im=0.0f;
                  
                  _vector_mulc_f(s1.c1,z,s0.c1);
                  _vector_mulc_f(s1.c2,z,s0.c2);
                  _vector_mulc_f(s1.c3,z,s0.c3);
                  _vector_mulc_f(s1.c4,z,s0.c4);
                  
                  for (mu=0;mu<4;mu++)
                  {
                     s2=mul_gamma(mu,s0);
                     
                     z.re=0.0f;
                     z.im=sp[mu];

                     _vector_mulc_f(s3.c1,z,s2.c1);
                     _vector_mulc_f(s3.c2,z,s2.c2);
                     _vector_mulc_f(s3.c3,z,s2.c3);
                     _vector_mulc_f(s3.c4,z,s2.c4);
                     
                     _vector_add_assign_f(s1.c1,s3.c1);
                     _vector_add_assign_f(s1.c2,s3.c2);
                     _vector_add_assign_f(s1.c3,s3.c3);
                     _vector_add_assign_f(s1.c4,s3.c4);
                  }
                  ps1[ix]=s1;
               }
            }
         }
      }
      
      Dphi(hmass,ps2,ps0);
      
      spinor_field_mul_add_assign_f(ps1,-1.0,ps2);
      sig=spinor_field_sqnorm_f(ps1)/spinor_field_sqnorm_f(ps0);
         
      printf("Maximal normalized difference = %.2e at p=(%d,%d,%d,%d), bc=%+d\n",sqrt(sig),
             np[0],np[1],np[2],np[3],(int)(1.0-2.0*bc));
      printf("(should be around 1*10^(-6) or so)\n\n");
   }
   
   free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif

   exit(0);
}
