/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File random_su2.c
*
* Random SU(2) matrices functions
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "update.h"
#include "hr_omp.h"

#define NVEC (32)
#define NRAN (2*NVEC)
#define PI_HALF 1.5707964f
#define PI 3.1415927f
#define TWO_PI 6.2831854f

static int i_vec=NVEC,i_y=NRAN,i_v=NRAN;
static double vec1[NVEC],vec2[NVEC],vec3[NVEC];
static double r[NRAN],u[NRAN],v[NRAN],y[NRAN];


static void update_vec(void)
{
   int i;
   double r1,r2,rsq;
     
   ranlxd(r,NRAN);
      
   for (i=0;i<NVEC;i++)
   {
      r1=2.0*r[i]-1.0;
      r2=TWO_PI*r[NVEC+i]-PI;      
      rsq=sqrt(1.0-r1*r1);

      vec1[i]=r1;
      vec2[i]=rsq*sin(r2);
      vec3[i]=rsq*cos(r2);      
   }

   i_vec=0;
}


static void update_y(void)
{
   int i;
   double r1,r2,r3,r4,s,c;

   ranlxd(y,NRAN);
   ranlxd(u,NRAN);
   ranlxd(r,NRAN);
      
   for (i=0;i<NVEC;i++)
   {
      r1=-log(1.0-y[i]);
      r2=PI_HALF*y[NVEC+i];
      r3=log(1.0-u[i]);
      r4=log(1.0-u[NVEC+i]);

      s=sin(r2);
      s*=s;
      c=1.0-s;

      y[i]=r1*s-r3;
      y[NVEC+i]=r1*c-r4;

      r1=r[i]*r[i];
      r2=r[NVEC+i]*r[NVEC+i];
      u[i]=r1+r1;
      u[NVEC+i]=r2+r2;
   }
   
   i_y=0;   
}

void random_su2(double rho,double s[])
  /*
   *  Computes a random vector s[4] with probability density
   *  proportional to exp(rho*s[0])*delta(1-s^2) assuming rho>=0
   */
{
   double rhoinv,s0p1,ut,rt;
   double s0,s1,s2,s3,sq;

_OMP_PRAGMA( critical )
{
   if (i_vec==NVEC)
      update_vec();
   
   if (rho>1.5)
   {
      rhoinv=1.0/rho;

      for (;;)
      {         
         if (i_y==NRAN)
            update_y();

         s0p1=2.0-rhoinv*y[i_y];
         ut=u[i_y++];

         if (ut<=s0p1)
            break;
      }
   }
   else if (rho>0.3)
   {
      rhoinv=1.0/rho;
      rt=exp(rho+rho)-1.0;

      for (;;)
      {
         if (i_v==NRAN)
         {
            ranlxd(v,NRAN);
            i_v=0;
         }         

         s0p1=rhoinv*log(1.0+rt*v[i_v++]);
         ut=v[i_v++];
         
         if ((ut*ut)<=(s0p1*(2.0-s0p1)))
            break;
      }
   }
   else
   {
      for (;;)
      {
         if (i_v==NRAN)
         {
            ranlxd(v,NRAN);
            i_v=0;
         }

         s0p1=2.0*v[i_v++];
         rt=exp(rho*(s0p1-2.0));
         ut=v[i_v++];
         
         if ((ut*ut)<=(s0p1*(2.0-s0p1)*rt*rt))
            break;
      }
   }

   sq=sqrt(s0p1*(2.0-s0p1));
   s0=s0p1-1.0;
   s1=sq*vec1[i_vec];
   s2=sq*vec2[i_vec];
   s3=sq*vec3[i_vec];
   
   sq=s0*s0+s1*s1+s2*s2+s3*s3;
   sq=1.5-0.5*sq;

   s[0]=sq*s0;
   s[1]=sq*s1;
   s[2]=sq*s2;
   s[3]=sq*s3;

   i_vec+=1;
} //OMP Critical section
  
}
