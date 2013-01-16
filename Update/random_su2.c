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

// Maybe the following are already there

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"

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
      r1=2.0f*r[i]-1.0f;
      r2=TWO_PI*r[NVEC+i]-PI;      
      rsq=sqrt(1.0f-r1*r1);

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
      r1=-log(1.0f-y[i]);
      r2=PI_HALF*y[NVEC+i];
      r3=log(1.0f-u[i]);      
      r4=log(1.0f-u[NVEC+i]);

      s=sin(r2);
      s*=s;
      c=1.0f-s;

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

   if (i_vec==NVEC)
      update_vec();
   
   if (rho>1.5f)
   {
      rhoinv=1.0f/rho;

      for (;;)
      {         
         if (i_y==NRAN)
            update_y();

         s0p1=2.0f-rhoinv*y[i_y];
         ut=u[i_y++];

         if (ut<=s0p1)
            break;
      }
   }
   else if (rho>0.3f)
   {
      rhoinv=1.0f/rho;
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
         
         if ((ut*ut)<=(s0p1*(2.0f-s0p1)))
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

         s0p1=2.0f*v[i_v++];
         rt=exp(rho*(s0p1-2.0f));
         ut=v[i_v++];
         
         if ((ut*ut)<=(s0p1*(2.0f-s0p1)*rt*rt))
            break;
      }
   }

   sq=sqrt(s0p1*(2.0f-s0p1));   
   s0=s0p1-1.0f;   
   s1=sq*vec1[i_vec];
   s2=sq*vec2[i_vec];
   s3=sq*vec3[i_vec];
   
   sq=s0*s0+s1*s1+s2*s2+s3*s3;
   sq=1.5f-0.5f*sq;

   s[0]=sq*s0;
   s[1]=sq*s1;
   s[2]=sq*s2;
   s[3]=sq*s3;

   i_vec+=1;
}


double random_so2(double a)
//
//  Computes a random angle theta \in [-\pi,\pi] with probability density
//  proportional to exp(a*\cos(theta)) assuming rho>=0
//
{

double x,xprime,epsi,astar,mas,alpha,beta,delta,betafrac,hx,gx;
 
 //  Method used here:
 // arXiv:hep-lat/9210016v1
 // "IMPROVEMENT OF EFFICIENCY IN GENERATING RANDOM U(1) 
 //                    VARIABLES WITH BOLTZMANN DISTRIBUTION"
 // By T.Hattori and H.Nakajima
 /*
 %\cite{Hattori:1992qk}
\bibitem{Hattori:1992qk} 
  T.~Hattori and H.~Nakajima,
  %``Improvement of efficiency in generating random U(1) variables with Boltzmann distribution,''
  Nucl.\ Phys.\ Proc.\ Suppl.\  {\bf 26}, 635 (1992)
  [hep-lat/9210016].
  %%CITATION = HEP-LAT/9210016;%%
  */      
      
      for (;;)
      {
         if (i_v==NRAN)	// Making sure we have enough random numbers
         {
            ranlxd(v,NRAN);
            i_v=0;
         }

		epsi	= 0.001;
		astar	= 0.798953686083986;
		mas		= max(0,a-astar);
		delta=0.35*mas+1.03*sqrt(mas);
		
		// Step 1
		alpha=min(sqrt(a*(2-epsi)),max(sqrt(epsi*a),delta));
		beta=max(alpha*alpha/a,(cosh(PI*alpha)-1)/(exp(2*a)-1))-1;
		
		// Step 2
		betafrac=sqrt((1+beta)/(1-beta));
		x=v[i_v++];
        xprime=v[i_v++];				// Random number in [0 ; 1[
        hx=(2/alpha)*atanh(   betafrac*tan(  (2*x-1)*atan(  tanh(PI*alpha/2)/betafrac  )  )  );
		gx=exp(-a*(1-cos(hx)))*(cosh(alpha*hx)+beta)/(1+beta);
         
         if (xprime <= gx)
         {
         	return hx;
            break;
    	 }
      }
   

}

double random_so2_2(double a,double b)
//
//  Computes a random angle theta \in [-\pi,\pi] with probability density
//  proportional to exp(a*\cos(theta)) assuming rho>=0
//
{
     
      double maxV=exp(sqrt(a*a+b*b));
         
      
      for (;;)
      {
         if (i_v==NRAN)	// Making sure we have enough random numbers
         {
            ranlxd(v,NRAN);
            i_v=0;
         }
         
         double y=maxV*v[i_v++];
         double x=2*PI*v[i_v++]-PI;
         if (y<exp(a*cos(x)+b*sin(x))) return x;
      }
}         

