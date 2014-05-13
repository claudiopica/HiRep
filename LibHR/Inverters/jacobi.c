/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Diagonalization of a square matrix using the Jacobi method
*
* See chapter 11 in
* 
*   W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery,
*   Numerical Recipes in FORTRAN, 2nd Edition 
*   (Cambridge University Press, Cambridge, 1992)
*
* for a description of the algorithm in the real case
*
* The externally accessible functions are
*
*   void jacobi1(int n,float a[],float d[],float v[])
*     Computes the eigenvalues and eigenvectors of a real symmetric matrix
*
*   void jacobi2(int n,complex a[],float d[],complex v[])
*     Computes the eigenvalues and eigenvectors of a complex hermitian matrix
*
* In both cases the matrix elements are assumed to be given by a[n*i+j],
* where the indices i and j range from 0 to n-1. On output the eigenvalues 
* and the associated orthonormal eigenvectors are d[j] and v[n*i+j] such
* that sum_k a[n*i+k]*v[n*k+j]=d[j]*v[n*i+j]. The eigenvalues are sorted
* in ascending order
*
* Revision information:
* $Id: jacobi.c,v 1.2 2005/07/20 10:10:13 luscher Exp $
*
* Author: Martin Luescher <luscher@mail.cern.ch>
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "error.h"
#include "complex.h"
#include "suN.h"
#include "inverters.h"

#define MAX_SWEEP 100


static void sort1(int n,double d[],double v[])
{
   int i,j,k;
   double p;
    
   for (i=0;i<n-1;i++)
   {
      k=i;
      p=d[i];
      
      for (j=i+1;j<n;j++)
      {
	 if (d[j]<p)
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
	    p=v[n*j+i];
	    v[n*j+i]=v[n*j+k];
	    v[n*j+k]=p;
	 }
      }
   }
}


static void sort2(int n,double d[],complex v[])
{
   int i,j,k;
   double p;
   complex q;
   
   for (i=0;i<n-1;i++)
   {
      k=i;
      p=d[i];
      
      for (j=i+1;j<n;j++)
      {
	 if (d[j]<p)
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


void jacobi1(int n,double a[],double d[],double v[])
{
   int k,l,j,sweep;
   double tol,abs_sum,thresh_factor,sd_factor,thresh;
   double r1,r2,r3,r4;
   double t,e,s,c,tau;
   double xn,xd0,xdh,xd1;
   
   xd0=0.0f;
   xdh=0.5f;
   xd1=1.0f;
   xn=(double)n;
   sd_factor=100.0f;
   thresh_factor=0.2f/(xn*xn);

   tol=xd0;
   abs_sum=xd0;
   
   for (k=0;k<n;k++)
   {
      v[n*k+k]=xd1;
      d[k]=a[n*k+k];
      tol+=fabs(d[k]);
      
      for (l=k+1;l<n;l++)
      {
	 v[n*k+l]=xd0;
	 v[n*l+k]=xd0;

	 error(a[n*k+l]!=a[n*l+k],1,
                   "jacobi1 [jacobi.c]","Matrix is not symmetric");

	 abs_sum+=fabs(a[n*k+l]);
      }
   }

   tol+=2.0f*abs_sum;
   tol*=DBL_EPSILON;
   
   for (sweep=0;(abs_sum>tol)&&(sweep<MAX_SWEEP);sweep++)
   {
      thresh=xd0;
      if (sweep<=2)
	 thresh=thresh_factor*abs_sum;
      
      for (k=0;k<n-1;k++)
      {
	 for (l=k+1;l<n;l++)
	 {
	    r1=sd_factor*(fabs(a[n*k+l]));
	    r2=fabs(d[k]);
	    r3=fabs(d[l]);

	    if ((sweep>3)&&(r1<=(r2*DBL_EPSILON))&&(r1<=(r3*DBL_EPSILON)))
	       a[n*k+l]=xd0;
	    
	    r1=fabs(a[n*k+l]);
	    if (r1<=thresh)
	       continue;

	    r2=d[l]-d[k];
	    r3=fabs(r2);

	    if ((sd_factor*r1)<(r3*DBL_EPSILON))
	    {
	       t=r1/r2;
	    }
	    else
	    {
	       r4=xdh*r2/r1;
	       if (r4<xd0)
	       {
		  t=xd1/(r4-(sqrt(xd1+r4*r4)));
	       }
	       else
	       {
		  t=xd1/(r4+(sqrt(xd1+r4*r4)));
	       }
	    }

	    e=xd1;
	    if (a[n*k+l]<xd0)
	       e=-xd1;
	    a[n*k+l]=xd0;

	    c=xd1/(sqrt(xd1+t*t));
	    s=t*c;
	    tau=s/(xd1+c);
	    
	    r2=t*r1;
	    d[k]-=r2;
	    d[l]+=r2;
	    
	    for (j=0;j<n;j++)
	    {
	       if (j<k)
	       {
		  r1=a[n*j+k];
		  r2=a[n*j+l];
		  a[n*j+k]=-s*( tau*r1+e*r2)+r1;
		  a[n*j+l]= s*(-tau*r2+e*r1)+r2;
	       }
	       else if ((j>k)&&(j<l))
	       {
		  r1=a[n*k+j];
		  r2=a[n*j+l];
		  a[n*k+j]=-s*( tau*r1+e*r2)+r1;
		  a[n*j+l]= s*(-tau*r2+e*r1)+r2;
	       }
	       else if (j>l)
	       {
		  r1=a[n*k+j];
		  r2=a[n*l+j];
		  a[n*k+j]=-s*( tau*r1+e*r2)+r1;
		  a[n*l+j]= s*(-tau*r2+e*r1)+r2;
	       }
		  
	       r1=v[n*j+k];
	       r2=v[n*j+l];
	       v[n*j+k]=-s*( tau*r1+e*r2)+r1;
	       v[n*j+l]= s*(-tau*r2+e*r1)+r2;
	    }
	 }
      }
      
      abs_sum=xd0;
      
      for (k=0;k<n-1;k++)
      {
	 for (l=k+1;l<n;l++)
	 {
	    abs_sum+=fabs(a[n*k+l]);
	 }
      }
   }

   error(sweep>=MAX_SWEEP,1,"jacobi1 [jacobi.c]",
             "Maximum number of sweeps exceeded");
   
   for (k=0;k<n-1;k++)
   {
      for (l=k+1;l<n;l++)
      {
	 a[n*k+l]=a[n*l+k];
      }
   }

   sort1(n,d,v);
}


void jacobi2(int n,complex a[],double d[],complex v[])
{
   int k,l,j,sweep;
   double tol,abs_sum,thresh_factor,sd_factor,thresh;
   double r1,r2,r3,r4;
   double t,s,c,tau;
   double xn,xd0,xdh,xd1;
   complex z1,z2;
   complex e;
   complex zd0,zd1;

   xd0=0.0f;
   xdh=0.5f;
   xd1=1.0f;
   xn=(double)n;
   sd_factor=100.0f;
   thresh_factor=0.2f/(xn*xn);

   zd0.re=xd0;
   zd0.im=xd0;
   zd1.re=xd1;
   zd1.im=xd0;

   tol=xd0;
   abs_sum=xd0;

   for (k=0;k<n;k++)
   {
      v[n*k+k]=zd1;
      d[k]=a[n*k+k].re;
      tol+=fabs(d[k]);

      for (l=k+1;l<n;l++)
      {
	 v[n*k+l]=zd0;
	 v[n*l+k]=zd0;

	 error((a[n*k+l].re!=a[n*l+k].re)||(a[n*k+l].im!=-a[n*l+k].im),1,
                   "jacobi2 [jacobi.c]","Matrix is not hermitian");
	  
         abs_sum+=(fabs(a[n*k+l].re)+fabs(a[n*k+l].im));
      }
   }

   tol+=2.0f*abs_sum;
   tol*=DBL_EPSILON;   
   
   for (sweep=0;(abs_sum>tol)&&(sweep<MAX_SWEEP);sweep++)
   {
      thresh=xd0;
      if (sweep<=2)
	 thresh=thresh_factor*abs_sum;

      for (k=0;k<n-1;k++)
      {
	 for (l=k+1;l<n;l++)
	 {
	    r1=sd_factor*
               (fabs(a[n*k+l].re)+fabs(a[n*k+l].im));
	    r2=fabs(d[k]);
	    r3=fabs(d[l]);

	    if ((sweep>3)&&(r1<=(r2*DBL_EPSILON))&&(r1<=(r3*DBL_EPSILON)))
	       a[n*k+l]=zd0;
	
	    r2=fabs(a[n*k+l].re);
	    r3=fabs(a[n*k+l].im);

	    if (r2>r3)
	    {
	       r3/=r2;
	       r1=r2*(sqrt(xd1+r3*r3));
	    }
	    else if (r2<r3)
	    {
	       r2/=r3;
	       r1=r3*(sqrt(xd1+r2*r2));
	    }
	    else
	    {
	       r1=r2*(sqrt(xd1+xd1));
	    }

	    if (r1<=thresh)
	       continue;

	    r2=d[l]-d[k];
	    r3=fabs(r2);

	    if ((sd_factor*r1)<(r3*DBL_EPSILON))
	    {
	       t=r1/r2;
	    }
	    else
	    {
	       r4=xdh*r2/r1;
	       if (r4<xd0)
	       {
		  t=xd1/(r4-(sqrt(xd1+r4*r4)));
	       }
	       else
	       {
		  t=xd1/(r4+(sqrt(xd1+r4*r4)));
	       }
	    }

	    e.re=a[n*k+l].re/r1;
	    e.im=a[n*k+l].im/r1;
	    a[n*k+l]=zd0;

	    c=xd1/(sqrt(xd1+t*t));
	    s=t*c;
	    tau=s/(xd1+c);

	    r2=t*r1;
	    d[k]-=r2;
	    d[l]+=r2;
	      
	    for (j=0;j<n;j++)
	    {
	       if (j<k)
	       {
		  z1=a[n*j+k];
		  z2=a[n*j+l];
		  a[n*j+k].re=-s*( tau*z1.re+e.re*z2.re+e.im*z2.im)+z1.re;
		  a[n*j+k].im=-s*( tau*z1.im-e.im*z2.re+e.re*z2.im)+z1.im;
		  a[n*j+l].re= s*(-tau*z2.re+e.re*z1.re-e.im*z1.im)+z2.re;
		  a[n*j+l].im= s*(-tau*z2.im+e.im*z1.re+e.re*z1.im)+z2.im;
	       }
	       else if ((j>k)&&(j<l))
	       {
		  z1=a[n*k+j];
		  z2=a[n*j+l];
		  a[n*k+j].re=-s*( tau*z1.re+e.re*z2.re+e.im*z2.im)+z1.re;
		  a[n*k+j].im=-s*( tau*z1.im+e.im*z2.re-e.re*z2.im)+z1.im;
		  a[n*j+l].re= s*(-tau*z2.re+e.re*z1.re+e.im*z1.im)+z2.re;
		  a[n*j+l].im= s*(-tau*z2.im+e.im*z1.re-e.re*z1.im)+z2.im;
	       }
	       else if (j>l)
	       {
		  z1=a[n*k+j];
		  z2=a[n*l+j];
		  a[n*k+j].re=-s*( tau*z1.re+e.re*z2.re-e.im*z2.im)+z1.re;
		  a[n*k+j].im=-s*( tau*z1.im+e.im*z2.re+e.re*z2.im)+z1.im;
		  a[n*l+j].re= s*(-tau*z2.re+e.re*z1.re+e.im*z1.im)+z2.re;
		  a[n*l+j].im= s*(-tau*z2.im-e.im*z1.re+e.re*z1.im)+z2.im;
	       }

	       z1=v[n*j+k];
	       z2=v[n*j+l];
	       v[n*j+k].re=-s*( tau*z1.re+e.re*z2.re+e.im*z2.im)+z1.re;
	       v[n*j+k].im=-s*( tau*z1.im-e.im*z2.re+e.re*z2.im)+z1.im;
	       v[n*j+l].re= s*(-tau*z2.re+e.re*z1.re-e.im*z1.im)+z2.re;
	       v[n*j+l].im= s*(-tau*z2.im+e.im*z1.re+e.re*z1.im)+z2.im;
	    }
	 }
      }

      abs_sum=xd0;
      
      for (k=0;k<n-1;k++)
      {
	 for (l=k+1;l<n;l++)
	 {
	    abs_sum+=
               (fabs(a[n*k+l].re)+fabs(a[n*k+l].im));
	 }
      }
   }

   error(sweep>=MAX_SWEEP,1,"jacobi2 [jacobi.c]",
             "Maximum number of sweeps exceeded");

   for (k=0;k<n-1;k++)
   {
      for (l=k+1;l<n;l++)
      {
	 a[n*k+l].re=a[n*l+k].re;
	 a[n*k+l].im=-a[n*l+k].im;
      }
   }

   sort2(n,d,v);  
}

