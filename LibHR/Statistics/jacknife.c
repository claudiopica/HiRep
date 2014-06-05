/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File jacknife.c
*
* Functions for jacknife error estimate
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "statistics.h"


double auto_corr_time(int n,int tmax,double g[],int *flag)
  /*
   *     Returns the auto-corelation time associated with auto-correlation
   *     function g[tmax] formed from n measurements
   *     On exit flag=0 if the estimation of the auto-correlation time was
   *     stable  and flag=1 otherwise
   */
{
   int i,j,itest;
   double par1,par2;   
   double *t,dt,del,taumax;
      
   par1=5.0;
   par2=3.0;
         
   t=malloc(tmax*sizeof(*t));
         
   t[0]=0.5;  
         
   for (i=1;i<tmax;i++)
   {
      t[i]=t[i-1]+g[i]/g[0];

      if (t[i]<=0.0)
         tmax=i;
   }

   taumax=0.0;
      
   for (i=0;i<tmax;i++)
   {
      if (t[i]>taumax)
         taumax=t[i];
         
      if (i>=(int)(par1*t[i]))
      {  
         itest=0;
    
         for (j=i+1;j<tmax;j++)
         {
            if (j>(i+(int)(par2*t[i])))
               break;
    
            dt=t[j]-t[i];
            del=(double)(2*(2*j+1))/(double)(n);
      
            if ((dt*dt)>(del*t[i]*t[i]))
               itest=1;
         }
         
         if (itest==0)
         {
            *flag=0;
           dt=t[i];
           free(t);
            return(dt);
         }
      }
   }
            
   *flag=1;
  free(t);
  
   return taumax;
}


double sigma_bin(int n,int binsize,double a[])
  /*
   *     Forms from the data in the array a[n], a new set of data b[numbin]
   *     consisting of bin averages of a[] of bin-length binsize. 
   *     Returns the standard deviation of b[] from its mean value 
   *     divided by sqrt(numbin)
   */
{
   int i,j,icount,numbin;
   double *b,s0;
      
   numbin=n/binsize;
   b=malloc(numbin*sizeof(*b));
   
   icount=0;
   for (i=0;i<numbin;i++)   
   {
      b[i]=0.0;
      for (j=0;j<binsize;j++)
      {
         b[i]+=a[icount];
         icount++;
      }
      b[i]/=(double)(binsize);
   }
 
   s0=sigma0(numbin,b);
  free(b);
  
  return s0;
}


double sigma_replicas(int n,int r,double a[],double *tau,int *flag)
  /*
   *     Returns the statistical error associated with the data series a[n]
   *     containing the measurements from r replicas.
   *     It is assumed that the series a[] is arranged 
   *     starting with all data from replica 1,
   *     followed by all those from replica 2, etc.
   *     The auto-correlation functions from each replica are computed
   *     separately and then an average is made over them.
   *     The calculated integrated auto-correlation time from this 
   *     averaged auto-correlation function is assigned to the parameter tau
   *     On exit flag=0 if the error estimation was stable and flag=1 otherwise
   */
{
   int tmax,i,j,icount,ipar,nr;
   double *ar,*gr,*g,var,abar,sig0;

   ipar=30;

   tmax=n/ipar+1;

   nr=n*r;

   g=malloc(tmax*sizeof(double));
   for (i=0;i<tmax;i++)
      g[i]=0.0;   

   gr=malloc(tmax*sizeof(double));
   ar=malloc(n*sizeof(double)); 

   icount=0;
   for (j=0;j<r;j++)
   {
      for (i=0;i<n;i++)
      {
         ar[i]=a[icount];
         icount++;
      }   

      auto_corr(n,ar,tmax,gr);
 
      for (i=0;i<tmax;i++)
         g[i]+=gr[i];
   } 
          
   for (i=0;i<tmax;i++)
      g[i]/=(double)(r);

   abar=average(nr,a);
   sig0=sigma0(nr,a);
   
   if (((fabs(abar)+sig0)==fabs(abar))||(g[0]==0.0))
   {
      *tau=0.5;
      *flag=0;
     free(g);
     free(gr);
     free(ar);

      return(sig0);
   }

   *tau=auto_corr_time(nr,tmax,g,flag);

   var=2.0*(*tau)*g[0]/(double)(nr);

  free(g);
  free(gr);
  free(ar);
  
   return(sqrt(var));
}


double sigma_jackknife(int nobs,int n,double a[],double *ave_j,
                       double (*pobs)(double v[]))
  /*
   *     Here a[] stores n measurements of nobs different observables
   *     a[0],...,a[n-1]  contains the measurements of obervable1
   *     a[n],...,a[2n-1] contains the measurements of obervable2
   *     etc
   *     *pobs is a pointer to a function F of the average of the obervables 
   *     *ave_j gives the best estimate of F, 
   *     and sigma_jackknife returns the estimated jackknife error.
   */
{
   int i,j,k;
   double fact,*f,*atot,*aj;
   
   atot=malloc(nobs*sizeof(double));
   aj=malloc(nobs*sizeof(double));
   f=malloc(n*sizeof(double));

   fact=1.0/(double)(n);
   
   k=0;
   for (i=0;i<nobs;i++)
   {
      atot[i]=0.0;

      for (j=0;j<n;j++)
      {
         atot[i]+=a[k];
         k++; 
      } 

      aj[i]=atot[i]*fact;
   }
   
   *ave_j=pobs(aj); 
      
   fact=1.0/(double)(n-1);

   for (j=0;j<n;j++)   
   {  
      k=j;
      for (i=0;i<nobs;i++)
      {   
         aj[i]=(atot[i]-a[k])*fact;
         k+=n;
      }

      f[j]=pobs(aj);
   } 
   
   fact=sigma0(n,f)*(double)(n-1);
  
  free(atot);
  free(aj);
  free(f);
  
  return fact;
}

