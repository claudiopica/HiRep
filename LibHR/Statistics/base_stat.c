/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File base_stat.c
*
* Basic function for statistical analisys of data series
*
*******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "statistics.h"
#include "error.h"


double average(int n,double a[])
  /*
   *     Returns the average of the data stored in the array a[n]
   */
{
   int i;
   double abar;

   abar=0.0;
   
   for (i=0;i<n;i++)
      abar+=a[i];

   abar/=(double)(n);
   
   return(abar);   
}


double sigma0(int n,double a[])
  /*
   *     Returns the standard deviation of the data in the array a[n]
   *     from their mean value divided by sqrt(n)
   */
{
   int i;
   double abar,var,fact;

   abar=0.0;
   var=0.0;
   
   for (i=0;i<n;i++)
   {
      abar+=a[i];
      var+=(a[i]*a[i]);
   }

   fact=1.0/(double)(n);
   abar*=fact;
   var=fact*var-abar*abar;

   return(sqrt(fact*fabs(var)));
}


void auto_corr(int n,double a[],int tmax,double gamma[])
  /*
   *     Computes the auto-correlation function gamma[tmax] of the data
   *     in the array a[n]
   */
{
   int i,t;
   double fact,sh,sl,s;

   error((n<1)||(tmax<1)||(n<tmax),1,"auto_corr [stat.c]",
         "Arguments out of range");
   sl=0.0;

   for (i=0;i<n;i++)
      sl+=a[i];

   sh=sl;
      
   for (t=0;t<tmax;t++)
   {
      fact=1.0/(double)(n-t);
      s=-fact*sl*sh;

      for (i=t;i<n;i++)
         s+=a[i-t]*a[i];
      
      gamma[t]=fact*s;
      
      sl-=a[n-t-1];
      sh-=a[t];      
   }
}


double sigma(int n,double a[],double *tau,int *flag)
  /*
   *     Returns the statistical error associated with the data series a[n] 
   *     taking auto-correlations into account. The calculated integrated
   *     auto-correlation time is assigned to the parameter tau. On exit
   *     flag=0 if the error estimation was stable and flag=1 otherwise 
   */
{
   int tmax,i,j,itest;
   double abar,sig0;
   double *g,*t,dt,del,var,taumax;

   tmax=n/30+1;
   g=malloc(tmax*sizeof(double));
   t=malloc(tmax*sizeof(double));   

   auto_corr(n,a,tmax,g);

   abar=average(n,a);
   sig0=sigma0(n,a);
   
   if (((fabs(abar)+sig0)==fabs(abar))||(g[0]==0.0))
   {
      *tau=0.5;
      *flag=0;
     free(g); free(t);
      return(sig0);
   }
   
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
      
      if (i>=(int)(5.0*t[i]))
      {
         itest=0;

         for (j=i+1;j<tmax;j++)
         {
            if (j>(i+(int)(3.0*t[i])))
               break;
            
            dt=t[j]-t[i];
            del=(double)(2*(2*j+1))/(double)(n);

            if ((dt*dt)>(del*t[i]*t[i]))
               itest=1;
         }

         if (itest==0)
         {
            var=2.0*t[i]*g[0]/(double)(n);
            *tau=t[i];
            *flag=0;
           
           free(g); free(t);
            return(sqrt(var));
         }
      }
   }

   var=2.0*taumax*g[0]/(double)(n);
   *tau=taumax;
   *flag=1;

  free(g); free(t);
   return(sqrt(var));
}

