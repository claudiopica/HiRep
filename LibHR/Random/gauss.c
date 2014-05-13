/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File gauss.c
* 
* Generation of Gaussian random numbers
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327 
#endif

void gauss(double r[],int n)
{
   int k;
   double u[2];
   double x1,x2,rho,y1,y2;

   for (k=0;k<n;)
   {
      ranlxd(u,2);
      x1=u[0];
      x2=u[1];

      rho=-log(1.0-x1);
      rho=sqrt(2.*rho);
      x2*=2.0*M_PI;
      y1=rho*sin(x2);
      y2=rho*cos(x2);
      
      r[k++]=y1;
      if (k<n)
         r[k++]=y2;
   }
}


void gauss_flt(float rd[],int n)
{
   int k;
   double ud[2];
   double x1,x2,rho,y1,y2;

   for (k=0;k<n;)
   {
      ranlxd(ud,2);
      x1=ud[0];
      x2=ud[1];

      rho=-log(1.0-x1);
      rho=sqrt(2.*rho);
      x2*=2.0*M_PI;
      y1=rho*sin(x2);
      y2=rho*cos(x2);
      
      rd[k++]=(float)y1;
      if (k<n)
         rd[k++]=(float)y2;
   }
}

