/*******************************************************************************
*
* Random number generator "ranlxs"
*
* See the notes 
*
*   "User's guide for ranlxs and ranlxd [C programs]" (December 1997)
*
*   "Double precision implementation of the random number 
*    generator ranlux" (December 1997)
*
* for a detailed description
*
* The externally accessible functions are 
*
*   void ranlxs(float r[],int n)
*     Computes the next n single-precision random numbers and 
*     assigns them to the elements r[0],...,r[n-1] of the array r[]
* 
*   void rlxs_init(int level,int seed)
*     Initialization of the generator
*
*   void rlxs_get(int state[])
*     Extracts the current state of the generator and stores the 
*     information in the array state[25]
*
*   void rlxs_reset(int state[])
*     Resets the generator to the state defined by the array state[25]
*
* Version: 2.2
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 15.07.1999
*
*******************************************************************************/

#include <climits>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "ranlxs.h"

static int pr,ir,jr,is,is_old,init=0;
static int next[12],snext[24]; 
static float xflt[24];
static double zero,one,carry;
static double sbase,sone_bit,base,one_bit,shift;
static double xdbl[12],ydbl[12];


namespace //Local linkage functions
{

  void error(int no)
  {
    switch(no)
      {
      case 0:
	std::cerr<<"Error in rlxs_init\n";
	std::cerr<<"Arithmetic on this machine is not suitable for ranlxs\n";
	break;
      case 1:
	std::cerr<<"Error in subroutine rlxs_init\n";
	std::cerr<<"Bad choice of luxury level (should be 0,1 or 2)\n";
	break;
      case 2:
	std::cerr<<"Error in subroutine rlxs_init\n";
	std::cerr<<"Bad choice of seed (should be between 1 and 2^31-1)\n";
	break;
      case 3:
	std::cerr<<"Error in rlxs_get\n";
	std::cerr<<"Undefined state\n";
	break;
      case 4:
	std::cerr<<"Error in rlxs_reset\n";
	std::cerr<<"Arithmetic on this machine is not suitable for ranlxs\n";
	break;
      case 5:
	std::cerr<<"Error in rlxs_reset\n";
	std::cerr<<"Unexpected input data\n";
	break;
      }         
    std::cerr<<"Program aborted\n";
    exit(0);
  }
  
  extern inline void ranlux_step(double &x1, double &x2, const int i1, const int i2, const int i3)
  {
    x1=xdbl[i1]-xdbl[i2];          
    if(x2<zero){                              
      x1-=one_bit;                 
      x2+=one;                     
    }                              
    xdbl[i3]=x2;
  }


  void update()
  {
    int k,kmax,l;
    double x,y1,y2,y3;

    for (k=0;ir>0;++k) 
      {
	y1=xdbl[jr]-xdbl[ir];
	y2=y1-carry;
	if (y2<zero) 
	  { 
	    carry=one_bit;
	    y2+=one;
	  }
	else 
	  carry=zero;
	xdbl[ir]=y2;
	ir=next[ir];
	jr=next[jr];
      }

    kmax=pr-12;

    for (;k<=kmax;k+=12)
      {
	y1=xdbl[7]-xdbl[0];
	y1-=carry;

	ranlux_step(y2,y1, 8, 1, 0);
	ranlux_step(y3,y2, 9, 2, 1);
	ranlux_step(y1,y3,10, 3, 2);
	ranlux_step(y2,y1,11, 4, 3);
	ranlux_step(y3,y2, 0, 5, 4);
	ranlux_step(y1,y3, 1, 6, 5);
	ranlux_step(y2,y1, 2, 7, 6);
	ranlux_step(y3,y2, 3, 8, 7);
	ranlux_step(y1,y3, 4, 9, 8);
	ranlux_step(y2,y1, 5,10, 9);
	ranlux_step(y3,y2, 6,11,10);
      
	if (y3<zero)
	  {
	    carry=one_bit;
	    y3+=one;
	  }
	else
	  carry=zero;
	xdbl[11]=y3;
      }  

    kmax=pr;

    for (;k<kmax;++k) 
      {
	y1=xdbl[jr]-xdbl[ir];
	y2=y1-carry;
	if (y2<zero) 
	  { 
	    carry=one_bit;
	    y2+=one;
	  }
	else 
	  carry=zero;
	xdbl[ir]=y2;
	ydbl[ir]=y2+shift;
	ir=next[ir];
	jr=next[jr];
      }

    ydbl[ir]=xdbl[ir]+shift;

    for (k=next[ir];k>0;)
      {
	ydbl[k]=xdbl[k]+shift;
	k=next[k];
      }

    for (k=0,l=0;k<12;++k) 
      {
	x=xdbl[k];
	y2=ydbl[k]-shift;
	if (y2>x)
	  y2-=sone_bit;
	y1=(x-y2)*sbase;

	xflt[l++]=(float)y1;
	xflt[l++]=(float)y2;
      }

    is=ir+ir;
    is_old=is;
  }


  void define_constants()
  {
    int k;

    init=1;
    zero=0.0;
    one=1.0;
    sbase=ldexp(one,24);
    sone_bit=ldexp(one,-24);
    base=ldexp(one,48);
    one_bit=ldexp(one,-48);
    shift=ldexp(one,DBL_MANT_DIG-25);

    for (k=0;k<12;++k) 
      {
	next[k]=(k+1)%12;
	snext[2*k]=(2*k+1)%24;
	snext[2*k+1]=(2*k+2)%24;
      }
  }

}

void rlxs_init(int level,int seed)
{
   int ibit,jbit,i,k,l,xbit[31];
   double x,y;

   if ((INT_MAX<2147483647)||
       (FLT_RADIX!=2)||(FLT_MANT_DIG<24)||(DBL_MANT_DIG<48))
      error(0);

   if      (level==0)
      pr=109;
   else if (level==1)
      pr=202;
   else if (level==2)
      pr=397;
   else
      error(1);

   define_constants();
   i=seed;

   for (k=0;k<31;++k) 
   {
      xbit[k]=i%2;
      i/=2;
   }

   if ((seed<=0)||(i!=0))
      error(2);

   ibit=0;
   jbit=18;
 
   for (k=0;k<12;++k) 
   {
      x=zero;

      for (l=1;l<=48;++l) 
      {
         y=(double)xbit[ibit];
         x+=x+y;
         xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
         ibit=(ibit+1)%31;
         jbit=(jbit+1)%31;
      }
      xdbl[k]=one_bit*x;
   }

   carry=zero;
   ir=0;
   jr=7;
   is=23;
   is_old=0;
}

#include <ctime>

void ranlxs(float r[],int n)
{
   int k;

   if (init==0) {
      rlxs_init(2,time(0));
      //rlxs_init(0,1);
   }

   for (k=0;k<n;++k) 
   {
      is=snext[is];
      if (is==is_old)
         update();
      r[k]=xflt[is];
   }
}

void rlxs_get(int state[])
{
   int k;
   double x,y1,y2;

   if (init==0)
      error(3);

   for (k=0;k<12;++k) 
   {
      x=sbase*xdbl[k];
      y1=sbase*modf(x,&y2);
      state[2*k]=(int)y1;
      state[2*k+1]=(int)y2;
   }

   k=12*pr+ir;
   k=12*k+jr;
   k=24*k+is;
   state[24]=2*k+(int)(carry*base);
}


void rlxs_reset(int state[])
{
   int k;
   double y1,y2;

   if ((INT_MAX<2147483647)||
       (FLT_RADIX!=2)||(FLT_MANT_DIG<24)||(DBL_MANT_DIG<48))
      error(4);

   define_constants();

   for (k=0;k<24;++k) 
   {
      if ((state[k]>=(int)sbase)||(state[k]<0))
         error(5);
   }

   k=state[24];
   if (k<0)
      error(5);
   carry=one_bit*(double)(k%2);
   k/=2;
   is=k%24;
   k/=24;
   jr=k%12;
   k/=12;
   ir=k%12;
   pr=k/12;
   is_old=2*ir;

   if (((pr!=109)&&(pr!=202)&&(pr!=397))||
       (jr!=((ir+7)%12)))
      error(5);

   for (k=0;k<12;++k) 
   {
      y1=(double)state[2*k];
      y2=(double)state[2*k+1];
      xdbl[k]=one_bit*(y1+y2*sbase);
      xflt[2*k]=(float)(sone_bit*y1);
      xflt[2*k+1]=(float)(sone_bit*y2);
   }
}
