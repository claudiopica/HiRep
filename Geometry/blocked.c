/*******************************************************************************
*
* File geometry.c
*
* Definition of the lattice geometry
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "geometry.h"
#include "global.h"
#include "safe_mod.h"

static int init=0;
static int bl=1,nl=L;

static void set_block_size(void)
{
   if (L%4==0)
      bl=4;
   else if (L%3==0)
      bl=3;
   else if (L%2==0)
      bl=2;

   nl=L/bl;
   init=1;
}


static int index2(int x0,int x1,int x2,int x3)
{
   int y0,y1,y2,y3;
   int xb1,xb2,xb3;
   int xn1,xn2,xn3;
   int ib,in;
   
   if (init==0)
      set_block_size();
   
   y0=safe_mod(x0,T);
   y1=safe_mod(x1,L);
   y2=safe_mod(x2,L);
   y3=safe_mod(x3,L);

   xb1=y1%bl;
   xb2=y2%bl;
   xb3=y3%bl;

   xn1=y1/bl;
   xn2=y2/bl;
   xn3=y3/bl;   

   if ((y0%2)==0)
      ib=y0*bl*bl*bl+2*(xb3+xb2*bl+xb1*bl*bl);
   else
      ib=(y0-1)*bl*bl*bl+2*(xb3+xb2*bl+xb1*bl*bl)+1;

   in=xn3+xn2*nl+xn1*nl*nl;

   return(ib+in*T*bl*bl*bl);
}


void geometry_blocked(void)
{
   int x0,x1,x2,x3,ix,iy;

   for (x0=0;x0<T;x0++){
     iy=-1;
     for (x1=0;x1<L;x1++){
       for (x2=0;x2<L;x2++){
	 for (x3=0;x3<L;x3++){
	   iy++;
	   ix=index2(x0,x1,x2,x3);
	   ipt[x0][x1][x2][x3]=ix;
	   /*ipt_4d[x0][iy]=ix;*/
	   
	   iup[ix][0]=index2(x0+1,x1,x2,x3);
	   idn[ix][0]=index2(x0-1,x1,x2,x3);
	   iup[ix][1]=index2(x0,x1+1,x2,x3);
	   idn[ix][1]=index2(x0,x1-1,x2,x3);
	   iup[ix][2]=index2(x0,x1,x2+1,x3);
	   idn[ix][2]=index2(x0,x1,x2-1,x3);
	   iup[ix][3]=index2(x0,x1,x2,x3+1);
	   idn[ix][3]=index2(x0,x1,x2,x3-1);
	   /* tslice[ix]=x0; */
	 }
       }
     }
   }
}

static int index_noT(int x0,int x1,int x2,int x3)
{
   int y0,y1,y2,y3;
   int xb1,xb2,xb3;
   int xn1,xn2,xn3;
   int ib,in;
   
   if (init==0)
      set_block_size();
   
   y0=safe_mod(x0,T);
   y1=safe_mod(x1,L);
   y2=safe_mod(x2,L);
   y3=safe_mod(x3,L);

   xb1=y1%bl;
   xb2=y2%bl;
   xb3=y3%bl;

   xn1=y1/bl;
   xn2=y2/bl;
   xn3=y3/bl;   

   ib=y0*bl*bl*bl+(xb3+xb2*bl+xb1*bl*bl);

   in=xn3+xn2*nl+xn1*nl*nl;

   return(ib+in*T*bl*bl*bl);
}

void geometry_blocked_noT(void)
{
   int x0,x1,x2,x3,ix;

   for (x0=0;x0<T;x0++){
     for (x1=0;x1<L;x1++){
       for (x2=0;x2<L;x2++){
	 for (x3=0;x3<L;x3++){
	   ix=index_noT(x0,x1,x2,x3);
	   ipt[x0][x1][x2][x3]=ix;
	   
	   iup[ix][0]=index_noT(x0+1,x1,x2,x3);
	   idn[ix][0]=index_noT(x0-1,x1,x2,x3);
	   iup[ix][1]=index_noT(x0,x1+1,x2,x3);
	   idn[ix][1]=index_noT(x0,x1-1,x2,x3);
	   iup[ix][2]=index_noT(x0,x1,x2+1,x3);
	   idn[ix][2]=index_noT(x0,x1,x2-1,x3);
	   iup[ix][3]=index_noT(x0,x1,x2,x3+1);
	   idn[ix][3]=index_noT(x0,x1,x2,x3-1);
	   /* tslice[ix]=x0; */
	 }
       }
     }
   }
}

