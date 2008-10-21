/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File geometry.c
*
* Definition of the lattice geometry
*
*******************************************************************************/

#include "geometry.h"
#include "global.h"
#include "safe_mod.h"

static int init=0;
static int xbl,ybl,zbl,xnl,ynl,znl;

static int block_size(int L)
{
	int bl=1;
   if (L%4==0)
      bl=4;
   else if (L%3==0)
      bl=3;
   else if (L%2==0)
      bl=2;

	 return bl;

}

static void init_block_size() 
{
	xbl=block_size(X);
	xnl=X/xbl;
	ybl=block_size(Y);
	ynl=Y/ybl;
	zbl=block_size(Z);
	znl=Z/zbl;
}


static int index2(int x0,int x1,int x2,int x3)
{
   int y0,y1,y2,y3;
   int xb1,xb2,xb3;
   int xn1,xn2,xn3;
   int ib,in;
   
   if (init==0) {
      init_block_size();
			init=1;
	 }
   
   y0=safe_mod(x0,T);
   y1=safe_mod(x1,X);
   y2=safe_mod(x2,Y);
   y3=safe_mod(x3,Z);

   xb1=y1%xbl;
   xb2=y2%ybl;
   xb3=y3%zbl;

   xn1=y1/xbl;
   xn2=y2/ybl;
   xn3=y3/zbl;   

   if ((y0%2)==0)
      ib=y0*xbl*ybl*zbl+2*(xb3+xb2*zbl+xb1*zbl*ybl);
   else
      ib=(y0-1)*xbl*ybl*zbl+2*(xb3+xb2*zbl+xb1*zbl*ybl)+1;

   in=xn3+xn2*znl+xn1*znl*ynl;

   return(ib+in*T*xbl*ybl*zbl);
}


void geometry_blocked(void)
{
   int x0,x1,x2,x3,ix,iy;

	 geometry_init();

   for (x0=0;x0<T;x0++){
     iy=0;
     for (x1=0;x1<X;x1++){
       for (x2=0;x2<Y;x2++){
				 for (x3=0;x3<Z;x3++){
					 ix=index2(x0,x1,x2,x3);
					 ipt(x0,x1,x2,x3)=ix;
					 ipt_4d(x0,iy)=ix;

					 iup(ix,0)=index2(x0+1,x1,x2,x3);
					 idn(ix,0)=index2(x0-1,x1,x2,x3);
					 iup(ix,1)=index2(x0,x1+1,x2,x3);
					 idn(ix,1)=index2(x0,x1-1,x2,x3);
					 iup(ix,2)=index2(x0,x1,x2+1,x3);
					 idn(ix,2)=index2(x0,x1,x2-1,x3);
					 iup(ix,3)=index2(x0,x1,x2,x3+1);
					 idn(ix,3)=index2(x0,x1,x2,x3-1);

					 ++iy;
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
   
	 geometry_init();

   if (init==0) {
      init_block_size();
			init=1;
	 }
   
   y0=safe_mod(x0,T);
   y1=safe_mod(x1,X);
   y2=safe_mod(x2,Y);
   y3=safe_mod(x3,Z);

   xb1=y1%xbl;
   xb2=y2%ybl;
   xb3=y3%zbl;

   xn1=y1/xbl;
   xn2=y2/ybl;
   xn3=y3/zbl;   

   ib=y0*xbl*ybl*zbl+(xb3+xb2*zbl+xb1*ybl*zbl);

   in=xn3+xn2*znl+xn1*ynl*znl;

   return(ib+in*T*xbl*ybl*zbl);
}

void geometry_blocked_noT(void)
{
   int x0,x1,x2,x3,ix,iy;

   for (x0=0;x0<T;x0++){
		 iy=0;
     for (x1=0;x1<X;x1++){
			 for (x2=0;x2<Y;x2++){
				 for (x3=0;x3<Z;x3++){
					 ix=index_noT(x0,x1,x2,x3);
					 ipt(x0,x1,x2,x3)=ix;
					 ipt_4d(x0,iy)=ix;

					 iup(ix,0)=index_noT(x0+1,x1,x2,x3);
					 idn(ix,0)=index_noT(x0-1,x1,x2,x3);
					 iup(ix,1)=index_noT(x0,x1+1,x2,x3);
					 idn(ix,1)=index_noT(x0,x1-1,x2,x3);
					 iup(ix,2)=index_noT(x0,x1,x2+1,x3);
					 idn(ix,2)=index_noT(x0,x1,x2-1,x3);
					 iup(ix,3)=index_noT(x0,x1,x2,x3+1);
					 idn(ix,3)=index_noT(x0,x1,x2,x3-1);

					 ++iy;
				 }
			 }
     }
   }
}

