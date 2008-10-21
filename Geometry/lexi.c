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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "geometry.h"
#include "global.h"
#include "safe_mod.h"

static int index_lexi(int x0,int x1,int x2,int x3)
{
   int y0,y1,y2,y3;

   y0=safe_mod(x0,T);
   y1=safe_mod(x1,X);
   y2=safe_mod(x2,Y);
   y3=safe_mod(x3,Z);

   return(y3+Z*(y2+Y*(y1+X*y0)));
}

void geometry_lexi(void)
{
   int x0,x1,x2,x3,ix,iy;

	 geometry_init();

   for (x0=0;x0<T;x0++){
		 iy=0;
     for (x1=0;x1<X;x1++){
       for (x2=0;x2<Y;x2++){
				 for (x3=0;x3<Z;x3++){
					 ix=index_lexi(x0,x1,x2,x3);
					 ipt(x0,x1,x2,x3)=ix;
					 ipt_4d(x0,iy)=ix;

					 iup(ix,0)=index_lexi(x0+1,x1,x2,x3);
					 idn(ix,0)=index_lexi(x0-1,x1,x2,x3);
					 iup(ix,1)=index_lexi(x0,x1+1,x2,x3);
					 idn(ix,1)=index_lexi(x0,x1-1,x2,x3);
					 iup(ix,2)=index_lexi(x0,x1,x2+1,x3);
					 idn(ix,2)=index_lexi(x0,x1,x2-1,x3);
					 iup(ix,3)=index_lexi(x0,x1,x2,x3+1);
					 idn(ix,3)=index_lexi(x0,x1,x2,x3-1);

					 ++iy;
				 }
			 }
     }
   }
}

