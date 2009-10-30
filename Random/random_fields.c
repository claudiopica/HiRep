/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File random_fields.c
*
* Pseudorandom generation of fields
*
*******************************************************************************/

#include <stdlib.h>

#include "global.h"
#include "random.h"
#include "error.h"
#include "field_ordering.h"
#include "spinor_field.h"
#include "communications.h"
#include <math.h>

void random_u(suNg_field *gf)
{
  _DECLARE_INT_ITERATOR(ix);

  error(gf==NULL,1,"random_u [random_fields.c]",
	"Attempt to access unallocated memory space");   
  
  _MASTER_FOR(gf->type,ix) {
    /* unroll 4 directions */
    suNg *ptr=(gf->ptr)+coord_to_index(ix,0);
    random_suNg(ptr++);
    random_suNg(ptr++);
    random_suNg(ptr++);
    random_suNg(ptr);
  }

  start_gf_sendrecv(gf);

}

void unit_u(suNg_field *gf)
{
  _DECLARE_INT_ITERATOR(ix);
  suNg unity;

  error(gf==NULL,1,"unit_u [random_fields.c]",
	"Attempt to access unallocated memory space");   
  
   _suNg_unit(unity);
  _MASTER_FOR(gf->type,ix) {
    /* unroll 4 directions */
    suNg *ptr=(gf->ptr)+coord_to_index(ix,0);
    *(ptr++)=unity;
    *(ptr++)=unity;
    *(ptr++)=unity;
    *(ptr)=unity;
  }

  start_gf_sendrecv(gf);

}

void SF_gauge_bcs(suNg_field *gf, int strength)
{
  _DECLARE_INT_ITERATOR(i);
  suNg unity;
  _suNg_unit(unity);
  int ix, iy, iz;

  error(gf==NULL,1,"SF_gauge_bcs [random_fields.c]",
	"Attempt to access unallocated memory space");   

  /*Boundary gauge fields*/
  	double pi = 3.14159265358979323846264338327950288419;
	double eta; /*This paramterises the boundary fields*/
	suNg Bound0, BoundT;
   if(strength==1) /*SF bcs*/
   {
        if(NG==4)
        {
           /*SU(4) specific fields:*/
           eta = 0.0;
           _suNg_zero(Bound0);
           ((Bound0).c[0]).re = cos((-0.5*eta - (sqrt(2)*pi/4.0))/(double)GLB_X);
           ((Bound0).c[0]).im = sin((-0.5*eta - (sqrt(2)*pi/4.0))/(double)GLB_X);
           ((Bound0).c[5]).re = cos((-0.5*eta - ((2-sqrt(2))*pi/4.0))/(double)GLB_X);
           ((Bound0).c[5]).im = sin((-0.5*eta - ((2-sqrt(2))*pi/4.0))/(double)GLB_X);
           ((Bound0).c[10]).re = cos(( 0.5*eta + ((2-sqrt(2))*pi/4.0))/(double)GLB_X);
           ((Bound0).c[10]).im = sin(( 0.5*eta + ((2-sqrt(2))*pi/4.0))/(double)GLB_X);
           ((Bound0).c[15]).re = cos(( 0.5*eta + (sqrt(2)*pi/4.0))/(double)GLB_X);
       	   ((Bound0).c[15]).im = sin(( 0.5*eta + (sqrt(2)*pi/4.0))/(double)GLB_X);

           _suNg_zero(BoundT);
           ((BoundT).c[0]).re = cos((0.5*eta - ((2+sqrt(2))*pi/4.0))/(double)GLB_X);
           ((BoundT).c[0]).im = sin((0.5*eta - ((2+sqrt(2))*pi/4.0))/(double)GLB_X);
           ((BoundT).c[5]).re = cos((0.5*eta - ((4-sqrt(2))*pi/4.0))/(double)GLB_X);
           ((BoundT).c[5]).im = sin((0.5*eta - ((4-sqrt(2))*pi/4.0))/(double)GLB_X);
           ((BoundT).c[10]).re = cos((- 0.5*eta + ((4-sqrt(2))*pi/4.0))/(double)GLB_X);
       	   ((BoundT).c[10]).im = sin((- 0.5*eta + ((4-sqrt(2))*pi/4.0))/(double)GLB_X);
           ((BoundT).c[15]).re = cos((- 0.5*eta + ((2+sqrt(2))*pi/4.0))/(double)GLB_X);
           ((BoundT).c[15]).im = sin((- 0.5*eta + ((2+sqrt(2))*pi/4.0))/(double)GLB_X);
        }
	else if(NG==3)
	{
	   /*SU(3) specific fields:*/
	   eta = 0.0;
	   _suNg_zero(Bound0);
	   ((Bound0).c[0]).re = cos((eta - (pi/3.0))/(double)GLB_X);
	   ((Bound0).c[0]).im = sin((eta - (pi/3.0))/(double)GLB_X);
	   ((Bound0).c[4]).re = cos(- 0.5*eta/(double)GLB_X);
	   ((Bound0).c[4]).im = sin(- 0.5*eta/(double)GLB_X);
	   ((Bound0).c[8]).re = cos((- 0.5*eta + (pi/3.0))/(double)GLB_X);
	   ((Bound0).c[8]).im = sin((- 0.5*eta + (pi/3.0))/(double)GLB_X);
	
	   _suNg_zero(BoundT);
	   ((BoundT).c[0]).re = cos((- eta - pi)/(double)GLB_X);
	   ((BoundT).c[0]).im = sin((- eta - pi)/(double)GLB_X);
	   ((BoundT).c[4]).re = cos((0.5*eta + (pi/3.0))/(double)GLB_X);
	   ((BoundT).c[4]).im = sin((0.5*eta + (pi/3.0))/(double)GLB_X);
	   ((BoundT).c[8]).re = cos((0.5*eta + (2.0*pi/3.0))/(double)GLB_X);
	   ((BoundT).c[8]).im = sin((0.5*eta + (2.0*pi/3.0))/(double)GLB_X);
	}
	else if(NG==2)
	{
	   /*SU(2) specific fields*/
	   eta = pi/4.0;
	   _suNg_zero(Bound0);
	   ((Bound0).c[0]).re = cos(-eta/(double)GLB_X);
	   ((Bound0).c[0]).im = sin(-eta/(double)GLB_X);
	   ((Bound0).c[3]).re = cos(eta/(double)GLB_X);
	   ((Bound0).c[3]).im = sin(eta/(double)GLB_X);
	
	   _suNg_zero(BoundT);
	   ((BoundT).c[0]).re = cos((eta-pi)/(double)GLB_X);
	   ((BoundT).c[0]).im = sin((eta-pi)/(double)GLB_X);
	   ((BoundT).c[3]).re = cos((pi-eta)/(double)GLB_X);
	   ((BoundT).c[3]).im = sin((pi-eta)/(double)GLB_X);
	}
	else
	{
	   error(gf==NULL,1,"SF_gauge_bcs [random_fields.c]",
	         "No SF boundary gauge fields defined for NG");
	}
   }
   else /*UNIT bcs*/
   {
	_suNg_unit(Bound0);
	_suNg_unit(BoundT);
   }	
  
if(COORD[0]==0)
{
  _MASTER_FOR(gf->type,i) {
		       	for (ix=0; ix<GLB_X/NP_X; ++ix)
		        for (iy=0; iy<GLB_Y/NP_Y; ++iy)
		        for (iz=0; iz<GLB_Z/NP_Z; ++iz)
		        {
			{
			{
			if (ipt(0,ix,iy,iz)==i)
			{
			    suNg *ptr=(gf->ptr)+coord_to_index(i,0);
			    *(ptr++)=unity;
			    *(ptr++)=unity;
			    *(ptr++)=unity;
			    *(ptr)=unity;
			}
			if (ipt(1,ix,iy,iz)==i)
			{
			    suNg *ptr=(gf->ptr)+coord_to_index(i,0);
			    ptr++; /*don't change timelike link*/
			    *(ptr++)=Bound0;
			    *(ptr++)=Bound0;
			    *(ptr)=Bound0;
			}
			}
			}
			}
  }
}
if(COORD[0]==NP_T-1)
{
  _MASTER_FOR(gf->type,i) {
		       	for (ix=0; ix<GLB_X/NP_X; ++ix)
		        for (iy=0; iy<GLB_Y/NP_Y; ++iy)
		        for (iz=0; iz<GLB_Z/NP_Z; ++iz)
		        {
			{
			{
			if (ipt((GLB_T/NP_T)-1,ix,iy,iz)==i)
			{
			    suNg *ptr=(gf->ptr)+coord_to_index(i,0);
			    *(ptr++)=unity;
			    *(ptr++)=BoundT;
			    *(ptr++)=BoundT;
			    *(ptr)=BoundT;
			}
			}
			}
			}
  }
}

  start_gf_sendrecv(gf);

}

