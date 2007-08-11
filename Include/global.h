/*******************************************************************************
*
* File global.h
*
* Global parameters and arrays
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#define VOLUME (T*L*L*L)
#define VOL3   (L*L*L)

#if ((T%2)!=0)
#error : The temporal lattice size T must be even
#endif

#if ((T<2)||(L<2))
#error : The lattice size should be at least 2 in all directions
#endif

#include "suN.h"

#ifdef MAIN_PROGRAM
#  define EXTERN
#else
#  define EXTERN extern
#endif

/* Geometry indexes */
EXTERN int ipt[T][L][L][L];
EXTERN int ipt_4d[T][VOL3];
EXTERN int iup[VOLUME][4];
EXTERN int idn[VOLUME][4];

/* Gauge field pointer */
EXTERN suNg *u_gauge;
EXTERN suNg_flt *u_gauge_flt;
EXTERN suNf *u_gauge_f;
EXTERN suNf_flt *u_gauge_f_flt;

#define gfield_ordering(ptr,ix,mu) (ptr+((ix)*4+(mu)))
#define index_to_coord(i,ix,mu) (mu)=i&3;(ix)=i>>2

#define pu_gauge(ix,mu) gfield_ordering(u_gauge,ix,mu)
#define pu_gauge_old(ix,mu) gfield_ordering(u_gauge_old,ix,mu)
#define pu_gauge_dble(ix,mu) gfield_ordering(u_gauge_dble,ix,mu)
#define pu_gauge_f(ix,mu) gfield_ordering(u_gauge_f,ix,mu)
#define pu_gauge_f_dble(ix,mu) gfield_ordering(u_gauge_f_dble,ix,mu)

#endif

