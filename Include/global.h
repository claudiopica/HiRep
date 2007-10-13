/*******************************************************************************
*
* File global.h
*
* Global parameters and arrays
*
*******************************************************************************/

#ifndef GLOBAL_H
#define GLOBAL_H

#include <stddef.h>

#define VOLUME (T*L*L*L)
#define VOL3   (L*L*L)

#if ((T%2)!=0)
#error : The temporal lattice size T must be even
#endif

#if ((T<2)||(L<2))
#error : The lattice size should be at least 2 in all directions
#endif

#ifdef MAIN_PROGRAM
#  define GLB_PTR(type,name) type *name=NULL;
#define EXTERN
#else
#  define GLB_PTR(type,name) extern type *name;
#define EXTERN extern
#endif

/* Geometry indexes */
GLB_PTR(int, ipt);
GLB_PTR(int, ipt_4d);
GLB_PTR(int, iup);
GLB_PTR(int, idn);

#define ipt(t,x,y,z) ipt[(t)*(VOL3)+(x)*(L*L)+(y)*(L)+(z)]
#define ipt_4d(t,x) ipt_4d[(t)*(VOL3)+(x)]
#define iup(site,dir) iup[(site)*4+(dir)]
#define idn(site,dir) idn[(site)*4+(dir)]

/* Gauge field */
#include "suN_types.h"

GLB_PTR(suNg, u_gauge);
GLB_PTR(suNg_flt, u_gauge_flt);
GLB_PTR(suNf, u_gauge_f);
GLB_PTR(suNf_flt, u_gauge_f_flt);

#define coord_to_index(ix,mu) ((ix)*4+(mu))
#define index_to_coord(i,ix,mu) (mu)=(i&3);(ix)=(i>>2)

#define pu_gauge(ix,mu) (u_gauge+coord_to_index(ix,mu))
#define pu_gauge_flt(ix,mu) (u_gauge_flt+coord_to_index(ix,mu))
#define pu_gauge_f(ix,mu) (u_gauge_f+coord_to_index(ix,mu))
#define pu_gauge_f_flt(ix,mu) (u_gauge_f_flt+coord_to_index(ix,mu))

#undef GLB_PTR
#undef EXTERN
#undef INIT

#endif

