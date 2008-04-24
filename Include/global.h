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

#ifdef MAIN_PROGRAM
#  define GLB_PTR(type,name) type *name=NULL
#  define GLB_VAR(type,name,value) type name=(value)
#  define GLB_STR(type,name) type name
#else
#  define GLB_PTR(type,name) extern type *name
#  define GLB_VAR(type,name,value) extern type name
#  define GLB_STR(type,name) extern type name
#endif

GLB_VAR(int,T,0);
GLB_VAR(int,X,0);
GLB_VAR(int,Y,0);
GLB_VAR(int,Z,0);
/*GLB_VAR(int,L,0);*/
GLB_VAR(int,VOL3,0);
GLB_VAR(int,VOLUME,0);

GLB_VAR(int,T_BORDER,1);
GLB_VAR(int,X_BORDER,1);
GLB_VAR(int,Y_BORDER,1);
GLB_VAR(int,Z_BORDER,1);

/*
#define VOLUME (T*L*L*L)
#define VOL3   (L*L*L)

#if ((T%2)!=0)
#error : The temporal lattice size T must be even
#endif

#if ((T<2)||(L<2))
#error : The lattice size should be at least 2 in all directions
#endif
*/

/* Geometry indexes */
GLB_PTR(int, ipt);
GLB_PTR(int, ipt_4d);
GLB_PTR(int, iup);
GLB_PTR(int, idn);

/* Geometry structures */
#define ipt(t,x,y,z) ipt[(t)*(VOL3)+(x)*(Y*Z)+(y)*(Z)+(z)]
#define ipt_4d(t,x) ipt_4d[(t)*(VOL3)+(x)]
#define iup(site,dir) iup[(site)*4+(dir)]
#define idn(site,dir) idn[(site)*4+(dir)]

/* Geometry structures */
#include "geometry.h"

GLB_STR(geometry_descriptor,glattice); /* global lattice */
GLB_STR(geometry_descriptor,glat_even); /* global even lattice */
GLB_STR(geometry_descriptor,glat_odd); /* global odd lattice */


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
#undef GLB_STR
#undef EXTERN
#undef INIT


#endif


