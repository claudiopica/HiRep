/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

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
#  define GLB_VAR(type,name,init...) type name init
#else
#  define GLB_VAR(type,name,init...) extern type name
#endif

/* local lattice attributes */
GLB_VAR(int,T,=0); /* local lattice size in direction T */ 
GLB_VAR(int,X,=0); /* local lattice size in direction X */
GLB_VAR(int,Y,=0); /* local lattice size in direction Y */
GLB_VAR(int,Z,=0); /* local lattice size in direction Z */
/* this two probably are not more needed... */
GLB_VAR(int,VOL3,=0); 
GLB_VAR(int,VOLUME,=0);

/* Nodes attributes
 * NP = number of processes in each direction 
 * 1 => local direction
 */
GLB_VAR(int,NP_T,=1); /* number of processes in direction T */
GLB_VAR(int,NP_X,=1); /* number of processes in direction X */
GLB_VAR(int,NP_Y,=1); /* number of processes in direction Y */
GLB_VAR(int,NP_Z,=1); /* number of processes in direction Z */

/* global lattice attributes */
GLB_VAR(int,GLB_T,=0); /* global size of the lattice in direction T */
GLB_VAR(int,GLB_X,=0); /* global size of the lattice in direction X */
GLB_VAR(int,GLB_Y,=0); /* global size of the lattice in direction Y */
GLB_VAR(int,GLB_Z,=0); /* global size of the lattice in direction Z */

GLB_VAR(int,T_BORDER,=0);
GLB_VAR(int,X_BORDER,=0);
GLB_VAR(int,Y_BORDER,=0);
GLB_VAR(int,Z_BORDER,=0);

GLB_VAR(int,T_EXT,=0);
GLB_VAR(int,X_EXT,=0);
GLB_VAR(int,Y_EXT,=0);
GLB_VAR(int,Z_EXT,=0);

/* MPI stuff */
GLB_VAR(int,WORLD_SIZE,=1); /* mpi rank for this process */
GLB_VAR(int,CART_SIZE,=1); /* mpi rank for this process */
#ifdef WITH_MPI
#include <mpi.h>
GLB_VAR(MPI_Comm,cart_comm,=MPI_COMM_NULL);
#endif

/* ID for this process */
GLB_VAR(int,PID,=0); /* ID of this process */

GLB_VAR(int,CID,=0); /* cartesian ID for this process */
GLB_VAR(int,COORD[4],={0}); /* cartesian coordinates for this process */
GLB_VAR(int,PSIGN,=0); /* parity of this process */

/* Geometry indexes */
GLB_VAR(int,*ipt, =NULL);
GLB_VAR(int,*ipt_4d,=NULL);
GLB_VAR(int,*iup,=NULL);
GLB_VAR(int,*idn,=NULL);

/* Geometry structures */
#define ipt(t,x,y,z) ipt[((((t)+T_BORDER)*(X_EXT)+((x)+X_BORDER))*(Y_EXT)+((y)+Y_BORDER))*(Z_EXT)+((z)+Z_BORDER)]
#define ipt_ext(t,x,y,z) ipt[(((t)*(X_EXT)+(x))*(Y_EXT)+(y))*(Z_EXT)+(z)]
#define ipt_4d(t,x) ipt_4d[(t)*(VOL3)+(x)]
#define iup(site,dir) iup[(site)*4+(dir)]
#define idn(site,dir) idn[(site)*4+(dir)]

/* Geometry structures */
#include "geometry.h"

GLB_VAR(geometry_descriptor,glattice,={0}); /* global lattice */
GLB_VAR(geometry_descriptor,glat_even,={0}); /* global even lattice */
GLB_VAR(geometry_descriptor,glat_odd,={0}); /* global odd lattice */
GLB_VAR(mem_t,std_mem_t, /* default memory allocation type for fields */
#ifdef WITH_GPU
					=(CPU_MEM | GPU_MEM)
#else
        	=(CPU_MEM)
#endif
        );
GLB_VAR(mem_t,alloc_mem_t, std_mem_t); /* memory type requested for allocating fields */


/* Gauge field */
#include "field_ordering.h"
#include "suN_types.h"
#include "spinor_field.h"

GLB_VAR(suNg_field,*u_gauge,=NULL);
GLB_VAR(suNg_field_flt,*u_gauge_flt,=NULL);
GLB_VAR(suNf_field,*u_gauge_f,=NULL);
GLB_VAR(suNf_field_flt,*u_gauge_f_flt,=NULL);

#define pu_gauge(ix,mu) ((u_gauge->ptr)+coord_to_index(ix,mu))
#define pu_gauge_flt(ix,mu) ((u_gauge_flt->ptr)+coord_to_index(ix,mu))
#define pu_gauge_f(ix,mu) ((u_gauge_f->ptr)+coord_to_index(ix,mu))
#define pu_gauge_f_flt(ix,mu) ((u_gauge_f_flt->ptr)+coord_to_index(ix,mu))


/* Boundary conditions */
#ifdef ANTIPERIODIC_BC_T
#define BC_T 1.
#else
#define BC_T 0.
#endif

#ifdef ANTIPERIODIC_BC_X
#define BC_X 1.
#else
#define BC_X 0.
#endif

#ifdef ANTIPERIODIC_BC_Y
#define BC_Y 1.
#else
#define BC_Y 0.
#endif

#ifdef ANTIPERIODIC_BC_Z
#define BC_Z 1.
#else
#define BC_Z 0.
#endif

GLB_VAR(double,bc[4],={BC_T,BC_X,BC_Y,BC_Z});

#undef BC_T
#undef BC_x
#undef BC_Y
#undef BC_Z

/* input parameters */
#include "input_par.h"
GLB_VAR(input_glb,glb_var,=init_input_glb(glb_var));



#ifdef WITH_GPU
#define BLOCK_SIZE 256
GLB_VAR(int,*iup_gpu,=NULL);
GLB_VAR(int,*idn_gpu,=NULL);
#endif //WITH_GPU

#undef GLB_VAR
#endif


