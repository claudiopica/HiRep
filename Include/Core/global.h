/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file global.h
 * @brief Global parameters and arrays
 */

#ifndef GLOBAL_H
#define GLOBAL_H

//this file needs check_options to be included
//because check_options makes some other definitions
//based on the compilation options
#include "check_options.h"
#include "core_utils.h"

#include "geometry.h"
#include <stddef.h>
#include "flags.h"

#ifndef M_PI
// Define M_PI if its not part of the C standard
#define M_PI 3.14159265358979323846264338327
#endif

// Allow also to use PI for M_PI
#define PI M_PI

/* local lattice attributes */
GLB_VAR(int, T, = 0); /* local lattice size in direction T */
GLB_VAR(int, X, = 0); /* local lattice size in direction X */
GLB_VAR(int, Y, = 0); /* local lattice size in direction Y */
GLB_VAR(int, Z, = 0); /* local lattice size in direction Z */
GLB_VAR(long int, GLB_VOL3, = 0);
GLB_VAR(long int, GLB_VOLUME, = 0);
/* this two probably are not more needed... */
GLB_VAR(long int, VOL3, = 0);
GLB_VAR(long int, VOLUME, = 0);

/* Nodes attributes
 * NP = number of processes in each direction
 * 1 => local direction
 */
GLB_VAR(int, NP_T, = 1); /* number of processes in direction T */
GLB_VAR(int, NP_X, = 1); /* number of processes in direction X */
GLB_VAR(int, NP_Y, = 1); /* number of processes in direction Y */
GLB_VAR(int, NP_Z, = 1); /* number of processes in direction Z */

/* global MPI arrangment attributes */
GLB_VAR(int, MPI_BLK_T, = 1); /* MPI arrangment size of the sub-block in direction T */
GLB_VAR(int, MPI_BLK_X, = 1); /* MPI arrangment size of the sub-block in direction X */
GLB_VAR(int, MPI_BLK_Y, = 1); /* MPI arrangment size of the sub-block in direction Y */
GLB_VAR(int, MPI_BLK_Z, = 1); /* MPI arrangment size of the sub-block in direction Z */

/* global lattice attributes */
GLB_VAR(int, GLB_T, = 0); /* global size of the lattice in direction T */
GLB_VAR(int, GLB_X, = 0); /* global size of the lattice in direction X */
GLB_VAR(int, GLB_Y, = 0); /* global size of the lattice in direction Y */
GLB_VAR(int, GLB_Z, = 0); /* global size of the lattice in direction Z */

GLB_VAR(int, T_BORDER, = 0);
GLB_VAR(int, X_BORDER, = 0);
GLB_VAR(int, Y_BORDER, = 0);
GLB_VAR(int, Z_BORDER, = 0);

GLB_VAR(int, T_EXT, = 0);
GLB_VAR(int, X_EXT, = 0);
GLB_VAR(int, Y_EXT, = 0);
GLB_VAR(int, Z_EXT, = 0);

/*path blocking size*/
GLB_VAR(int, PB_T, = 2);
GLB_VAR(int, PB_X, = 2);
GLB_VAR(int, PB_Y, = 2);
GLB_VAR(int, PB_Z, = 2);

/* MPI stuff */
GLB_VAR(int, WORLD_SIZE, = 1); /* mpi rank for this process */
GLB_VAR(int, CART_SIZE, = 1); /* mpi rank for this process */
GLB_VAR(int, N_REP, = 1); /* number of replicas*/
GLB_VAR(int, MPI_WORLD_SIZE, = 1); /* mpi rank for this process */
GLB_VAR(int, MPI_PID, = 0); /* mpi rank inside MPI_COMM_WORLD (unique across replicas) */
#ifdef WITH_MPI
#include <hr_mpi.h>
GLB_VAR(MPI_Comm, GLB_COMM, = MPI_COMM_WORLD); /* this is the global communicator for a replica */
GLB_VAR(MPI_Comm, cart_comm, = MPI_COMM_NULL); /* cartesian communicator for the replica */
#endif

GLB_VAR(int, RID, = 0); /* Replica ID of this process */
GLB_VAR(int, PID, = 0); /* Process ID inside a replica */
GLB_VAR(int, LID, = 0); /* Process ID inside a replica local to the node */

GLB_VAR(int, CID, = 0); /* Cartesian ID inside a replica */
GLB_VAR(int, COORD[4], = { 0 }); /* cartesian coordinates for this process */
GLB_VAR(int, PSIGN, = 0); /* parity of this process */

/* Geometry indexes */
GLB_VAR(int, *ipt, = NULL);
GLB_VAR(int, *ipt_4d, = NULL); //CP: is this used?
GLB_VAR(int, *iup, = NULL);
GLB_VAR(int, *idn, = NULL);
GLB_VAR(char, *imask, = NULL);
GLB_VAR(int, zerocoord[4], = { 0, 0, 0, 0 });
GLB_VAR(int, *timeslices, = NULL);

GLB_VAR(int, BLK_T, = 4);
GLB_VAR(int, BLK_X, = 4);
GLB_VAR(int, BLK_Y, = 4);
GLB_VAR(int, BLK_Z, = 4);

#define BLK_T_GPU 4
#define BLK_X_GPU 4
#define BLK_Y_GPU 4
#define BLK_Z_GPU 4

/* Geometry structures */
GLB_VAR(geometry_descriptor, glattice, = { 0 }); /* global lattice */
GLB_VAR(geometry_descriptor, glat_even, = { 0 }); /* global even lattice */
GLB_VAR(geometry_descriptor, glat_odd, = { 0 }); /* global odd lattice */
GLB_VAR(geometry_descriptor, glat_even_red, = { 0 });
GLB_VAR(geometry_descriptor, glat_odd_red, = { 0 });
GLB_VAR(geometry_descriptor, glat_even_black, = { 0 });
GLB_VAR(geometry_descriptor, glat_odd_black, = { 0 });
GLB_VAR(geometry_descriptor, glat_red, = { 0 });
GLB_VAR(geometry_descriptor, glat_black, = { 0 });

// new geometry structure
GLB_VAR(box_t *, geometryBoxes, );
GLB_VAR(box_t *, geometryBoxes_gpu, );
GLB_VAR(coord4 *, icoord_gpu, );
GLB_VAR(coord4 *, sb_icoord_gpu, );

#ifdef UPDATE_EO
#define glat_default glat_even
#else
#define glat_default glattice
#endif

#ifdef WITH_GPU
#define STD_MEM_TYPE (CPU_MEM | GPU_MEM)
#include "gpu.h"
#define BLOCK_SIZE 256
#define BLOCK_SIZE_LINEAR_ALGEBRA 256
#define BLOCK_SIZE_GLOBAL_SUM 512
#define BLOCK_SIZE_DIRAC 256
#define BLOCK_SIZE_CLOVER 256
#define BLOCK_SIZE_DIRAC_FLT 512
#define BLOCK_SIZE_SYNC 32

GLB_VAR(cudaStream_t, non_default_stream, = NULL);
GLB_VAR(cudaStream_t, memory_streams[16]);
GLB_VAR(kernel_field_input, **input, = NULL);

GLB_VAR(input_gpu, gpu_var, = init_input_gpu(gpu_var));
GLB_VAR(int, gpu_id, = 0);
GLB_VAR(int, *ipt_gpu, = NULL);
GLB_VAR(int, *iup_gpu, = NULL);
GLB_VAR(int, *idn_gpu, = NULL);
GLB_VAR(char, *imask_gpu, = NULL);
GLB_VAR(unsigned int, grid_size_max_gpu, = 65535);
GLB_VAR(int, *timeslices_gpu, = NULL);
#else
#define STD_MEM_TYPE (CPU_MEM)
#endif
GLB_VAR(mem_t, std_mem_t, = STD_MEM_TYPE); /* default memory allocation type for fields */
GLB_VAR(mem_t, alloc_mem_t, = STD_MEM_TYPE); /* memory type requested for allocating fields */

/* Communicate only the CPU or only the GPU copy depending 
   on compilation variables. This allows also to communicate
   both for testing purposes.*/
#ifndef WITH_GPU
GLB_VAR(comm_t, std_comm_t, = CPU_COMM);
#else
GLB_VAR(comm_t, std_comm_t, = GPU_COMM);
#endif

/* Gauge field */
#include "suN_types.h"
#include "spinor_field.h"

GLB_VAR(suNg_field, *u_gauge, = NULL);
GLB_VAR(suNg_scalar_field, *u_scalar, = NULL);
GLB_VAR(suNg_field_flt, *u_gauge_flt, = NULL);
GLB_VAR(suNf_field, *u_gauge_f, = NULL);
GLB_VAR(suNg_field, *u_gauge_s, = NULL);
GLB_VAR(suNf_field_flt, *u_gauge_f_flt, = NULL);
GLB_VAR(clover_term, *cl_term, = NULL);
#if defined(WITH_EXPCLOVER) && defined(WITH_GPU)
GLB_VAR(clover_term, *cl_term_expAplus, = NULL);
GLB_VAR(clover_term, *cl_term_expAminus, = NULL);
GLB_VAR(clover_term, *cl_term_expAplusinv, = NULL);
GLB_VAR(clover_term, *cl_term_expAminusinv, = NULL);
#endif
GLB_VAR(clover_force, *cl_force, = NULL);
GLB_VAR(ldl_field, *cl_ldl, = NULL);
GLB_VAR(suNg_av_field, *suN_momenta, = NULL);
GLB_VAR(suNg_scalar_field, *scalar_momenta, = NULL);
GLB_VAR(int, gauge_field_active, = 0); // whether gauge field interactions is active

#define pu_gauge(ix, mu) ((u_gauge->ptr) + coord_to_index(ix, mu))
#define pu_scalar(ix) ((u_scalar->ptr) + ix)
#define pu_gauge_flt(ix, mu) ((u_gauge_flt->ptr) + coord_to_index(ix, mu))
#define pu_gauge_f(ix, mu) ((u_gauge_f->ptr) + coord_to_index(ix, mu))
#define pu_gauge_f_flt(ix, mu) ((u_gauge_f_flt->ptr) + coord_to_index(ix, mu))

/* input parameters */
#include "IO/input_par.h"
GLB_VAR(input_glb, glb_var, = init_input_glb(glb_var));

/* Random number generator parameters */
GLB_VAR(input_rlx, rlx_var, = init_input_rlx(rlx_var));

/* logger parameters */
GLB_VAR(input_logger, logger_var, = init_input_logger(logger_var));

/* Does the represented field need to be allocated? */
#if (!defined(REPR_FUNDAMENTAL) && !defined(WITH_QUATERNIONS)) || defined(BC_T_SF_ROTATED)
#define ALLOCATE_REPR_GAUGE_FIELD
#endif

GLB_VAR(double, *plaq_weight, = NULL);
GLB_VAR(double, *rect_weight, = NULL);
GLB_VAR(double, *plaq_weight_gpu, = NULL);
GLB_VAR(double, *rect_weight_gpu, = NULL);

/* Theta Boundary conditions */
#ifdef FERMION_THETA
GLB_VAR(hr_complex, eitheta[4], = { 1.0, 1.0, 1.0, 1.0 });
#endif

/* Fields four fermion interactions */
/* Auxiliary fields for four fermion interactions */
GLB_VAR(scalar_field, *ff_sigma, = NULL);
GLB_VAR(scalar_field, *ff_pi, = NULL);
GLB_VAR(scalar_field, *ff_sigma_mom, = NULL);
GLB_VAR(scalar_field, *ff_pi_mom, = NULL);

GLB_VAR(int, four_fermion_active, = 0); // whether four fermion interactions are active

#ifndef THREADSIZE
#define THREADSIZE 32
#endif
#endif
