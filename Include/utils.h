/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File utils.h
* 
* Some useful functions
*
*******************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include "suN_types.h"
#include "spinor_field.h"
#include "inverters.h"
#include "geometry.h"
#include "gpu.h"

void ExpX(double dt, suNg_algebra_vector *h, suNg *u);

void apply_bc();
void apply_bc_flt();

#ifdef BASIC_SF
double SF_test_spinor_bcs(spinor_field *sp);
#endif /* BASIC_SF */

#if defined(BASIC_SF) || defined(ROTATED_SF)
void SF_spinor_bcs(spinor_field *sp);
void SF_spinor_bcs_flt(spinor_field_flt *sp);
void SF_force_bcs(suNg_av_field *force);
double SF_test_force_bcs(suNg_av_field *force);

void SF_gauge_bcs(suNg_field *gf, int strength);
double SF_test_gauge_bcs();
#endif /* BASIC_SF || ROTATED_SF */

#ifdef TWISTED_BC
void init_twbc();
void free_twbc();
#endif


void cross_prod(suNg_vector *v1,suNg_vector *v2,suNg_vector *v3);
void cross_prod_flt(suNg_vector_flt *v1,suNg_vector_flt *v2,suNg_vector_flt *v3);
void project_to_suNg(suNg *u);
void project_to_suNg_flt(suNg_flt *u);
void project_cooling_to_suNg(suNg* g_out, suNg* g_in, int cooling);



#ifdef __cplusplus
extern "C" {
#endif
void assign_u2ud_cpu(void);
void assign_ud2u_cpu(void);
void assign_ud2u_f_cpu(void);
void assign_u2ud_f_cpu(void);
void assign_s2sd_cpu(spinor_field *out, spinor_field_flt *in);
void assign_sd2s_cpu(spinor_field_flt *out, spinor_field *in);

#ifdef __cplusplus
}
#endif

#ifdef WITH_GPU
#ifdef __cplusplus
extern "C" {
#endif
void assign_u2ud(void);
void assign_ud2u(void);
void assign_ud2u_f(void);
void assign_u2ud_f(void);

void assign_s2sd(spinor_field *out, spinor_field_flt *in);
void assign_sd2s(spinor_field_flt *out, spinor_field *in);
#ifdef __cplusplus
}
#endif
#else
extern void (*assign_u2ud) (void);
extern void (*assign_ud2u) (void);
extern void (*assign_ud2u_f) (void);
extern void (*assign_u2ud_f) (void);
extern void (*assign_s2sd) (spinor_field *out, spinor_field_flt *in);
extern void (*assign_sd2s) (spinor_field_flt *out, spinor_field *in);
#endif

/* use power method to find max eigvalue of H2 */
int max_H(spinor_operator H, geometry_descriptor *type, double *max);

/* EVA preconditioning */
typedef struct _eva_prec {
  /* EVA parameters */
  int nevt; /* search space dimension */
  int nev; /* number of accurate eigenvalues */
  int kmax; /* max degree of polynomial */
  int maxiter; /* max number of subiterations */
  double omega1; /* absolute precision */
  double omega2; /* relative precision */
} eva_prec;
void set_def_matrix(eva_prec *e_par, spinor_operator H, geometry_descriptor *type);
void eva_def(spinor_field *out, spinor_field *in);
void eva_def_inv(spinor_field *out, spinor_field *in, double m);



/* HYP smearing */
void spatialHYP_smearing(suNg_field* out, suNg_field* in, double weight[3]);
void HYP_smearing(suNg_field* out, suNg_field* in, double weight[3]);
double min_tplaq(suNg_field* g);
void HYP_span_parameters(double mtp[6859]);
int HYP_best_parameters(double mtp[6859], double w[3]);


/* Timing */
#include <sys/time.h>
int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y);

/* gpu Timer */
#ifdef WITH_GPU
typedef cudaEvent_t gpu_timer;
cudaEvent_t gpuTimerStart();
float gpuTimerMeasure(cudaEvent_t e);
float gpuTimerStop(cudaEvent_t e);

#ifdef __cplusplus
extern "C" {
#endif
  void init_gpu(input_gpu gpu_var);
#ifdef __cplusplus
}
#endif
#endif 


#endif 
