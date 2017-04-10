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

void ExpX(double dt, suNg_algebra_vector *h, suNg *u);

void vector_star(suNg_vector*,suNg_vector*);


typedef struct {
  double gauge_boundary_improvement_cs;
  double gauge_boundary_improvement_ct;
  double chiSF_boundary_improvement_ds;
  double fermion_twisting_theta[4];
  int SF_BCs;
  suNg gauge_boundary_up;
  suNg gauge_boundary_dn;
} BCs_pars_t;  

void init_BCs(BCs_pars_t *pars);
void free_BCs();
void apply_BCs_on_represented_gauge_field();
void apply_BCs_on_fundamental_gauge_field();
void apply_BCs_on_momentum_field(suNg_av_field *force);
void apply_BCs_on_spinor_field(spinor_field *sp);
void apply_BCs_on_spinor_field_flt(spinor_field_flt *sp);
void apply_background_field_zdir(suNg_field* V,double Q,int n);
void apply_BCs_on_clover_term(suNfc_field*);

void cross_prod(suNg_vector *v1,suNg_vector *v2,suNg_vector *v3);
void cross_prod_flt(suNg_vector_flt *v1,suNg_vector_flt *v2,suNg_vector_flt *v3);
void project_to_suNg(suNg *u);
void project_to_suNg_flt(suNg_flt *u);
void project_cooling_to_suNg(suNg* g_out, suNg* g_in, int cooling);

#ifndef GAUGE_SON
void ludcmp(complex* a, int* indx, double* d,int N);
void lubksb(complex* a, int* indx, complex* b,int N);
void inv_suNg(suNg* a);
void det_suNg(complex* res, suNg* a);
#else
int project_to_suNg_real(suNg *out, suNg *in);
void det_suNg(double* res, suNg *a);
void diag_hmat(suNg *hmat, double *dag);
#endif


void assign_u2ud(void);
void assign_ud2u(void);
void assign_ud2u_f(void);

/* void assign_s2sd(int len, suNf_spinor *out, suNf_spinor_flt *in); */
/* void assign_sd2s(int len, suNf_spinor_flt *out, suNf_spinor *in); */

void assign_s2sd(spinor_field *out, spinor_field_flt *in);
void assign_sd2s(spinor_field_flt *out, spinor_field *in);

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

#endif 
