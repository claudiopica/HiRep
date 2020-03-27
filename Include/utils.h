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
#include <stdlib.h>


/*SUN exp matrix*/

void WF_Exp(suNg *u, suNg *X);
void WF_Exp_Taylor(suNg *u, suNg *X);
void ExpX(double dt, suNg_algebra_vector *h, suNg *u);

void vector_star(suNg_vector *, suNg_vector *);

typedef struct
{
  double gauge_boundary_improvement_cs;
  double gauge_boundary_improvement_ct;
  double chiSF_boundary_improvement_ds;
  double fermion_twisting_theta[4];
  int SF_BCs;
  suNg gauge_boundary_up;
  suNg gauge_boundary_dn;
} BCs_pars_t;

void init_BCs(BCs_pars_t *pars);
void init_plaq_open_BCs(double * plaq_weight,double * rect_weight,double ct, double cs);

void free_BCs();
void apply_BCs_on_represented_gauge_field();
void apply_BCs_on_fundamental_gauge_field();
void apply_BCs_on_momentum_field(suNg_av_field *force);
void apply_BCs_on_spinor_field(spinor_field *sp);
void apply_BCs_on_spinor_field_flt(spinor_field_flt *sp);
void apply_background_field_zdir(suNg_field *V, double Q, int n);
void apply_BCs_on_clover_term(suNfc_field *);

void init_pure_gauge_anisotropy(double *chi);


inline int safe_mod(int x,int y)
{
   if (x>=0)
      return(x%y);
   else
      return((y-(abs(x)%y))%y);
}

/*Global shift for fields, the routine accepts also NULL entries in which case it does nothing*/
void shift_fields(int *shift, spinor_field *sin, suNg_field *uin, spinor_field *sout, suNg_field *uout);


void cross_prod(suNg_vector *v1, suNg_vector *v2, suNg_vector *v3);
void cross_prod_flt(suNg_vector_flt *v1, suNg_vector_flt *v2, suNg_vector_flt *v3);
void project_to_suNg(suNg *u);
void project_to_suNg_flt(suNg_flt *u);
void project_cooling_to_suNg(suNg *g_out, suNg *g_in, int cooling);
void covariant_project_to_suNg(suNg *u);

#ifndef GAUGE_SON
void ludcmp(double complex *a, int *indx, double *d, int N);
void lubksb(double complex *a, int *indx, double complex *b, int N);
void inv_hermNg(suNg *a);
void det_hermNg(double complex *res, suNg *a);
#else
int project_to_suNg_real(suNg *out, suNg *in);
void inv_hermNg(suNg *a);
void det_hermNg(double *res, suNg *a);
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
typedef struct _eva_prec
{
  /* EVA parameters */
  int nevt;      /* search space dimension */
  int nev;       /* number of accurate eigenvalues */
  int kmax;      /* max degree of polynomial */
  int maxiter;   /* max number of subiterations */
  double omega1; /* absolute precision */
  double omega2; /* relative precision */
} eva_prec;
void set_def_matrix(eva_prec *e_par, spinor_operator H, geometry_descriptor *type);
void eva_def(spinor_field *out, spinor_field *in);
void eva_def_inv(spinor_field *out, spinor_field *in, double m);

/* HYP smearing */
void spatialHYP_smearing(suNg_field *out, suNg_field *in, double weight[3]);
void HYP_smearing(suNg_field *out, suNg_field *in, double weight[3]);
double min_tplaq(suNg_field *g);
void HYP_span_parameters(double mtp[6859]);
int HYP_best_parameters(double mtp[6859], double w[3]);

/* Timing */
#include <sys/time.h>
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

/* CINFO */
void print_compiling_info();
void print_compiling_info_short();

/* Spatial Trasformations*/
void initialize_spatial_active_slices(int *tlist);
void free_spatial_active_slices();
#ifdef MAIN_PROGRAM
int *active_slices_list = NULL;
int *glbT_to_active_slices = NULL;
int n_active_slices;
#else
extern int *active_slices_list;
extern int *glbT_to_active_slices;
extern int n_active_slices;
#endif

/* Spatial blocking */
typedef enum
{
  NEW_SBLK = 1,
  CONT_SBLK = 0
} eval_spat_block;

int spatial_blocking_wrkspace(eval_spat_block eval, unsigned int level);
int single_level_spatial_blocking_wrkspace(int wrk_in);

/* Spatial rotation*/
void assign_spatial_rotated_wrkspace(int *map, int idx_wrkspace);

/* Spatial APE smearing*/
int spatial_APE_smear_wrkspace(double *smear_val, int wrkspace_in);

/* Workspace database*/
int iup_wrk(int site, int dir);
int idn_wrk(int site, int dir);
suNg *pu_gauge_wrk(int site, int dir);
suNg_field *u_gauge_wrk();
void reset_wrk_pointers();
void set_wrk_space(int i);
void set_wrk_space_and_pointers(int i, suNg_field **g_wrk_out, int **i_up_wrk_out, int **i_dn_wrk_out);
int reserve_wrk_space();
int reserve_wrk_space_with_pointers(suNg_field **g_wrk_out, int **i_up_wrk_out, int **i_dn_wrk_out);
void release_wrk_space(int id_release);
void free_wrk_space();

#endif
