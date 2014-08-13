/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef UPDATE_H
#define UPDATE_H

#include "suN.h"
#include "inverters.h"
#include "rational_functions.h"

void staples(int ix,int mu,suNg *v);
void test_staples();
void cabmar(double beta,suNg *u, suNg *v,int type);

void project_gauge_field(void);

void update(double beta,int nhb,int nor);
void random_su2(double rho,double s[]);

/* forces for the update */
void force0(double dt, suNg_av_field *force, void *par);

typedef struct {
  int n_pf;
  spinor_field *pf;
  double mass;
  rational_app *ratio;
  double inv_err2;
} force_rhmc_par;
void init_force_rhmc();
void free_force_rhmc();
void force_rhmc(double dt, suNg_av_field *force, void *par);

typedef struct {
  int n_pf;
  spinor_field *pf;
  int hasenbusch; /* 0 = force with Y = Ddag^{-1} phi (standard) ; 2 = force with Y = Ddag^{-1}(phi+dm*X) */
  double mass;
  double b;
  double inv_err2, inv_err2_flt;
} force_hmc_par;
void init_force_hmc();
void free_force_hmc();
void force_hmc(double dt, suNg_av_field *force, void *par);



typedef struct _integrator_par {
  int nsteps;
  void (*force)(double,suNg_av_field*,void*);
  void *force_par;
  void (*integrator)(suNg_av_field*, double, struct _integrator_par*);
  struct _integrator_par *next;
  int level;
} integrator_par;
void gauge_integrator(suNg_av_field *momenta, double tlen, integrator_par *int_par);
void O2MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par);
void O4MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par);



void gaussian_momenta(suNg_av_field *momenta);
void gaussian_spinor_field(spinor_field *s);
void gaussian_spinor_field_flt(spinor_field_flt *s);
void z2_spinor_field(spinor_field *s);


typedef struct _puregauge_par {
  /* sim parameters */
  double beta;
    
  double MT_prec; /* metropolis test precision */
  double MD_prec; /* molecular dynamics precision */
  
  double tlen; /* trajectory lenght */
  unsigned int gsteps; /* number of substeps for the gauge part every step */
} puregauge_par;
void init_puregauge(puregauge_par *par);
void free_puregauge();


typedef struct _rhmc_par {
  /* sim parameters */
  double beta;
  int nf;
  double mass;
  
  double theta[4];
  
  double SF_zf;
  double SF_ds;
  int SF_sign;
  double SF_ct;
  double SF_prec;
  
  double MT_prec; /* metropolis test precision */
  double MD_prec; /* molecular dynamics precision */
  double HB_prec; /* heatbath precision for pseudofermions */
  double force_prec; /* precision used in the inversions in the force */
  double force_prec_flt; /* precision used in single-precision acceleration for the inversions in the force (HMC only)*/
  unsigned int n_pf; /* number of psudofermions used in the evolution */
  
  double tlen; /* trajectory lenght */
  unsigned int nsteps; /* number of step in the integration */
  unsigned int gsteps; /* number of substeps for the gauge part every step */
} rhmc_par;
void init_rhmc(rhmc_par *par);
void free_rhmc();

typedef struct _hmc_par {
  /* sim parameters */
  double beta;
  int nf;
  double mass;
  
  double theta[4];

  double SF_zf;
  double SF_ds;
  int SF_sign;
  double SF_ct;
	
  double n_MT_prec; /* n metropolis test precision */
  double h_MT_prec; /* h metropolis test precision */
  double n_MT_prec_flt; /* n metropolis test precision float*/
  double h_MT_prec_flt; /* h metropolis test precision float*/
  double MD_prec; /* molecular dynamics precision */
  double HB_prec; /* heatbath precision for pseudofermions */
  double n_force_prec; /* precision used in the inversions in the force F1 (higher mass) of the Hasenbush acceleration*/
  double n_force_prec_flt; /* precision used in single-precision acceleration for the inversions in the force F1 (HMC only)*/
  double h_force_prec; /* precision used in the inversions in the force F2 (lower mass) of the Hasenbush acceleration*/
  double h_force_prec_flt; /* precision used in single-precision acceleration for the inversions in the force F2 (HMC only)*/

  double tlen; /* trajectory lenght */
  unsigned int nsteps; /* number of steps in the integration */
  unsigned int hsteps; /* number of substeps for the lower mass of the Hasenbush acceleration */
  unsigned int gsteps; /* number of substeps for the gauge part every step */
  
  int hasenbusch; /* 0=no hasenbusch ; 1=hasenbusch */
  double hasen_dm;
} hmc_par;

void init_hmc(hmc_par *par);
void free_hmc();

/* update the gauge field using RHMC algorithm
 * return code: (<0 means an error has occurred)
 * 0 => conf correctly generated but has not passed metropolis test
 * 1 => conf accepted
 * 
 * -1 => rhmc has not been initialized. call init_rhmc first.
 */
int update_rhmc();
int update_hmc();
int update_puregauge();
int update_rhmc_o();

/* this is the basic operator used in the update */
/* defined in update_rhmc.c */
void H2(spinor_field *out, spinor_field *in);
void H(spinor_field *out, spinor_field *in);
void H_flt(spinor_field_flt *out, spinor_field_flt *in);


/* local action */
typedef enum {
   NEW=1,
   DELTA=2
} local_action_type;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
typedef struct _action_par {
  /* sim parameters */
  double beta;
  int n_pf;
#ifdef ROTATED_SF
  double SF_ct;
#endif
} action_par;
void local_hmc_action(local_action_type type,
                      action_par *par,
                      scalar_field *loc_action,
                      suNg_av_field *momenta,
                      spinor_field *phi1,
                      spinor_field *phi2);

void suNg_field_copy(suNg_field *g1, suNg_field *g2);
void suNf_field_copy(suNf_field *g1, suNf_field *g2);

/* find spectral interval using eva */
void find_spec_H2(double *max, double *min, double mass);

#endif
