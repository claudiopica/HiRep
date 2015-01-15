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

/* functions and structures for the MRE algorithm */
typedef struct {
	spinor_field *s[2];
	int num[2];
	int max;
	int init;
} mre_par;

void mre_guess(mre_par*, int, spinor_field*, spinor_operator, spinor_field*);
void mre_store(mre_par*, int, spinor_field*);
void mre_init(mre_par*, int, double);

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
  int id;
  int n_pf;
  spinor_field *pf;
  int hasenbusch;
  double mass;
  double b;
  double mu;
  double inv_err2, inv_err2_flt;
  mre_par mpar;
} force_hmc_par;

void force_fermion_core(spinor_field* Xs, spinor_field* Ys, suNg_av_field* force, double dt, double* forcestat, int type);
void force_hmc(double dt, suNg_av_field *force, void *par);
void force_hmc_tm(double dt, suNg_av_field *force, void *par);


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
  unsigned int* hsteps; /* List of number of substeps for the heavier mass of the Hasenbush acceleration */
  unsigned int gsteps; /* number of substeps for the gauge part every step */
  
  int hasenbusch; /* 0=no hasenbusch ; 1=hasenbusch */
  int n_hasen; /* Number of Hasenbusch levels*/
  double* hasen_dm; /* List of differences of heavier mass of the Hasenbush acceleration */
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

void corret_pf_dist_hmc();
void calc_one_force(int n_force);


/* Action structures */
typedef enum {
  PureGauge,
  HMC,
  RHMC,
  TM,
  TM_alt,
  Hasenbusch,
  Hasenbusch_tm,
  Hasenbusch_tm_alt
} mon_type;

typedef struct _mon_pg_par {
  double beta;
} mon_pg_par;

typedef struct _mon_hmc_par {
  double mass;
  int mre_past;
  force_hmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_hmc_par;

typedef struct _mon_rhmc_par {
  double mass;
  rational_app ratio;
  force_rhmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_rhmc_par;

typedef struct _mon_tm_par {
  double mass;
  double mu;
  int mre_past;
  force_hmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_tm_par;


typedef struct _mon_hasenbusch_par {
  double mass;
  double dm;
  int mre_past;
  force_hmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_hasenbusch_par;

typedef struct _mon_hasenbusch_tm_par {
  double mass;
  double mu;
  double dmu;
  int mre_past;
  force_hmc_par fpar;
  spinor_field *pf; /* pseudofermion field */
} mon_hasenbusch_tm_par;


typedef struct _monomial {
  int id; /* monomial id */
  mon_type type; /* type of monomial */
  void *par; /* parameters */
  double MT_prec; /* metropolis precision */
  double MD_prec; /* molecular dynamics precision */
  double force_prec; /* force precision */
} monomial;

const monomial *add_mon(monomial *mon);
int num_mon();
const monomial *mon_n(int i);


typedef struct _integrator_par {
  int nsteps;
  int nmon;
  const monomial **mon_list;
  void (*integrator)(suNg_av_field*, double, struct _integrator_par*);
  struct _integrator_par *next;
  int level;
} integrator_par;

void gauge_integrator(suNg_av_field *momenta, double tlen, integrator_par *int_par);
void leapfrog_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par);
void O2MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par);
void O4MN_multistep(suNg_av_field *momenta, double tlen, integrator_par *int_par);


typedef struct _ghmc_par {
	/* action parameter */
 // mon_list *action;
  
  /* integrator */
  integrator_par *integrator;
  double tlen;

  /* Fermion Theta angles */
  double theta[4];
  
  /* SF stuff */
  double SF_zf;
  double SF_ds;
  int SF_sign;
  double SF_ct;
  
  /* other parameters */
  int mre_past; /* chronological guesser */

  /* dummy variables */
  /* NEED TO DELETE */
  double beta;
  int nf;
  double mass;

  
} ghmc_par;

void init_ghmc(ghmc_par *par);
void free_ghmc();
int update_ghmc();



/* local action */
typedef enum {
   NEW=1,
   DELTA=2
} local_action_type_old;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
typedef struct _action_par_old {
  /* sim parameters */
  double beta;
  int n_pf;
#ifdef ROTATED_SF
  double SF_ct;
#endif
} action_par_old;


void local_hmc_action(local_action_type_old type,
                      scalar_field *loc_action,
                      suNg_av_field *momenta);

void suNg_field_copy(suNg_field *g1, suNg_field *g2);
void suNf_field_copy(suNf_field *g1, suNf_field *g2);

/* find spectral interval using eva */
void find_spec_H2(double *max, double *min);

#endif
