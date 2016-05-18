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
#ifdef TLSYM
void rect_staples_1x2(int ix,int mu,suNg *v);
void test_rect_staples_1x2();
#endif

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

typedef struct {
  int id;
  int n_pf;
  spinor_field *pf;
  double mass;
  rational_app *ratio;
  double inv_err2;
} force_rhmc_par;

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

//parameters for the four fermion force
typedef struct {
  double gamma;
} force_auxfield_par;


void force_measure_begin();
void force_measure_end(int, const char*, double, int);
void force_fermion_core(spinor_field*, spinor_field*, suNg_av_field*, int, double, double);

void force_hmc(double, suNg_av_field*, void*);
void force_hmc_tm(double, suNg_av_field*, void*);
void force_rhmc(double, suNg_av_field*, void*);
void force0(double, suNg_av_field*, void*);
void force_hmc_auxfields( double dt, suNg_av_field *force, void *par ); //Force from a four_fermion monomial
void force_hmc_ff(double, suNg_av_field*, void*); //Force from a HMC_ff or Hasenbusch_ff monomial

void gaussian_momenta(suNg_av_field *momenta);
void gaussian_spinor_field(spinor_field *s);
void gaussian_spinor_field_flt(spinor_field_flt *s);
void z2_spinor_field(spinor_field *s);


/* For the fermion force ? */
void corret_pf_dist_hmc();
void calc_one_force(int n_force);


/* Action structures */
typedef enum {
  PureGauge,
  FourFermion,
  HMC,
  RHMC,
  TM,
  TM_alt,
  Hasenbusch,
  Hasenbusch_tm,
  Hasenbusch_tm_alt,
  HMC_ff,
  Hasenbusch_ff
} mon_type;

typedef struct _mon_pg_par {
  double beta;
} mon_pg_par;

//Parameters for the four fermion auxiliary field monomial
typedef struct _mon_ff_par {
  double gamma;
  char *start_config;  //cold -> set sigmafield to start_value, pi to 0
                       //Anything else -> use configuration file
  double start_value;
  force_auxfield_par fpar;
} mon_ff_par;

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


typedef struct _monomial_data {
  int id; /* monomial id */
  mon_type type; /* type of monomial */
  void *par; /* parameters */
  double MT_prec; /* metropolis precision */
  double MD_prec; /* molecular dynamics precision */
  double force_prec; /* force precision */
} monomial_data;

typedef struct _monomial {
  monomial_data data;
  
  /* Functions */
  void (*free)(struct _monomial *m); /* free memory */
  
  void (*force_f)(double dt, suNg_av_field *force, void *par); /* force function */
  void *force_par; /* parameters for the force function */

  void (*init_traj)(const struct _monomial *m);
  void (*gaussian_pf)(const struct _monomial *m);
  void (*correct_pf)(const struct _monomial *m);
  void (*correct_la_pf)(const struct _monomial *m);
  const spinor_field *(*pseudofermion)(const struct _monomial *m); /* returns ps field pointer */
  void (*add_local_action)(const struct _monomial *m, scalar_field *loc_action);
  
} monomial;

struct _monomial* pg_create(const monomial_data *data);
struct _monomial* hmc_create(const monomial_data *data);
struct _monomial* rhmc_create(const monomial_data *data);
struct _monomial* tm_create(const monomial_data *data);
struct _monomial* tm_alt_create(const monomial_data *data);
struct _monomial* hasen_create(const monomial_data *data);
struct _monomial* hasen_tm_create(const monomial_data *data);
struct _monomial* hasen_tm_alt_create(const monomial_data *data);

/* Four fermion monomials */
struct _monomial* ff_create(const monomial_data *data);
struct _monomial* hmc_ff_create(const monomial_data *data);
struct _monomial* hasen_ff_create(const monomial_data *data);


const monomial *add_mon(monomial_data *mon);
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
  
  /* integrator */
  integrator_par *integrator;
  double tlen;

  /* Fermion Theta angles */
  double theta[4];
  
  /* Probably not needed anymore */
  /* SF stuff */
  double SF_zf;
  double SF_ds;
  int SF_sign;
  double SF_ct;
  
} ghmc_par;

void init_ghmc(ghmc_par *par);
void free_ghmc();
int update_ghmc();



/* local action */
typedef enum {
   NEW=1,
   DELTA=2
} local_action_type;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
void local_hmc_action(local_action_type type,
                      scalar_field *loc_action,
                      suNg_av_field *momenta);
void pf_local_action(scalar_field *loc_action,
                     spinor_field *pf);


void suNg_field_copy(suNg_field *g1, suNg_field *g2);
void suNf_field_copy(suNf_field *g1, suNf_field *g2);

/* find spectral interval using eva */
void find_spec_H2(double *max, double *min);


/* Utility functions for four fermion interactions */
void scalar_field_copy(scalar_field *s1, scalar_field *s2);
void flip_scalar_field(scalar_field *s);
void set_scalar_field(scalar_field *s, double c);
void gaussian_scalar_field(scalar_field *s);

void update_auxfields( double dt );



#endif
