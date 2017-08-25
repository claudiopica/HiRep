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

typedef struct {
	int id;
	int n_pf;
	spinor_field *pf;
	double mass;
	rational_app *ratio;
	double inv_err2;
	suNg_av_field **momenta;
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
	int logdet;
	suNg_av_field **momenta;
} force_hmc_par;

typedef struct {
	double beta;
	double c0;
	double c1;
	suNg_av_field **momenta;
} force_gauge_par;

typedef struct {
	double mass;
	double lambda;
	suNg_scalar_field **momenta;
	suNg_av_field **g_momenta;
} force_scalar_par;

typedef struct {
	double gamma;
} force_auxfield_par;

typedef struct {
	suNg_field **field;
	suNg_av_field **momenta;
} field_gauge_par;

typedef struct {
	suNg_scalar_field **field;
	suNg_scalar_field **momenta;
} field_scalar_par;

void update_gauge_field(double, void*);
void update_auxfields(double, void*);

void update_scalar_field(double, void*);
void force_scalar(double, void*);

void lw_force(double, void*);
void lw_local_action(scalar_field*, double, double, double);

void fermion_force_begin();
void fermion_force_end(double dt, suNg_av_field*);
void force_fermion_core(spinor_field*, spinor_field*, int, double, double);
void force_clover_logdet(double, double);

void force_hmc(double, void*);
void force_hmc_tm(double, void*);
void force_rhmc(double, void*);
void force0(double, void*);
void force_hmc_auxfields(double, void*); //Force from a four_fermion monomial
void force_hmc_ff(double, void*); //Force from a HMC_ff or Hasenbusch_ff monomial

void gaussian_momenta(suNg_av_field *momenta);
void gaussian_scalar_momenta(suNg_scalar_field *momenta);
void gaussian_spinor_field(spinor_field *s);
void gaussian_spinor_field_flt(spinor_field_flt *s);
void z2_spinor_field(spinor_field *s);

/* For the fermion force ? */
void corret_pf_dist_hmc();
void calc_one_force(int n_force);

#include "monomials.h"

typedef struct _integrator_par {
  int nsteps;
  int nmon;
  const monomial **mon_list;
  void (*integrator)(double, struct _integrator_par*);
  struct _integrator_par *next;
  int level;
} integrator_par;

void leapfrog_multistep(double tlen, integrator_par *int_par);
void O2MN_multistep(double tlen, integrator_par *int_par);
void O4MN_multistep(double tlen, integrator_par *int_par);


typedef struct _ghmc_par {
  
  /* integrator */
  integrator_par *integrator;
  double tlen;
  double csw;
  double rho_s;
  double rho_t;

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

/* stout smearing */
void init_smearing(double, double);
double avr_smeared_plaquette();
void smear_gauge_field();
void smeared_gauge_force(suNg_av_field*,suNg_av_field*);

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
                      suNg_av_field *momenta,
                      suNg_scalar_field *momenta_s);
void pf_local_action(scalar_field *loc_action,
                     spinor_field *pf);


void suNg_field_copy(suNg_field *g1, suNg_field *g2);
void suNf_field_copy(suNf_field *g1, suNf_field *g2);
void suNg_scalar_field_copy(suNg_scalar_field *g1, suNg_scalar_field *g2);

/* find spectral interval using eva */
void find_spec_H2(double *max, double *min);


/* Utility functions for four fermion interactions */
void scalar_field_copy(scalar_field *s1, scalar_field *s2);
void flip_scalar_field(scalar_field *s);
void set_scalar_field(scalar_field *s, double c);
void gaussian_scalar_field(scalar_field *s);



#endif
