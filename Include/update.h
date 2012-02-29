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
void Force0(double dt, suNg_av_field *force, void *par);

typedef struct {
  spinor_field *pf;
  rational_app *ratio;
} force_rhmc_par;
void Force_rhmc_f(double dt, suNg_av_field *force, void *par);

typedef struct {
  spinor_field *pf;
  int hasenbusch; /* 0 = no hasenbusch ; 1 = force with Dtilde = a*D+b ; 2 = force with Y = Ddag^{-1}(a*phi+b*X) */
  double aD, bD, aY, bY;
} force_hmc_par;
void Force_hmc_f(double dt, suNg_av_field *force, void *par);



typedef struct _integrator_par {
  double tlen;
  int nsteps;
  void (*force)(double,suNg_av_field*,void*);
  void *force_par;
  void (*integrator)(suNg_av_field*, struct _integrator_par*);
  struct _integrator_par *next;
  int level;
} integrator_par;
void gauge_integrator(suNg_av_field *momenta, integrator_par *int_par);
void O2MN_multistep(suNg_av_field *momenta, integrator_par *int_par);



void gaussian_momenta(suNg_av_field *momenta);
void gaussian_spinor_field(spinor_field *s);
void gaussian_spinor_field_flt(spinor_field_flt *s);
void z2_spinor_field(spinor_field *s);


typedef struct _rhmc_par {
  /* sim parameters */
  double beta;
  int nf;
  double mass;
  
	double SF_zf;
	double SF_ds;
	int SF_sign;
	double SF_ct;
	
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
void init_hmc(rhmc_par *par);
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
void local_hmc_action(local_action_type type,
                      scalar_field *loc_action,
                      suNg_av_field *momenta,
                      spinor_field *phi1,
                      spinor_field *phi2);

void suNg_field_copy(suNg_field *g1, suNg_field *g2);

/* find spectral interval using eva */
void find_spec_H2(double *max, double *min, double mass);

#endif
