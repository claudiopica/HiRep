#ifndef UPDATE_H
#define UPDATE_H

#include "suN.h"
#include "inverters.h"

void staples(int ix,int mu,suNg *v);
void test_staples();
void cabmar(double beta,suNg *u, suNg *v,int type);

void project_gauge_field(void);

void update(double beta,int nhb,int nor);
void random_su2(double rho,double s[]);


void Force0(double dt, suNg_algebra_vector *force);
void Force(double dt, suNg_algebra_vector *force);

typedef struct {
	double tlen; /* trajectory lenght */
	unsigned int nsteps; /* number of step in the integration */
	unsigned int gsteps; /* number of substeps for the gauge part every step */
} int_par;

void leapfrog(suNg_algebra_vector *momenta, int_par *traj_par);
void O2MN_multistep(suNg_algebra_vector *momenta, int_par *traj_par);

void gaussian_momenta(suNg_algebra_vector *momenta);
void gaussian_spinor_field(spinor_field *s);
void gaussian_spinor_field_flt(spinor_field_flt *s);

typedef struct {
  /* sim parameters */
  double beta;
  int nf;
  double mass;
	
	double MT_prec; /* metropolis test precision */
	double MD_prec; /* molecular dynamics precision */
	double HB_prec; /* heatbath precision for pseudofermions */
	double force_prec; /* precision used in the inversions in the force */
	unsigned int n_pf; /* number of psudofermions used in the evolution */
	void (*integrator)(suNg_algebra_vector *, int_par *); /* integrator used in MD */
	int_par *MD_par;
	int (*mshift_solver)(mshift_par *, spinor_operator, spinor_field *, spinor_field *);
} rhmc_par;
void init_rhmc(rhmc_par *par);
void free_rhmc();

/* update the gauge field using RHMC algorithm
 * return code: (<0 means an error has occurred)
 * 0 => conf correctly generated but has not passed metropolis test
 * 1 => conf accepted
 * 
 * -1 => rhmc has not been initialized. call init_rhmc first.
 */
int update_rhmc();

/* this is the basic operator used in the update */
void H2(spinor_field *out, spinor_field *in);
void Force_rhmc_f(double dt, suNg_algebra_vector *force);

typedef enum {
   NEW=1,
   DELTA=2
} local_action_type;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
void local_hmc_action(local_action_type type,
                      double *loc_action,
                      suNg_algebra_vector *momenta,
                      suNf_spinor **phi1,
                      suNf_spinor **phi2);

void suNg_field_copy(suNg *g1, suNg *g2);

/* use power method to find min eigvalue of H2 */
int max_H2(double *min, double mass);

/* find spectral interval using eva */
void find_spec_H2(double *max, double *min, double mass);

#endif
