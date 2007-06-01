#ifndef UPDATE_H
#define UPDATE_H

#include "suN.h"

void staples(int ix,int mu,suNg *v);
void test_staples();
void cabmar(float beta,suNg *u, suNg *v,int type);

void project_gauge_field(void);

void update(float beta,int nhb,int nor);
void random_su2(float rho,float s[]);

typedef struct {
  /* sim parameters */
  float beta;
  int nf;
  float mass;
  /* MD parametes */
  float tlen;
  int nsteps;
} rhmc_par;

void Force0(float dt, suNg_algebra_vector *force);
void Force(float dt, suNg_algebra_vector *force);

void leapfrog(suNg_algebra_vector *momenta, float tlen, unsigned int nsteps);
void O2MN_multistep(suNg_algebra_vector *momenta, float tlen, unsigned int nsteps, unsigned int gsteps);

void gaussian_momenta(suNg_algebra_vector *momenta);
void gaussian_spinor_field(suNf_spinor *s);

void init_rhmc(rhmc_par *par);
void free_rhmc();
int update_rhmc();
/* this is the basic operator used in the update */
void H2(suNf_spinor *out, suNf_spinor *in);
void Force_rhmc_f(float dt, suNg_algebra_vector *force);

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
                      suNf_spinor *phi1,
                      suNf_spinor *phi2);

void suNg_field_copy(suNg *g1, suNg *g2);

/* use power method to find min eigvalue of H2 */
void max_H2(double *min, double mass);
/* use power method to find min eigvalue of H2 */
void min_H2(double *min, double max, double mass);

void find_spec_H2(double *max, double *min, double mass);

#endif
