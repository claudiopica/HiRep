/// Header file for:
/// - update_hb.c
/// - update_hb_multilevel.c
/// - update_mt.c

#ifndef UPDATE_ALGO_H
#define UPDATE_ALGO_H

#include "integrators.h"
#include "Observables/glueballs.h"

#ifdef __cplusplus
	extern "C" {
#endif

//update_hb.c
void project_gauge_field(void);
void update(double *beta, int nhb, int nor);

//update_hb_multilevel.c
void set_max_mh_level(int max_lev);
void update_hb_multilevel_gb_measure(int lev, double *beta, int nhb, int nor, int *ml_up, int *ml_skip, int nblockingstart, int nblockingsend, double *smear_val, cor_list *lcor);

//update_mt.c
typedef struct _ghmc_par {
  /* integrator */
  integrator_par *integrator;
  double tlen;
  double csw;
  double rho_s;
  double rho_t;

  /* Fermion Theta angles */
  double theta[4];

  /* SF stuff */
  double SF_zf;
  double SF_ds;
  int SF_sign;
  double SF_ct;
  int SF_background;

} ghmc_par;

void init_ghmc(ghmc_par *par);
void free_ghmc();
int update_ghmc();
int reverse_update_ghmc();
#ifdef MEASURE_FORCEHMC
//TODO: these two are only used in
// in HMC/hmc_forces.c 
// move to local header and move the implementation there
void corret_pf_dist_hmc();
void calc_one_force(int n_force);
#endif

#ifdef __cplusplus
	}
#endif
#endif //UPDATE_ALGO_H
