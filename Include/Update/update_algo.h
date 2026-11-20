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
void exec_project(void); // GPU helper

//update_hb_multilevel.c
void init_hb_multilevel(int lev, double lbeta, int lnhb, int lnor, int *lml_up, int *lml_skip, int lnblockingstart,
                        int lnblockingend, double lsmear_val, cor_list *llcor);
void update_hb_multilevel_gb_measure(int lev);
void update_hb_multilevel_gb_tune(int tuning_level);

//update_mt.c
typedef struct ghmc_par {
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
#ifdef __cplusplus
}
#endif
#endif //UPDATE_ALGO_H
