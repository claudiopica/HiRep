/***************************************************************************\
* Copyright (c) 2022, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/// Headerfile for:
/// - force0.c
/// - force0_gpu.cu
/// - luscherweisz.c
/// - force_fermion_core.c
/// - force_hmc.c
/// - force_hmc_tm.c
/// - force_rhmc.c
/// - force_hmc_ff.c
/// - force_4f.c
/// - force_scalar.c

#ifndef FORCE0_H
#define FORCE0_H

#include "spinor_field.h"
#include "mre.h"
#include "rational_functions.h"
#include "field_update.h"

#ifdef __cplusplus
extern "C" {
#endif

//force0.c
typedef struct force_gauge_par {
    double beta;
    double c0;
    double c1;
    suNg_av_field **momenta;
} force_gauge_par;

void force0(double, void *);

//luscherweisz.c
void lw_force_gpu(double dt, void *vpar);
void lw_force_cpu(double dt, void *vpar);

double lw_action_gpu(double, double, double);
double lw_action_cpu(double, double, double);

void lw_local_action_gpu(scalar_field *, double, double, double);
void lw_local_action_cpu(scalar_field *, double, double, double);

extern double (*lw_action)(double beta, double c0, double c1);
extern void (*lw_local_action)(scalar_field *loc_action, double beta, double c0, double c1);
extern void (*lw_force)(double dt, void *vpar);
extern void (*calculate_stfld)(int comm);
double lw_action_density(int ix, double beta, double c0, double c1);

//fermion_force_core.c
#ifdef WITH_CLOVER
void force_clover_logdet(double mass, double residue); //TODO: this simply forwards to compute_force_logdet. can we remove it?
#endif
#ifdef WITH_EXPCLOVER
extern void (*force_clover_fermion)(spinor_field *Xs, spinor_field *Ys, double residue);
void force_clover_fermion_taylor(spinor_field *Xs, spinor_field *Ys, double residue);
#endif

void force_fermion_core_gpu(spinor_field *, spinor_field *, int, double, double);
void force_fermion_core_cpu(spinor_field *, spinor_field *, int, double, double);

void fermion_force_begin_gpu(void);
void fermion_force_begin_cpu(void);

void fermion_force_end_gpu(double, suNg_av_field *);
void fermion_force_end_cpu(double, suNg_av_field *);

extern void (*force_fermion_core)(spinor_field *Xs, spinor_field *Ys, int auto_fill_odd, double dt, double residue);
extern void (*fermion_force_begin)(void);
extern void (*fermion_force_end)(double dt, suNg_av_field *force);
#ifdef WITH_GPU
void call_fermion_kernel(spinor_field *Xs, spinor_field *Ys, suNg_av_field *force_sum, double coeff);
void exec_calculate_stfld(suNg_field *stfld[], int comm);
void exec_lw_force(suNg_field **stfld, suNg_av_field *force, double dt, double beta, double c0, double c1);
#endif

//force_hmc.c
typedef struct force_hmc_par {
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

void free_force_hmc(void); //TODO: this should be static but are used in force_hmc_ff.c
void init_force_hmc(void); //TODO: this should be static but are used in force_hmc_ff.c
void force_hmc(double, void *);

//force_hmc_tm.c
void force_hmc_tm(double, void *); //uses force_hmc_par

//force_rhmc.c
typedef struct force_rhmc_par {
    int id;
    int n_pf;
    spinor_field *pf;
    double mass;
    rational_app *ratio;
    double inv_err2;
    suNg_av_field **momenta;
} force_rhmc_par;

void force_rhmc(double, void *);

//force_hmc_ff.c
void force_hmc_ff(double dt, void *vpar); //Force from a HMC_ff or Hasenbusch_ff monomial. Uses force_hmc_par

//force_4f.c
typedef struct force_auxfield_par {
    double gamma;
} force_auxfield_par;

void force_hmc_auxfields(double dt, void *vpar); //Force from a four_fermion monomial
void force_hmc_auxfields_fermion(double dt, void *vpar, scalar_field *sigma_mom, scalar_field *pi_mom, spinor_field *Xs,
                                 spinor_field *Ys, int hasenbusch);
void update_auxfields(double dt, void *vpar);

//force_scalar.c
typedef struct force_scalar_par {
    double mass;
    double lambda;
    suNg_scalar_field **momenta;
    suNg_av_field **g_momenta;
} force_scalar_par;

void force_scalar(double dt, void *par);

#ifdef WITH_GPU
void force0_kernel_gpu(suNg_av_field *force, double coeff);
#endif

#ifdef __cplusplus
}
#endif
#endif //FORCE0_H
