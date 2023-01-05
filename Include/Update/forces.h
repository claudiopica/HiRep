/***************************************************************************\
* Copyright (c) 2022, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/// Headerfile for:
/// - force0.c
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
typedef struct {
  double beta;
  double c0;
  double c1;
  suNg_av_field **momenta;
} force_gauge_par;

void force0(double, void *);

//luscherweisz.c 
double lw_action(double beta, double c0, double c1);
void lw_local_action(scalar_field *loc_action, double beta, double c0, double c1);
void lw_force(double dt, void *vpar);
void calculate_stfld(int comm);
double lw_action_density(int ix, double beta, double c0, double c1);

//fermion_force_core.c
#ifdef WITH_CLOVER
void force_clover_logdet(double mass,double residue); //TODO: this simply forwards to compute_force_logdet. can we remove it?
#endif
#ifdef WITH_EXPCLOVER
void force_clover_fermion(spinor_field *Xs, spinor_field *Ys, double residue);
void force_clover_fermion_taylor(spinor_field *Xs, spinor_field *Ys, double residue);
#endif
void force_fermion_core(spinor_field *Xs, spinor_field *Ys, int auto_fill_odd, double dt, double residue);
void fermion_force_begin(void);
void fermion_force_end(double dt, suNg_av_field *force);

//force_hmc.c
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

void free_force_hmc(void); //TODO: this should be static but are used in force_hmc_ff.c
void init_force_hmc(void); //TODO: this should be static but are used in force_hmc_ff.c
void force_hmc(double, void *);

//force_hmc_tm.c
void force_hmc_tm(double, void *); //uses force_hmc_par

//force_rhmc.c
typedef struct {
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
typedef struct {
  double gamma;
} force_auxfield_par;

void force_hmc_auxfields(double dt, void *vpar); //Force from a four_fermion monomial
void force_hmc_auxfields_fermion(double dt, void *vpar, scalar_field *sigma_mom, scalar_field *pi_mom, spinor_field *Xs, spinor_field *Ys, int hasenbusch);
void update_auxfields(double dt, void *vpar);

//force_scalar.c
typedef struct {
  double mass;
  double lambda;
  suNg_scalar_field **momenta;
  suNg_av_field **g_momenta;
} force_scalar_par;

void force_scalar(double dt,void *par);

#ifdef __cplusplus
	}
#endif
#endif //FORCE0_H
