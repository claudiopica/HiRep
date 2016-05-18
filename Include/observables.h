/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File observables.h
* 
* Functions for measuring observables
*
*******************************************************************************/

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "suN.h"
#include "inverters.h"
#include "meson_observables.h"
#include <stdio.h>

#define SPIN_2D_INDEX(i,j) ( (i)*4 + (j) )

double plaq(int ix,int mu,int nu);
void cplaq(complex *ret,int ix,int mu,int nu);
double avr_plaquette();
double SF_action(double beta);
double local_plaq(int ix);
void full_plaquette();

double rect_1x2(int ix,int mu,int nu);
void crect_1x2(complex *ret,int ix,int mu,int nu);
double avr_rect_1x2();
void full_rect_1x2();
double local_rect_1x2(int ix);

void SF_PCAC_wall_mass(double mass, double acc);
void SF_quark_propagator(spinor_field *in, double mass, spinor_field *out, double acc);

void polyakov();

void pta_qprop_QMR_eo(int g0[4], spinor_field **pta_qprop, int nm, double *m, double acc);
void pta_qprop_QMR(int g0[4], spinor_field **pta_qprop, int nm, double *m, double acc);
void pta_qprop_MINRES(int g0[4], spinor_field **pta_qprop, int nm, double *m, double acc);


typedef enum {NO_DILUTION, TIME_DILUTION, TIME_SPIN_DILUTION, EXACT} dilution_mode;
typedef struct _ata_qprop_pars {
  int n_masses;
  double mass[256];
  int n_eigenvalues;
  int eva_nevt;
  double eva_omega1;
  double eva_omega2;
  int eva_imax;
  int eva_kmax;
  int hopping_order;
  int n_truncation_steps;
  int n_sources_truncation;
  int n_sources_correction;
  dilution_mode dilution;
  double inverter_precision;
} ata_qprop_pars;

void traced_ata_qprop(complex*** prop, int n_points);
void ata_qprop_init(ata_qprop_pars* p);
void ata_qprop_free();

void z2semwall_qprop_free();
void z2semwall_mesons(int conf, int nhits, int nm, double *m, double acc);

void z2semwall_qprop_free_new();
void z2semwall_mesons_new(int conf, int nhits, int nm, double *m, double acc);

void create_point_source(spinor_field *source,int tau,int color);
void create_full_point_source(spinor_field *source, int tau);
void create_point_source_loc(spinor_field *source, int t, int x, int y, int z, int color);
int create_diluted_source_equal_eo(spinor_field *source);
void create_diluted_source_equal_atau_eo(spinor_field *source, int tau);
int create_diluted_source_equal(spinor_field *source);
void create_diluted_source_equal_spinorfield1(spinor_field *source,int tau);
void create_diluted_source_equal_atau(spinor_field *source, int tau);
void create_noise_source_equal_eo(spinor_field *source);
void create_noise_source_equal_oe(spinor_field *source);
void create_diluted_source_equal_atau_col(spinor_field *source, int tau,int col);
void create_noise_source_equal_col_dil(spinor_field *source,int col);
void create_gauge_fixed_wall_source(spinor_field *source, int tau, int color);
void create_gauge_fixed_momentum_source(spinor_field *source, int pt, int px, int py, int pz, int color);
void create_sequential_source(spinor_field *source, int tf, spinor_field* prop);
void restrict_timeslice(spinor_field *source, int tf, spinor_field* prop);
void create_diluted_volume_source(spinor_field *source, int parity_component, int mod);
void add_momentum(spinor_field* out, spinor_field* in, int px, int py, int pz);

void init_propagator_eo(int nm, double *m, double acc);
void eig_init(int nev, int nevt, int kmax, int maxiter, double lbnd, double omega1, double omega2);
void free_propagator_eo();
void calc_propagator_eo(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator(spinor_field *psi, spinor_field* eta, int ndilute);
void calc_propagator_multisource(spinor_field *psi, spinor_field* eta, int ndilute);
void calc_deflated_propagator(spinor_field *psi, spinor_field* eta, int ndilute, int Nuse);
void copy_evec( int n, spinor_field* psi1, double *eval );

void init_meson_correlators(int meas_offdiag);
void init_discon_correlators();
void init_vcvl_correlators();
void init_cvc_correlators();
void free_meson_observables();

void measure_mesons_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, meson_observable* mo, int nm, int tau, int n_mom, int offset,int lt);
void measure_mesons(meson_observable* mo,spinor_field *psi0, spinor_field *eta, int nm,int tau);
void measure_diquarks(meson_observable* mo, spinor_field *psi0, spinor_field *psi1, spinor_field *eta, int nm, int tau);
void measure_conserved_currents(meson_observable* mo,spinor_field* psi0, spinor_field* eta, int nm, int tau);
void measure_mesons_ext(meson_observable* mo,spinor_field *psi0, spinor_field *eta, int nm,int tau,int begin);
void measure_point_mesons_momenta(meson_observable* mo,spinor_field* psi0, spinor_field* eta, int nm, int tau, int n_mom);
void measure_point_mesons_momenta_ext(meson_observable* mo,spinor_field* psi0, spinor_field* eta, int nm, int tau, int n_mom, int begin);
void measure_formfactors(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int ti, int tf, int n_mom, int *pt);
void measure_formfactors_ext(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int ti, int tf, int n_mom, int begin);
void print_mesons(meson_observable* mo,double norm, int conf, int nm, double* mass, int lt, int n_mom, char* label);
void print_formfactor(int conf, int nm, double* mass, int n_mom, char* label, int tf);
void print_formfactor_ext(int conf, int nm, double* mass, int n_mom, char* label, int tf);

void measure_scattering_AD_core(meson_observable* mo, spinor_field* psi0,spinor_field* psi1,spinor_field* psi2,spinor_field* psi3, int tau, int split, int n_mom, int p_tot_x, int p_tot_y, int p_tot_z);
void measure_scattering_BC_core(meson_observable* mo, spinor_field* psi0,spinor_field* psi1, spinor_field* psi2,spinor_field* psi3, int tau, int split, int n_mom, int p_tot_x, int p_tot_y, int p_tot_z);

void measure_renormalization(spinor_field* psi_in, spinor_field* psi_out, int nm, int pt_in, int px_in, int py_in, int pz_in, int pt_out, int px_out, int py_out, int pz_out);
void print_renormalization(int conf, int nm, double* mass, char* label, int pt_in, int px_in, int py_in, int pz_in, int pt_out, int px_out, int py_out, int pz_out);

void contract_baryons(spinor_field *psi0,int tau);
void measure_glueballs();

void id_correlator(double *out, int t0, spinor_field *qp);
void g0_correlator(double *out, int t0, spinor_field *qp);
void g5_correlator(double *out, int t0, spinor_field *qp);
void g0g5_correlator(double *out, int t0, spinor_field *qp);
void g1_correlator(double *out, int t0, spinor_field *qp);
void g2_correlator(double *out, int t0, spinor_field *qp);
void g3_correlator(double *out, int t0, spinor_field *qp);
void g0g1_correlator(double *out, int t0, spinor_field *qp);
void g0g2_correlator(double *out, int t0, spinor_field *qp);
void g0g3_correlator(double *out, int t0, spinor_field *qp);
void g5g1_correlator(double *out, int t0, spinor_field *qp);
void g5g2_correlator(double *out, int t0, spinor_field *qp);
void g5g3_correlator(double *out, int t0, spinor_field *qp);
void g0g5g1_correlator(double *out, int t0, spinor_field *qp);
void g0g5g2_correlator(double *out, int t0, spinor_field *qp);
void g0g5g3_correlator(double *out, int t0, spinor_field *qp);
void g5_g0g5_re_correlator(double *out, int t0, spinor_field *qp);
void g5_g0g5_im_correlator(double *out, int t0, spinor_field *qp);

void id_trace_H(complex* out, complex* smat);
void g0_trace_H(complex* out, complex* smat);
void g5_trace_H(complex* out, complex* smat);
void g0g5_trace_H(complex* out, complex* smat);
void g1_trace_H(complex* out, complex* smat);
void g2_trace_H(complex* out, complex* smat);
void g3_trace_H(complex* out, complex* smat);
void g0g1_trace_H(complex* out, complex* smat);
void g0g2_trace_H(complex* out, complex* smat);
void g0g3_trace_H(complex* out, complex* smat);
void g5g1_trace_H(complex* out, complex* smat);
void g5g2_trace_H(complex* out, complex* smat);
void g5g3_trace_H(complex* out, complex* smat);
void g0g5g1_trace_H(complex* out, complex* smat);
void g0g5g2_trace_H(complex* out, complex* smat);
void g0g5g3_trace_H(complex* out, complex* smat);

void id_debug(complex Gamma[4][4], int* sign);
void g0_debug(complex Gamma[4][4], int* sign);
void g5_debug(complex Gamma[4][4], int* sign);
void g0g5_debug(complex Gamma[4][4], int* sign);
void g1_debug(complex Gamma[4][4], int* sign);
void g2_debug(complex Gamma[4][4], int* sign);
void g3_debug(complex Gamma[4][4], int* sign);
void g0g1_debug(complex Gamma[4][4], int* sign);
void g0g2_debug(complex Gamma[4][4], int* sign);
void g0g3_debug(complex Gamma[4][4], int* sign);
void g5g1_debug(complex Gamma[4][4], int* sign);
void g5g2_debug(complex Gamma[4][4], int* sign);
void g5g3_debug(complex Gamma[4][4], int* sign);
void g0g5g1_debug(complex Gamma[4][4], int* sign);
void g0g5g2_debug(complex Gamma[4][4], int* sign);
void g0g5g3_debug(complex Gamma[4][4], int* sign);


void id_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g0_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g1_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g2_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g3_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g5_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g0g5_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g5g1_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g5g2_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g5g3_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g0g1_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g0g2_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g0g3_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g0g5g1_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g0g5g2_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);
void g0g5g3_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in);


#define _WL_3VOL_INDEX(x,y,z) ((x)+(y)*X+(z)*X*Y)
void WL_initialize();
void WL_free();
void WL_load_path(int c[3], int nsteps);
void WL_Hamiltonian_gauge(suNg_field* out, suNg_field* in);
void WL_broadcast_polyakov(suNg* poly, suNg_field* gf);
void WL_correlators(double** ret, const suNg_field* gf, const suNg* poly, const int nsteps, const int* path, const int length, const int perm[3], int sign[3]);
void WL_wilsonloops(double HYP_weight[3]);


void init_modenumber(double m, double inv, int nh, char *approxfile);
void free_modenumber();
double ModeNumber(double M2);


typedef struct {
  complex ***g1_ij, ***g2_ij, ***g3_ij, ***g4_ij, *g1, *g2, *g3, *g4, **M;

  complex *l11, *l12, *l13;
  complex *l21, *l22, *l23;
  complex *l31, *l32, *l33;
  complex *l41, *l42, *l43;

  complex ***l11_ij, ***l12_ij, ***l13_ij;
  complex ***l21_ij, ***l22_ij, ***l23_ij;
  complex ***l31_ij, ***l32_ij, ***l33_ij;
} chisf_mem;

chisf_mem *  init_rotated_corr_mem();

void rotated_gXuup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gXddp(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gXudp(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gXdup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gvtuup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gvtddp(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gvtdup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gvtudp(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_g1uup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_g1ddp(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_g1udp(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_g1dup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_lXuup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_lXddp(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_lXdup(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_lXudp(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gXuum(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gXddm(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gXudm(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gXdum(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gvtuum(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gvtddm(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gvtdum(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_gvtudm(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_g1uum(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_g1ddm(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_g1udm(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_g1dum(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_lXuum(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_lXddm(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_lXudm(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);
void rotated_lXdum(chisf_mem* corr_mem, suNf_spinor *chi, spinor_field * prop_uu,spinor_field * prop_dd);




/* Functions that include four fermion interactions */
void init_propagator_ff_eo(int nm, double *m, double acc);
void free_propagator_ff_eo();

void calc_propagator_ff(spinor_field *psi, spinor_field* eta, int ndilute);
void calc_propagator_ff_eo(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator_ff_oe(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator_ff_hopping_eo(spinor_field *psi, spinor_field *eta, int hopping, int ndilute);
void calc_propagator_ff_hopping_oe(spinor_field *psi, spinor_field *eta, int hopping, int ndilute);
void calc_propagator_ff_hopping_series(spinor_field *psi, spinor_field *eta, int hopping, int ndilute);

/* Disconnected part of triplet correlators */
void init_triplet_discon_correlators();

void ff_observables();


#endif 

