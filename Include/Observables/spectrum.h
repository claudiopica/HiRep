/***************************************************************************\
* Copyright (c) 2013, Rudy Arthur, Ari Hietanen                             *
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File spectrum.h
* 
* Functions for measuring spectrum
*
*******************************************************************************/

// Header file for:
// - meson_measurements.c
// - baryon_measurements.c
// - meson_measurements_ff.c
// - mesons.c

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "Utils/data_storage.h"
#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

//used in mesons.c and trunc_hairpin.c
#define SPIN_2D_INDEX(i, j) ((i) * 4 + (j))

//meson_measurements.c
void measure_spectrum_semwall(int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                              data_storage_array **ret);
void measure_spectrum_discon_semwall(int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                                     data_storage_array **ret);
void measure_spectrum_discon_gfwall(int nm, double *m, int conf_num, double precision, storage_switch swc,
                                    data_storage_array **ret);
void measure_spectrum_discon_volume(int nm, double *m, int conf_num, double precision, int dil, storage_switch swc,
                                    data_storage_array **ret);
void measure_spectrum_gfwall(int nm, double *m, int conf_num, double precision, storage_switch swc, data_storage_array **ret);
void measure_spectrum_pt(int tau, int nm, double *m, int n_mom, int conf_num, double precision, storage_switch swc,
                         data_storage_array **ret);
void measure_spectrum_semwall_ext(int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                                  data_storage_array **ret);
void measure_spectrum_pt_ext(int tau, int nm, double *m, int n_mom, int conf_num, double precision, storage_switch swc,
                             data_storage_array **ret);
void measure_spectrum_semwall_fixedbc(int dt, int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                                      data_storage_array **ret);
void measure_spectrum_pt_fixedbc(int tau, int dt, int nm, double *m, int n_mom, int conf_num, double precision,
                                 storage_switch swc, data_storage_array **ret);
void measure_spectrum_gfwall_fixedbc(int dt, int nm, double *m, int conf_num, double precision, storage_switch swc,
                                     data_storage_array **ret);
void measure_formfactor_pt(int ti, int tf, int nm, double *m, int n_mom, int conf_num, double precision, storage_switch swc,
                           data_storage_array **ret);
void measure_formfactor_fixed(int ti, int tf, int dt, int nm, double *m, int n_mom, int conf_num, double precision,
                              storage_switch swc, data_storage_array **ret);
void measure_diquark_semwall_background(int nm, double *m, int nhits, int conf_num, double precision, double Q, int n,
                                        storage_switch swc, data_storage_array **ret);
void measure_conserved_formfactor_fixed(int ti, int tf, int dt, int nm, double *m, int n_mom, int conf_num, double precision,
                                        storage_switch swc, data_storage_array **ret); //TODO: not defined in lib

//baryon_measurements.c
void measure_baryons(double *m, int conf_num, double precision, storage_switch swc, data_storage_array **ret);

/* For measuring spectrum with a four fermion interaction */
//meson_measurements_ff.c
void measure_spectrum_ff_semwall(int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                                 data_storage_array **ret);
void measure_spectrum_discon_ff_semwall(int nm, double *m, int nhits, int degree_hopping, int nhits_hopping, int conf_num,
                                        double precision, storage_switch swc, data_storage_array **ret);
void measure_spectrum_ff_pt(int tau, int nm, double *m, int n_mom, int conf_num, double precision, storage_switch swc,
                            data_storage_array **ret);
void measure_spectrum_semwall_ff_ext(int nm, double *m, int nhits, int conf_num, double precision, storage_switch swc,
                                     data_storage_array **ret);

//mesons.c
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

void id_trace_H(hr_complex *out, hr_complex *smat);
void g0_trace_H(hr_complex *out, hr_complex *smat);
void g5_trace_H(hr_complex *out, hr_complex *smat);
void g0g5_trace_H(hr_complex *out, hr_complex *smat);
void g1_trace_H(hr_complex *out, hr_complex *smat);
void g2_trace_H(hr_complex *out, hr_complex *smat);
void g3_trace_H(hr_complex *out, hr_complex *smat);
void g0g1_trace_H(hr_complex *out, hr_complex *smat);
void g0g2_trace_H(hr_complex *out, hr_complex *smat);
void g0g3_trace_H(hr_complex *out, hr_complex *smat);
void g5g1_trace_H(hr_complex *out, hr_complex *smat);
void g5g2_trace_H(hr_complex *out, hr_complex *smat);
void g5g3_trace_H(hr_complex *out, hr_complex *smat);
void g0g5g1_trace_H(hr_complex *out, hr_complex *smat);
void g0g5g2_trace_H(hr_complex *out, hr_complex *smat);
void g0g5g3_trace_H(hr_complex *out, hr_complex *smat);

void id_debug(hr_complex Gamma[4][4], int *sign);
void g0_debug(hr_complex Gamma[4][4], int *sign);
void g5_debug(hr_complex Gamma[4][4], int *sign);
void g0g5_debug(hr_complex Gamma[4][4], int *sign);
void g1_debug(hr_complex Gamma[4][4], int *sign);
void g2_debug(hr_complex Gamma[4][4], int *sign);
void g3_debug(hr_complex Gamma[4][4], int *sign);
void g0g1_debug(hr_complex Gamma[4][4], int *sign);
void g0g2_debug(hr_complex Gamma[4][4], int *sign);
void g0g3_debug(hr_complex Gamma[4][4], int *sign);
void g5g1_debug(hr_complex Gamma[4][4], int *sign);
void g5g2_debug(hr_complex Gamma[4][4], int *sign);
void g5g3_debug(hr_complex Gamma[4][4], int *sign);
void g0g5g1_debug(hr_complex Gamma[4][4], int *sign);
void g0g5g2_debug(hr_complex Gamma[4][4], int *sign);
void g0g5g3_debug(hr_complex Gamma[4][4], int *sign);

void id_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g0_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g1_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g2_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g3_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g5_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g0g5_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g5g1_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g5g2_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g5g3_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g0g1_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g0g2_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g0g3_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g0g5g1_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g0g5g2_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);
void g0g5g3_eval_g5GammaDag_times_spinor(suNf_spinor *out, suNf_spinor *in);

#ifdef __cplusplus
}
#endif
#endif
