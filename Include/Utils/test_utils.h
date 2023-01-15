/***************************************************************************\
* Copyright (c) 2022 Sofie Martins                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file test_utils.h
 * @brief utilities to simplify testing
 */

#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "spinor_field.h"

#ifdef __cplusplus
    extern "C" {
#endif

// Testing utils
void test_setup();
int check_diff_norm(double, double);
int check_diff_norm_zero(double);
int check_finiteness(double);

// COPY 
void copy_gfield_cpu(suNg_field*, suNg_field*);
void copy_suNg_scalar_field_cpu(suNg_scalar_field*, suNg_scalar_field*);
void copy_gfield_flt_cpu(suNg_field_flt*, suNg_field_flt*);
void copy_gfield_f_cpu(suNf_field*, suNf_field*);
void copy_gfield_f_flt_cpu(suNf_field_flt*, suNf_field_flt*);
void copy_avfield_cpu(suNg_av_field*, suNg_av_field*);
void copy_sfield_cpu(scalar_field*, scalar_field*);
void copy_clover_ldl_cpu(ldl_field*, ldl_field*);
void copy_gtransf_cpu(suNg_field*, suNg_field*);
void copy_clover_term_cpu(suNfc_field*, suNfc_field*);
void copy_clover_force_cpu(suNf_field*, suNf_field*);

// RANDOM
void random_spinor_field_f_cpu(spinor_field*);
void random_spinor_field_f_flt_cpu(spinor_field_flt*);
void random_gfield_cpu(suNg_field*);
void random_suNg_scalar_field_cpu(suNg_scalar_field*);
void random_gfield_flt_cpu(suNg_field_flt*);
void random_gfield_f_cpu(suNf_field*);
void random_gfield_f_flt_cpu(suNf_field_flt*);
void random_avfield_cpu(suNg_av_field*);
void random_sfield_cpu(scalar_field*);
void random_gtransf_cpu(suNg_field*);
void random_clover_ldl_cpu(ldl_field*);
void random_clover_term_cpu(suNfc_field*);
void random_clover_force_cpu(suNf_field*);

// SUB ASSIGN
void sub_assign_gfield_cpu(suNg_field*, suNg_field*);
void sub_assign_gfield_flt_cpu(suNg_field_flt*, suNg_field_flt*);
void sub_assign_gfield_f_cpu(suNf_field*, suNf_field*);
void sub_assign_gfield_f_flt_cpu(suNf_field_flt*, suNf_field_flt*);
void sub_assign_suNg_scalar_field_cpu(suNg_scalar_field*, suNg_scalar_field*);
void sub_assign_avfield_cpu(suNg_av_field*, suNg_av_field*);
void sub_assign_sfield_cpu(scalar_field*, scalar_field*);
void sub_assign_gtransf_cpu(suNg_field*, suNg_field*);
void sub_assign_clover_ldl_cpu(ldl_field*, ldl_field*);
void sub_assign_clover_term_cpu(suNfc_field*, suNfc_field*);
void sub_assign_clover_force_cpu(suNf_field*, suNf_field*);

// SQNORM
double sqnorm_spinor_field_f_cpu(spinor_field*);
float sqnorm_spinor_field_f_flt_cpu(spinor_field_flt*);
double sqnorm_gfield_cpu(suNg_field*);
double sqnorm_gfield_f_cpu(suNf_field*);
float sqnorm_gfield_flt_cpu(suNg_field_flt*);
float sqnorm_gfield_f_flt_cpu(suNf_field_flt*);
double sqnorm_suNg_scalar_field_cpu(suNg_scalar_field*);
double sqnorm_avfield_cpu(suNg_av_field*);
double sqnorm_gtransf_cpu(suNg_field*);
double sqnorm_sfield_cpu(scalar_field*);
double sqnorm_clover_ldl_cpu(ldl_field*);
double sqnorm_clover_term_cpu(suNfc_field*);
double sqnorm_clover_force_cpu(suNf_field*);

// SET ZERO
void zero_gfield_cpu(suNg_field*);
void zero_gfield_f_cpu(suNf_field*);
void zero_gfield_flt_cpu(suNg_field_flt*);
void zero_gfield_f_flt_cpu(suNf_field_flt*);
void zero_suNg_scalar_field_cpu(suNg_scalar_field*);
void zero_avfield_cpu(suNg_av_field*);
void zero_gtransf_cpu(suNg_field*);
void zero_clover_ldl_cpu(ldl_field*);
void zero_clover_term_cpu(suNfc_field*);
void zero_clover_force_cpu(suNf_field*);

// DOUBLE PRECISION TO SINGLE PRECISION
void sync_single_precision_gauge_field();

#ifdef __cplusplus
    }
#endif
#endif
