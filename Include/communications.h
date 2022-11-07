/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef COMMUNICATIONS_H
#define COMMUNICATIONS_H

#ifdef __cplusplus
    extern "C" {
#endif

void global_sum(double *d, int n);
void global_sum_int(int *d, int n);
void global_max(double *d, int n);
void global_min(double *d, int n);
void bcast(double *d, int n);
void bcast_int(int *i, int n);

#include "spinor_field.h"
void complete_gf_sendrecv(suNg_field *gf);
void start_gf_sendrecv(suNg_field *gf);
void complete_sf_sendrecv(spinor_field *gf);
void start_sf_sendrecv(spinor_field *gf);
void complete_sc_sendrecv(suNg_scalar_field *gf);
void start_sc_sendrecv(suNg_scalar_field *gf);

void complete_clover_force_sendrecv(suNf_field *gf);
void start_clover_force_sendrecv(suNf_field *gf);

void complete_gt_sendrecv(suNg_field *gf);
void start_gt_sendrecv(suNg_field *gf);

void test_spinor_field(spinor_field *p);

/* Floating point sendrecv */
void complete_gf_sendrecv_flt(suNg_field_flt *gf);
void start_gf_sendrecv_flt(suNg_field_flt *gf);
void complete_sf_sendrecv_flt(spinor_field_flt *gf);
void start_sf_sendrecv_flt(spinor_field_flt *gf);


#ifdef WITH_GPU
    int global_sum_gpu_int(int *vector, int size);
    float global_sum_gpu_float(float *vector, int size);
    double global_sum_gpu_double(double *vector, int size);
    hr_complex_flt global_sum_gpu_complex_flt(hr_complex_flt *vector, int size);
    hr_complex global_sum_gpu_complex(hr_complex *vector, int size);
#endif

#ifdef __cplusplus
}
#endif

#ifdef WITH_GPU
    #define _DECLARE_COMMS(_name, _field_type) \
        void sync_gpu_##_name(_field_type*); \
        void start_sendrecv_gpu_##_name(_field_type*); \
        void complete_sendrecv_gpu_##_name(_field_type*); 

    _DECLARE_COMMS(spinor_field_f, spinor_field);
    _DECLARE_COMMS(spinor_field_f_flt, spinor_field_flt);
    _DECLARE_COMMS(gfield, suNg_field);
    _DECLARE_COMMS(gfield_flt, suNg_field_flt);
    _DECLARE_COMMS(gfield_f, suNf_field);
    _DECLARE_COMMS(gfield_f_flt, suNg_field);
    _DECLARE_COMMS(scalar_field, suNg_scalar_field);
    _DECLARE_COMMS(avfield, suNg_av_field);
    _DECLARE_COMMS(gtransf, suNg_field);
    _DECLARE_COMMS(clover_term, suNfc_field);
    _DECLARE_COMMS(clover_force, suNf_field);

    #undef _DECLARE_COMMS
#endif 

#endif /* COMMUNICATIONS_H */
