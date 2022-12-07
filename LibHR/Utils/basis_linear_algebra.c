#include "communications.h"
#include "hr_complex.h"
#include "geometry_check.h"
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "hr_complex.h"
#include <string.h>
#include "error.h"
#include "ranlux.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// The following functions are primarily for testing purposes
// This is all for CPU
// TODO: Turn this into macros for easier extensibility

void rand_field_dbl(double* d, int n) {
    for (int i = 0; i < n; ++i) 
    {
        d[i] = (double)rand()/(double)RAND_MAX;
    }
}

void rand_field_flt(float* f, int n) {
    for (int i = 0; i < n; i++) 
    {
        f[i] = (float)rand()/(float)RAND_MAX;
    }
}

void copy_gfield_cpu(suNg_field* out, suNg_field* in) 
{
    memcpy(out->ptr, in->ptr, 4*out->type->gsize_gauge*sizeof(suNg));
}

void copy_scalar_field_cpu(suNg_scalar_field *out, suNg_scalar_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg_vector));
}

void copy_gfield_flt_cpu(suNg_field_flt *out, suNg_field_flt *in) 
{
    memcpy(out->ptr, in->ptr, 4*out->type->gsize_gauge*sizeof(suNg_flt));
}

void copy_gfield_f_cpu(suNf_field* out, suNf_field* in) 
{
    memcpy(out->ptr, in->ptr, 4*out->type->gsize_gauge*sizeof(suNf));
}

void copy_gfield_f_flt_cpu(suNf_field_flt *out, suNf_field_flt *in) 
{
    memcpy(out->ptr, in->ptr, 4*out->type->gsize_gauge*sizeof(suNf_flt));
}

void copy_avfield_cpu(suNg_av_field *out, suNg_av_field *in) 
{
    memcpy(out->ptr, in->ptr, 4*out->type->gsize_gauge*sizeof(suNg_algebra_vector));
}

void copy_sfield_cpu(scalar_field *out, scalar_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_spinor*sizeof(double));
}

void copy_gtransf_cpu(suNg_field *out, suNg_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_spinor*sizeof(suNg));
}

void copy_clover_term_cpu(suNfc_field *out, suNfc_field *in) 
{
    memcpy(out->ptr, in->ptr, 4*out->type->gsize_gauge*sizeof(suNfc));
}

void copy_clover_force_cpu(suNf_field *out, suNf_field *in) 
{
    memcpy(out->ptr, in->ptr, 6*out->type->gsize_gauge*sizeof(suNf));
}


// ** SUB ASSIGN **
void sub_assign_gfield_cpu(suNg_field *out, suNg_field *in) 
{
    suNg *site_out, *site_in;
    _MASTER_FOR(in->type, ix)
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNg_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_gfield_flt_cpu(suNg_field_flt *out, suNg_field_flt *in) 
{
    suNg_flt *site_out, *site_in;
    _MASTER_FOR(in->type, ix)
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNg_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_gfield_f_cpu(suNf_field *out, suNf_field *in) 
{
    suNf *site_out, *site_in;
    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }   
    }
}

void sub_assign_gfield_f_flt_cpu(suNf_field_flt *out, suNf_field_flt *in) 
{
    suNf_flt *site_out, *site_in;
    _MASTER_FOR(in->type, ix)
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_scalar_field_cpu(suNg_scalar_field *out, suNg_scalar_field *in) 
{
    suNg_vector *site_out, *site_in;
    _MASTER_FOR(in->type, ix) 
    {
        site_out = _FIELD_AT(out, ix);
        site_in = _FIELD_AT(in, ix);
        _vector_sub_assign_g((*site_out), (*site_in));
    }
}

void sub_assign_avfield_cpu(suNg_av_field *out, suNg_av_field *in) 
{
    suNg_algebra_vector *site_in, *site_out;
    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _algebra_vector_sub_assign_g((*site_out), (*site_in));
        } 
    }
}

void sub_assign_sfield_cpu(scalar_field *out, scalar_field *in) 
{
    double *site_in, *site_out;
    _MASTER_FOR(in->type, ix) 
    {
        site_out = _FIELD_AT(out, ix);
        site_in = _FIELD_AT(in, ix);
        (*site_out) -= (*site_in);
    }
}

void sub_assign_gtransf_cpu(suNg_field* out, suNg_field* in) 
{
    suNg *site_in, *site_out;
    _MASTER_FOR(in->type, ix) 
    {
        site_out = _FIELD_AT(out, ix);
        site_in = _FIELD_AT(in, ix);
        _suNg_sub_assign((*site_out), (*site_in));
    }
}

void sub_assign_clover_term_cpu(suNfc_field* out, suNfc_field* in) 
{
    suNfc *site_out, *site_in;
    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_clover_force_cpu(suNf_field* out, suNf_field* in) 
{
    suNf *site_out, *site_in;
    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < 6; comp++) 
        {
            site_out = _6FIELD_AT(out, ix, comp);
            site_in = _6FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }
    }
}

// ** SQNORM **
double sqnorm_gfield_cpu(suNg_field *f) 
{
    suNg *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            double tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNg_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return sqnorm;
}

float sqnorm_gfield_flt_cpu(suNg_field_flt *f) 
{
    suNg_flt *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            float tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNg_sqnorm(tmp, (*site));
            sqnorm += (double)tmp;
        }
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return (float)sqnorm;
}

double sqnorm_gfield_f_cpu(suNf_field *f) 
{
    suNf *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix)
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            double tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNf_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return sqnorm;
}

float sqnorm_gfield_f_flt_cpu(suNf_field_flt *f) 
{
    suNf_flt *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            float tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNf_sqnorm(tmp, (*site));
            sqnorm += (double)tmp;
        }
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return (float)sqnorm;
}

double sqnorm_scalar_field_cpu(suNg_scalar_field *f) 
{
    suNg_vector *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        site = _FIELD_AT(f, ix);
        _vector_prod_add_assign_re_g(sqnorm, (*site), (*site));
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return sqnorm;
}

double sqnorm_avfield_cpu(suNg_av_field *f) 
{
    suNg_algebra_vector *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            double tmp;
            site = _4FIELD_AT(f, ix, comp);
            _algebra_vector_sqnorm_g(tmp, (*site));
            sqnorm += tmp;
        }
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return sqnorm;
}

double sqnorm_sfield_cpu(scalar_field *f) 
{
    double *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        site = _FIELD_AT(f, ix);
        sqnorm += (*site)*(*site);
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return sqnorm;
}

double sqnorm_gtransf_cpu(suNg_field *f) 
{
    suNg *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        double tmp;
        site = _FIELD_AT(f, ix);
        _suNg_sqnorm(tmp, (*site));
        sqnorm += tmp;
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return sqnorm;
}

double sqnorm_clover_term_cpu(suNfc_field *f) 
{
    suNfc *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            double tmp = 0.0;
            site = _4FIELD_AT(f, ix, comp);
            _suNfc_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return sqnorm;
}

double sqnorm_clover_force_cpu(suNf_field *f) 
{
    suNf *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix)
    {
        for (int comp = 0; comp < 6; comp++) 
        {
            double tmp;
            site = _6FIELD_AT(f, ix, comp);
            _suNf_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return sqnorm;
}

double sqnorm_spinor_field_f_cpu(spinor_field *f) {
    suNf_spinor *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        site = _FIELD_AT(f, ix);
        double tmp;
        _spinor_prod_re_f(tmp, (*site), (*site));
        sqnorm += tmp;
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm,1);
    #endif
    return sqnorm;
}

float sqnorm_spinor_field_f_flt_cpu(spinor_field_flt *f) {
    suNf_spinor_flt *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        site = _FIELD_AT(f, ix);
        float tmp;
        _spinor_prod_re_f(tmp, (*site), (*site));
        sqnorm += (double)tmp;
    }

    #ifdef WITH_MPI
        global_sum(&sqnorm, 1);
    #endif 
    return (float)sqnorm;
}

// Set field to zero
void zero_gfield_cpu(suNg_field *f) 
{
    int len = 4*f->type->gsize_gauge*sizeof(suNg)/sizeof(double);
    double* dbl_ptr = (double*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        dbl_ptr[i] = 0.0;
    }
}

void zero_gfield_f_cpu(suNf_field *f) 
{
    int len = 4*f->type->gsize_gauge*sizeof(suNf)/sizeof(double);
    double* dbl_ptr = (double*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        dbl_ptr[i] = 0.0;
    }
}

void zero_gfield_flt_cpu(suNg_field_flt *f) 
{
    int len = 4*f->type->gsize_gauge*sizeof(suNg_flt)/sizeof(float);
    float* flt_ptr = (float*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        flt_ptr[i] = 0.0f;
    }
}

void zero_gfield_f_flt_cpu(suNf_field_flt *f) 
{
    int len = 4*f->type->gsize_gauge*sizeof(suNf_flt)/sizeof(float);
    float* flt_ptr = (float*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        flt_ptr[i] = 0.0f;
    }
}

void zero_scalar_field_cpu(suNg_scalar_field *f) 
{
    int len = f->type->gsize_gauge*sizeof(suNg_vector)/sizeof(double);
    double* dbl_ptr = (double*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        dbl_ptr[i] = 0.0;
    }
}

void zero_avfield_cpu(suNg_av_field *f) 
{
    int len = 4*f->type->gsize_gauge*sizeof(suNg_algebra_vector)/sizeof(double);
    double* dbl_ptr = (double*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        dbl_ptr[i] = 0.0;
    }
}

void zero_gtransf_cpu(suNg_field *f) 
{
    int len = f->type->gsize_gauge*sizeof(suNg_field)/sizeof(double);
    double* dbl_ptr = (double*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        dbl_ptr[i] = 0.0;
    }
}

void zero_clover_term_cpu(suNfc_field *f) 
{
    int len = 4*f->type->gsize_gauge*sizeof(suNfc)/sizeof(double);
    double* dbl_ptr = (double*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        dbl_ptr[i] = 0.0;
    }
}

void zero_clover_force_cpu(suNf_field *f) 
{
    int len = 6*f->type->gsize_gauge*sizeof(suNf)/sizeof(double);
    double* dbl_ptr = (double*)(f->ptr);
    for (int i = 0; i < len; ++i) 
    {
        dbl_ptr[i] = 0.0;
    }
}

// ** RANDOM FIELDS FOR TESTING **
void random_spinor_field_f_cpu(spinor_field* f) 
{
    int n = f->type->gsize_spinor*sizeof(suNf_spinor)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_spinor_field_f_flt_cpu(spinor_field_flt* f) 
{
    int n = f->type->gsize_spinor*sizeof(suNf_spinor_flt)/sizeof(float);
    ranlxs((float*)(f->ptr), n);
}

void random_gfield_cpu(suNg_field* f) 
{
    int n = 4*f->type->gsize_gauge*sizeof(suNg)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_scalar_field_cpu(suNg_scalar_field *f) 
{
    int n = f->type->gsize_gauge*sizeof(suNg_vector)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_gfield_flt_cpu(suNg_field_flt *f) 
{ 
    int n = 4*f->type->gsize_gauge*sizeof(suNg_flt)/sizeof(float);
    ranlxs((float*)(f->ptr), n);
}

void random_gfield_f_cpu(suNf_field* f) 
{
    int n = 4*f->type->gsize_gauge*sizeof(suNf)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_gfield_f_flt_cpu(suNf_field_flt *f) 
{
    int n = 4*f->type->gsize_gauge*sizeof(suNf_flt)/sizeof(float);
    ranlxs((float*)(f->ptr), n);
}

void random_avfield_cpu(suNg_av_field *f) 
{
    int n = 4*f->type->gsize_gauge*sizeof(suNg_algebra_vector)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_sfield_cpu(scalar_field *f) 
{
    int n = f->type->gsize_spinor;
    ranlxd((double*)(f->ptr), n);
}

void random_gtransf_cpu(suNg_field *f) 
{
    int n = f->type->gsize_gauge*sizeof(suNg)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_clover_term_cpu(suNfc_field *f) 
{
    int n = 4*f->type->gsize_gauge*sizeof(suNfc)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_clover_force_cpu(suNf_field *f) 
{
    int n = 6*f->type->gsize_gauge*sizeof(suNf)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}




