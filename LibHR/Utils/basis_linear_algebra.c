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


// The following functions are primarily for testing purposes
// This is all for CPU
// TODO: This needs to work for MPI, so that tests can be run for MPI
// TODO: Use randlxd and randlxs for random number generation

void rand_field_dbl(double* d, int n) {
    for (int i = 0; i < n; ++i) 
    {
        d[i] = (double)rand();
    }
}

void rand_field_flt(float* f, int n) {
    for (int i = 0; i < n; i++) 
    {
        f[i] = (double)rand();
    }
}

// ** COPY **
void copy_gfield_cpu(suNg_field* out, suNg_field* in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg));
}

void copy_scalar_field_cpu(suNg_scalar_field *out, suNg_scalar_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg_vector));
}

void copy_gfield_flt_cpu(suNg_field_flt *out, suNg_field_flt *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg_flt));
}

void copy_gfield_f_cpu(suNf_field* out, suNf_field* in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNf));
}

void copy_gfield_f_flt_cpu(suNf_field_flt *out, suNf_field_flt *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNf_flt));
}

void copy_avfield_cpu(suNg_av_field *out, suNg_av_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg_algebra_vector));
}

void copy_sfield_cpu(scalar_field *out, scalar_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_spinor*sizeof(double));
}

void copy_clover_term_cpu(suNfc_field *out, suNfc_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNfc));
}

void copy_clover_force_cpu(suNf_field *out, suNf_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNf));
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
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _vector_sub_assign_g((*site_out), (*site_in));
        }
    }
}

void sub_assign_avfield_cpu(suNg_av_field *out, suNg_av_field *in) 
{
    suNg_algebra_vector *site_in, *site_out;
    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < T*X*Y*Z; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _algebra_vector_sub_assign_g((*site_out), (*site_in));
        } 
    }
}

void sub_assign_gtransf(suNg_field* out, suNg_field* in) 
{
    sub_assign_gfield_cpu(out, in); 
}

void sub_assign_clover_term(suNfc_field* out, suNfc_field* in) 
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

void sub_assign_clover_force(suNf_field* out, suNf_field* in) 
{
    sub_assign_gfield_f_cpu(out, in);
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
    return sqnorm;
}

float sqnorm_gfield_flt_cpu(suNg_field_flt *f) 
{
    suNg_flt *site;
    float sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            float tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNg_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }
    
    return sqnorm;
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
    return sqnorm;
}

float sqnorm_gfield_f_flt_cpu(suNf_field_flt *f) 
{
    suNf_flt *site;
    float sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            float tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNf_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }
    return sqnorm;
}

double sqnorm_scalar_field_cpu(suNg_scalar_field *f) 
{
    suNg_vector *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            double tmp = 0.0;
            site = _4FIELD_AT(f, ix, comp);
            /*Sqnorm missing from suN.h macros */
            sqnorm += tmp;
        }
    }
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
    return sqnorm;
}

double sqnorm_gtransf_cpu(suNg_field *f) 
{
    return sqnorm_gfield_cpu(f);
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
}

double sqnorm_clover_force_cpu(suNf_field *f) 
{
    return sqnorm_gfield_f_cpu(f);
}

// Set field to zero
void zero_gfield_cpu(suNg_field *f) 
{
    suNg *site;
    _MASTER_FOR(f->type, ix) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site = _4FIELD_AT(f, ix, comp);
            _vector_zero_g((*site));
        }
    }
}

// ** RANDOM FIELDS FOR TESTING **
void random_gfield_cpu(suNg_field* f) 
{
    int n = f->type->gsize_gauge*sizeof(suNg)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_scalar_field_cpu(suNg_scalar_field *f) 
{
    int n = f->type->gsize_gauge*sizeof(suNg_vector)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_gfield_flt_cpu(suNg_field_flt *f) 
{
    int n = f->type->gsize_gauge*sizeof(suNg_flt)/sizeof(float);
    ranlxs((float*)(f->ptr), n);
}

void random_gfield_f_cpu(suNf_field* f) 
{
    int n = f->type->gsize_gauge*sizeof(suNf)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_gfield_f_flt_cpu(suNf_field_flt *f) 
{
    int n = f->type->gsize_gauge*sizeof(suNf_flt)/sizeof(float);
    ranlxs((float*)(f->ptr), n);
}

void random_avfield_cpu(suNg_av_field *f) 
{
    int n = f->type->gsize_gauge*sizeof(suNg_algebra_vector)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_sfield_cpu(scalar_field *f) 
{
    int n = f->type->gsize_gauge*sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_clover_term_cpu(suNfc_field *f) 
{
    int n = f->type->gsize_gauge*sizeof(suNfc)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}

void random_clover_force_cpu(suNf_field *f) 
{
    int n = f->type->gsize_gauge*sizeof(suNf)/sizeof(double);
    ranlxd((double*)(f->ptr), n);
}




