#include "communications.h"
#include "hr_complex.h"
#include "geometry_check.h"
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "hr_complex.h"
#include <string.h>

// The following functions are primarily for testing purposes
// This is all for CPU
// TODO: This needs to work for MPI, so that tests can be run for MPI

// ** COPY **
void copy_suNg_field_cpu(suNg_field* out, suNg_field* in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg));
}

void copy_suNg_scalar_field_cpu(suNg_scalar_field *out, suNg_scalar_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg_vector));
}

void copy_suNg_field_flt_cpu(suNg_field_flt *out, suNg_field_flt *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg_flt));
}

void copy_suNf_field_f_cpu(suNf_field* out, suNf_field* in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNf));
}

void copy_suNf_field_f_flt_cpu(suNf_field_flt *out, suNf_field_flt *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNf_flt));
}

void copy_suNg_av_field_cpu(suNg_av_field *out, suNg_av_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNg_algebra_vector));
}

void copy_scalar_field_cpu(scalar_field *out, scalar_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_spinor*sizeof(double));
}

void copy_ldl_field_cpu(ldl_field *out, ldl_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(ldl_t));
}

void copy_suNfc_field_cpu(suNfc_field *out, suNf_field *in) 
{
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge*sizeof(suNfc));
}

// ** SUB ASSIGN **
void sub_assign_suNg_field_cpu(suNg_field *out, suNg_field *in) 
{
    suNg *site_out, *site_in;
    for (int ix = 0; ix < T*X*Y*Z; ix++) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNg_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_suNg_scalar_field_cpu(suNg_scalar_field *out, suNg_scalar_field *in) 
{
    suNg_vector *site_out, *site_in;
    for (int ix = 0; ix < T*X*Y*Z; ix++) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _vector_sub_assign_g((*site_out), (*site_in));
        }
    }
}

void sub_assign_suNf_field_cpu(suNf_field *out, suNf_field *in) 
{
    suNf *site_out, *site_in;
    for (int ix = 0; ix < T*X*Y*Z; ix++) 
    {
        for (int comp = 0; comp < 4; comp++) 
        {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }   
    }
}

// ** SQNORM **
double sqnorm_suNg_field_cpu(suNg_field *f) 
{
    suNg *site;
    double sqnorm = 0.0;
    for (int ix = 0; ix < T*X*Y*Z; ix++) 
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

double sqnorm_suNf_field_cpu(suNf_field *f) 
{
    suNf *site;
    double sqnorm = 0.0;
    for (int ix = 0; ix < T*X*Y*Z; ix++) 
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

