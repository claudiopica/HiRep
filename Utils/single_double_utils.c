/***************************************************************************\
* Copyright (c) 2012, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File single_double_utils.c
*
* Functions for conversion from single to double precision and viceversa
*
*******************************************************************************/

#include <stdlib.h>
#include "utils.h"
#include "suN.h"
#include "error.h"
#include "global.h"
#include "spinor_field.h"


static inline void assign_f2d(double *out, float *in, unsigned int len)
{
    while(len-->0) { out[len]=(double)in[len]; }
}

static inline void assign_d2f(float *out, double *in, unsigned int len)
{
    while(len-->0) { out[len]=(float)in[len]; }
}

void assign_u2ud_cpu(void)
{
    for(int pidx=0; pidx<glattice.total_master_pieces;++pidx) {
        int ix=glattice.master_start[pidx]; //first site in the block
        double *r=(double*)(pu_gauge(ix,0));
        float *rf=(float*)(pu_gauge_flt(ix,0));
        
        int len=glattice.master_end[pidx]-glattice.master_start[pidx]+1;
        len*=4*sizeof(suNg)/sizeof(*r); //lenght of the block in real numbers
        
        assign_f2d(r, rf, len);
    }
}


void assign_ud2u_cpu(void)
{
    for(int pidx=0; pidx<glattice.total_master_pieces;++pidx) {
        int ix=glattice.master_start[pidx]; //first site in the block
        double *r=(double*)(pu_gauge(ix,0));
        float *rf=(float*)(pu_gauge_flt(ix,0));
        
        int len=glattice.master_end[pidx]-glattice.master_start[pidx]+1;
        len*=4*sizeof(suNg)/sizeof(*r); //lenght of the block in real numbers
        
        assign_d2f(rf, r, len);
    }
}

void assign_ud2u_f_cpu(void)
{
    for(int pidx=0; pidx<glattice.total_master_pieces;++pidx) {
        int ix=glattice.master_start[pidx]; //first site in the block
        double *r=(double*)(pu_gauge_f(ix,0));
        float *rf=(float*)(pu_gauge_f_flt(ix,0));
        
        int len=glattice.master_end[pidx]-glattice.master_start[pidx]+1;
        len*=4*sizeof(suNf)/sizeof(*r); //lenght of the block in real numbers
        
        assign_d2f(rf, r, len);
    }
}

void assign_u2ud_f_cpu(void)
{
    for(int pidx=0; pidx<glattice.total_master_pieces;++pidx) {
        int ix=glattice.master_start[pidx]; //first site in the block
        double *r=(double*)(pu_gauge_f(ix,0));
        float *rf=(float*)(pu_gauge_f_flt(ix,0));
        
        int len=glattice.master_end[pidx]-glattice.master_start[pidx]+1;
        len*=4*sizeof(suNf)/sizeof(*r); //lenght of the block in real numbers
        
        assign_f2d(r, rf, len);
    }
}


void assign_s2sd_cpu(spinor_field *out, spinor_field_flt *in) 
{
    _TWO_SPINORS_MATCHING(in,out);
    
    for(int pidx=0; pidx<in->type->total_master_pieces;++pidx) {
        int ix=in->type->master_start[pidx]; //first site in the block
        double *r=(double*)(_FIELD_AT(out,ix));
        float *rf=(float*)(_FIELD_AT(in,ix));
        
        int len=in->type->master_end[pidx]-in->type->master_start[pidx]+1;
        len*=sizeof(suNf_spinor)/sizeof(*r); //lenght of the block in real numbers
        
        assign_f2d(r, rf, len);
    }
}

void assign_sd2s_cpu(spinor_field_flt *out, spinor_field *in) 
{
    _TWO_SPINORS_MATCHING(in,out);
    
    for(int pidx=0; pidx<in->type->total_master_pieces;++pidx) {
        int ix=in->type->master_start[pidx]; //first site in the block
        double *r=(double*)(_FIELD_AT(in,ix));
        float *rf=(float*)(_FIELD_AT(out,ix));
        
        int len=in->type->master_end[pidx]-in->type->master_start[pidx]+1;
        len*=sizeof(suNf_spinor)/sizeof(*r); //lenght of the block in real numbers
        
        assign_d2f(rf, r, len);
    }
}

#ifndef WITH_GPU
void (*assign_s2sd)(spinor_field *out, spinor_field_flt *in)=assign_s2sd_cpu;
void (*assign_sd2s)(spinor_field_flt *out, spinor_field *in)=assign_sd2s_cpu;

void (*assign_u2ud) (void)=assign_u2ud_cpu;
void (*assign_ud2u) (void)=assign_ud2u_cpu;
void (*assign_u2ud_f)(void)=assign_u2ud_f_cpu;
void (*assign_ud2u_f) (void)=assign_ud2u_f_cpu;
#endif
