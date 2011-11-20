/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File field_pack_gpu.c
 *
 * Functions for fields allocation on GPUs
 *
 *******************************************************************************/

#include <stdlib.h>
#include "suN_types.h"
#include "error.h"
#include "memory.h"
#include "global.h"
#include "spinor_field.h"
#include "geometry.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include "gpu.h"


void spinor_field_togpuformat(spinor_field *out, spinor_field* in) {
    _DECLARE_INT_ITERATOR(ix);
    suNf_spinor *r=0;

    //check input and output type are the same
    error(out->type!=in->type,1,"spinor_field_togpuformat " __FILE__, "Spinors don't match!");
    
#ifdef UPDATE_EO
    if (in->type==&glattice) {
        // we call recursively this function twice
        // on the even and odd sublattices
        in->type=out->type=&glat_even;
        spinor_field_togpuformat(out, in);
    	in->type=out->type=&glat_odd;
        spinor_field_togpuformat(out, in);
        in->type=out->type=&glattice;
        return;
    }
#endif //UPDATE_EO
    
    _PIECE_FOR(in->type,ix) {
        const int start = in->type->master_start[_PIECE_INDEX(ix)];
        const int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1;
        complex *cout=(complex*)(_FIELD_AT(out,start));
        _SITE_FOR(in->type,ix) {
        
        	r=_FIELD_AT(in,ix);
            
            for (int j=0; j<sizeof(*r)/sizeof(complex); ++j) {
            	cout[j*N]=((complex*)(r))[j];
            }
            ++cout;
       
    	}
    }
}

void spinor_field_tocpuformat(spinor_field *out, spinor_field* in) {
    _DECLARE_INT_ITERATOR(ix);
    suNf_spinor *r=0;

    //check input and output type are the same
    error(out->type!=in->type,1,"spinor_field_togpuformat " __FILE__, "Spinors don't match!");
    
#ifdef UPDATE_EO
    if (in->type==&glattice) {
        // we call recursively this function twice
        // on the even and odd sublattices
        in->type=out->type=&glat_even;
        spinor_field_togpuformat(out, in);
    	in->type=out->type=&glat_odd;
        spinor_field_togpuformat(out, in);
        in->type=out->type=&glattice;
        return;
    }
#endif //UPDATE_EO    
    
    _PIECE_FOR(in->type,ix) {
        int start = in->type->master_start[_PIECE_INDEX(ix)];
        int N = in->type->master_end[_PIECE_INDEX(ix)]-in->type->master_start[_PIECE_INDEX(ix)]+1; 
        complex *cin=(complex*)(_FIELD_AT(in,start));
        _SITE_FOR(in->type,ix) {
            
        	r=_FIELD_AT(out,ix);
        	
            for (int j=0; j<sizeof(*r)/sizeof(complex); ++j) {
                ((complex*)(r))[j]=cin[j*N];
            }
            ++cin;
            
    	}
    }
}