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

    
    _PIECE_FOR(in->type,ix) {
        int start = in->type->master_start[_PIECE_INDEX(ix)];
        int N = in->type->master_end[_PIECE_INDEX(ix)] -  in->type->master_start[_PIECE_INDEX(ix)] + 1; 
        _SITE_FOR(in->type,ix) {
        
        	r=_FIELD_AT(in,ix);
        	
            for (int j=0; j<sizeof(*r)/sizeof(complex); ++j) {
            	((complex*)(out->ptr))[ix+j*N]=((complex*)(r))[j];
            }
       
    	}
    }
}

void spinor_field_tocpuformat(spinor_field *out, spinor_field* in) {
    _DECLARE_INT_ITERATOR(ix);
    suNf_spinor *r=0;
    
    
    _PIECE_FOR(in->type,ix) {
        int start = in->type->master_start[_PIECE_INDEX(ix)];
        int N = in->type->master_end[_PIECE_INDEX(ix)] -  in->type->master_start[_PIECE_INDEX(ix)] + 1; 
        _SITE_FOR(in->type,ix) {
            
        	r=_FIELD_AT(out,ix);
        	
            for (int j=0; j<sizeof(*r)/sizeof(complex); ++j) {
            	((complex*)(r))[j]=((complex*)(in->ptr))[ix+j*N];
            }
            
    	}
    }
}