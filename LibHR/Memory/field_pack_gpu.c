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

#define _DECLARE_COPY_TO_GPU_FUNCTION(_name, _field_type, _site_type, _site_vec, _dim) \
    void togpuformat_##_name(_field_type *out, _field_type *in) \
    { \
      error(out->type != in->type, 1, __FILE__, "Field geometries don't match!\n"); \
      \
      _site_type *site_in = 0; \
      \
      _PIECE_FOR(in->type, ixp)  \
      { \
        const int start = in->type->master_start[ixp]; \
        const int stride = in->type->master_end[ixp] - in->type->master_start[ixp]+1; \
        const int even_odd_offset = in->type->master_shift; \
        \
        _SITE_FOR(in->type, ixp, ix) \
        { \
          site_in = _FIELD_AT(in, ix); \
          \
          for (int comp = 0; comp < sizeof(_site_type)/sizeof(_site_vec); ++comp) \
          { \
            write_gpu_##_site_vec(even_odd_offset, (*site_in).c[comp], out, ix, comp); \
          } \
        } \
      } \
    } 


#define _DECLARE_COPY_TO_CPU_FUNCTION(_name, _field_type, _site_type, _site_vec, _dim) \
    void fromgpuformat_##_name(_field_type *out, _field_type *in) \
    { \
      error(out->type != in->type, 1, __FILE__, "Field geometries don't match!\n"); \
      \
      _site_type *site_out = 0; \
      \
      _PIECE_FOR(in->type, ixp) \
      { \
        const int start = in->type->master_start[ixp]; \
        const int stride = in->type->master_end[ixp] - in->type->master_start[ixp]+1; \
        const int even_odd_offset = in->type->master_shift; \
        _SITE_FOR(in->type, ixp, ix) \
        { \
          site_out = _FIELD_AT(in, ix); \
          \
          for (int comp = 0; comp < sizeof(_site_type)/sizeof(_site_vec); ++comp) \
          { \
            read_gpu_##_site_vec(even_odd_offset, (*site_out).c[comp], in, ix, comp); \
          } \
        } \
      } \
    } 

_DECLARE_COPY_TO_CPU_FUNCTION(spinor_field_f, spinor_field, suNf_spinor, suNf_vector, 1);
_DECLARE_COPY_TO_GPU_FUNCTION(spinor_field_f, spinor_field, suNf_spinor, suNf_vector, 1);

void spinor_field_togpuformat(spinor_field *out, spinor_field *in) {
    suNf_spinor *r=0;

    //check input and output type are the same
    error(out->type!=in->type,1,"spinor_field_togpuformat " __FILE__, "Spinors don't match!");
/*    
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
*/    
    _PIECE_FOR(in->type,ixp) {
        const int start = in->type->master_start[ixp];
        const int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1;
        hr_complex *cout=(hr_complex*)(_FIELD_AT(out,start));
        _SITE_FOR(in->type, ixp, ix) {
        
            r=_FIELD_AT(in,ix);
            
            for (int j=0; j<sizeof(*r)/sizeof(hr_complex); ++j) {
            	cout[j*N]=((hr_complex*)(r))[j];
            }
            ++cout;
       
    	}
    }
}



void spinor_field_togpuformat_flt(spinor_field_flt *out, spinor_field_flt *in) {
  suNf_spinor_flt *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"spinor_field_togpuformat_flt " __FILE__, "Spinors don't match!");

/*  
//#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    spinor_field_togpuformat_flt(out, in);
    in->type=out->type=&glat_odd;
    spinor_field_togpuformat_flt(out, in);
    in->type=out->type=&glattice;
    return;
  }
//#endif //UPDATE_EO
*/
  
  _PIECE_FOR(in->type,ixp) {
    const int start = in->type->master_start[ixp];
    const int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1;
    hr_complex_flt *cout=(hr_complex_flt*)(_FIELD_AT(out,start));
    _SITE_FOR(in->type,ixp,ix) {
      
      r=_FIELD_AT(in,ix);
      
      for (int j=0; j<sizeof(*r)/sizeof(hr_complex_flt); ++j) {
        cout[j*N]=((hr_complex_flt*)(r))[j];
      }
      ++cout;
      
    }
  }
}

void spinor_field_tocpuformat_flt(spinor_field_flt *out, spinor_field_flt *in) {
  suNf_spinor_flt *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"spinor_field_tocpuformat_flt " __FILE__, "Spinors don't match!");

/*  
//#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    spinor_field_tocpuformat_flt(out, in);
    in->type=out->type=&glat_odd;
    spinor_field_tocpuformat_flt(out, in);
    in->type=out->type=&glattice;
    return;
  }
//#endif //UPDATE_EO    
*/

  _PIECE_FOR(in->type,ixp) {
    int start = in->type->master_start[ixp];
    int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1; 
    hr_complex_flt *cin=(hr_complex_flt*)(_FIELD_AT(in,start));
    _SITE_FOR(in->type,ixp,ix) {
      
      r=_FIELD_AT(out,ix);
      
      for (int j=0; j<sizeof(*r)/sizeof(hr_complex_flt); ++j) {
        ((hr_complex_flt*)(r))[j]=cin[j*N];
      }
      ++cin;
      
    }
  }
}


void gfield_togpuformat(suNg_field *out, suNg_field *in) {
  suNg *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_togpuformat " __FILE__, "Gauge field types don't match!");
  
//#ifdef UPDATE_EO
  //if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    //in->type=out->type=&glat_even;
    //gfield_togpuformat(out, in);
    //in->type=out->type=&glat_odd;
    //gfield_togpuformat(out, in);
    //in->type=out->type=&glattice;
    //return;
  //}
//#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ixp) {
    const int start = in->type->master_start[ixp];
    const int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1;
    double *cout=(double*)(_4FIELD_AT(out,start,0));
    _SITE_FOR(in->type,ixp,ix) {
      
      r=_4FIELD_AT(in,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(double); ++j) {
        cout[j*N]=((double*)(r))[j];
      }
      ++cout;
    }
  }
}

void gfield_tocpuformat(suNg_field *out, suNg_field *in) {
  suNg *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_tocpuformat " __FILE__, "Gauge field types don't match!");
  
//#ifdef UPDATE_EO
  //if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    //in->type=out->type=&glat_even;
    //gfield_tocpuformat(out, in);
    //in->type=out->type=&glat_odd;
    //gfield_tocpuformat(out, in);
    //in->type=out->type=&glattice;
    //return;
 // }
//#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ixp) {
    const int start = in->type->master_start[ixp];
    const int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1;
    double *cin=(double*)(_4FIELD_AT(in,start,0));
    _SITE_FOR(in->type,ixp,ix) {
      
      r=_4FIELD_AT(out,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(double); ++j) {
        ((double*)(r))[j]=cin[j*N];
      }
      ++cin;
      
    }
  }
}

void gfield_togpuformat_f(suNf_field *out, suNf_field *in) {
  gfield_togpuformat((suNg_field *)out, (suNg_field *)in);
}

void gfield_tocpuformat_f(suNf_field *out, suNf_field *in) {
  gfield_tocpuformat((suNg_field *)out, (suNg_field *)in);
}

void gfield_togpuformat_flt(suNg_field_flt *out, suNg_field_flt *in) {
  suNg_flt *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_togpuformat " __FILE__, "Gauge field types don't match!");
  
//#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    gfield_togpuformat_flt(out, in);
    in->type=out->type=&glat_odd;
    gfield_togpuformat_flt(out, in);
    in->type=out->type=&glattice;
    return;
  }
//#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ixp) {
    const int start = in->type->master_start[ixp];
    const int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1;
    float *cout=(float*)(_GPU_4FIELD_AT(out,start,0));
    _SITE_FOR(in->type,ixp,ix) {
      
      r=_4FIELD_AT(in,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(float); ++j) {
        cout[j*N]=((float*)(r))[j];
      }
      ++cout;
      
    }
  }
}

void gfield_tocpuformat_flt(suNg_field_flt *out, suNg_field_flt *in) {
  suNg_flt *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_tocpuformat " __FILE__, "Gauge field types don't match!");
  
//#ifdef UPDATE_EO
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    gfield_tocpuformat_flt(out, in);
    in->type=out->type=&glat_odd;
    gfield_tocpuformat_flt(out, in);
    in->type=out->type=&glattice;
    return;
  }
//#endif //UPDATE_EO
  
  _PIECE_FOR(in->type,ixp) {
    const int start = in->type->master_start[ixp];
    const int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1;
    float *cin=(float*)(_GPU_4FIELD_AT(in,start,0));
    _SITE_FOR(in->type,ixp,ix) {
      
      r=_4FIELD_AT(out,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(double); ++j) {
        ((float*)(r))[j]=cin[j*N];
      }
      ++cin;
      
    }
  }
}

void gfield_togpuformat_f_flt(suNf_field_flt *out, suNf_field_flt *in) {
  gfield_togpuformat_flt((suNg_field_flt *)out, (suNg_field_flt *)in);
}

void gfield_tocpuformat_f_flt(suNf_field_flt *out, suNf_field_flt *in) {
  gfield_tocpuformat_flt((suNg_field_flt *)out, (suNg_field_flt *)in);
}



void avfield_togpuformat(suNg_av_field *out, suNg_av_field *in) {
  suNg_algebra_vector *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"avfield_togpuformat " __FILE__, "Algebra vector field types don't match!");
  
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    avfield_togpuformat(out, in);
    in->type=out->type=&glat_odd;
    avfield_togpuformat(out, in);
    in->type=out->type=&glattice;
    return;
  }
  
  _PIECE_FOR(in->type,ixp) {
    const int start = in->type->master_start[ixp];
    const int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1;
    double *cout=(double*)(_GPU_4FIELD_AT(out,start,0));
    _SITE_FOR(in->type,ixp,ix) {
      
      r=_4FIELD_AT(in,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(double); ++j) {
        cout[j*N]=((double*)(r))[j];
      }
      ++cout;
    }
  }
}

void avfield_tocpuformat(suNg_av_field *out, suNg_av_field *in) {
  suNg_algebra_vector *r=0;
  
  //check input and output type are the same
  error(out->type!=in->type,1,"gield_tocpuformat " __FILE__, "Gauge field types don't match!");
  
  if (in->type==&glattice) {
    // we call recursively this function twice
    // on the even and odd sublattices
    in->type=out->type=&glat_even;
    avfield_tocpuformat(out, in);
    in->type=out->type=&glat_odd;
    avfield_tocpuformat(out, in);
    in->type=out->type=&glattice;
    return;
  }
  
  _PIECE_FOR(in->type,ixp) {
    const int start = in->type->master_start[ixp];
    const int N = in->type->master_end[ixp]-in->type->master_start[ixp]+1;
    double *cin=(double*)(_GPU_4FIELD_AT(in,start,0));
    _SITE_FOR(in->type,ixp,ix) {
      
      r=_4FIELD_AT(out,ix,0);
      
      for (int j=0; j<4*sizeof(*r)/sizeof(double); ++j) {
        ((double*)(r))[j]=cin[j*N];
      }
      ++cin;
      
    }
  }
}
