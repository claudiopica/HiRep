#include "global.h"
#include "logger.h"
#include "error.h"
#include "geometry.h"
#include "spinor_field.h"
#include "suN_types.h"
#include "utils.h"
#include <string.h>
#ifdef WITH_MPI
    #include <mpi.h>
#endif


#if defined(WITH_MPI) && defined(WITH_GPU)
#define _DECLARE_SYNC(_name, _field_type, _site_type, _size, _geom)\
    void sync_gpu_##_name(_field_type *f) \
    { \
        for (int i = 0; i < f->type->ncopies_##_geom; ++i) \
        { \
            _site_type *target = _GPU_DFIELD_BLK(f, f->type->copy_to[i], (_size));\
            _site_type *source = _GPU_DFIELD_BLK(f, f->type->copy_from[i], (_size));\
            int mem_size = (_size)*(f->type->copy_len[i])*sizeof(*(f->gpu_ptr));\
            cudaMemcpy(target, source, mem_size, cudaMemcpyDeviceToDevice); \
        } \
    }

#define _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    void start_sendrecv_##_name_gpu(_field_type *f) \
    { \
        for (int i = 0; i < f->type->nbuffers_##_geom; ++i) \
        { \
            _site_type *target = _GPU_DFIELD_BLK(f, f->type->rbuf_start[i]); \
            _site_type *source = _GPU_DFIELD_BLK(f, f->type->sbuf_start[i]); \
            int mem_size = (_size)*(f->type->sbuf_len[i])*sizeof(*(f->gpu_ptr)); \
            cudaMemcpy(target, source, mem_size, cudaMemcpyDeviceToDevice);\
        } \
    }

/* Spinor fields */
_DECLARE_SYNC(spinor_field_f, spinor_field, suNf_spinor, 1, spinor);
_DECLARE_SYNC(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, spinor);
_DECLARE_SYNC(sfield, scalar_field, double, 1, spinor);

/* Gauge fields */
_DECLARE_SYNC(gfield, suNg_field, suNg, 4, gauge);
_DECLARE_SYNC(gfield_flt, suNg_field_flt, suNg_flt, 4, gauge);
_DECLARE_SYNC(gfield_f, suNg_field, suNg, 4, gauge);
_DECLARE_SYNC(gfield_f_flt, suNg_field_flt, suNg_flt, 4, gauge);
_DECLARE_SYNC(scalar_field, suNg_scalar_field, suNg_vector, 1, gauge);
_DECLARE_SYNC(avfield, suNg_av_field, suNg_algebra_vector, 4, gauge);
_DECLARE_SYNC(gtransf, suNg_field, suNg, 1, gauge);
_DECLARE_SYNC(clover_term, suNfc_field, suNfc, 4, gauge);
_DECLARE_SYNC(clover_force, suNf_field, suNf, 6, gauge);


#undef _DECLARE_SYNC

#endif