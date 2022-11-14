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


// TODO: Single precision will not work, because then we need to cast to flt
// TODO: put gpu as last suffix

#if defined(WITH_GPU) && defined(WITH_MPI)

#define _DECLARE_SYNC(_name, _field_type, _site_type, _size, _geom)\
    void sync_gpu_##_name(_field_type *f) \
    { \
        for (int i = 0; i < f->type->ncopies_##_geom; ++i) \
        { \
            _site_type *target = _GPU_DFIELD_BLK(f, f->type->copy_to[i], (_size));\
            _site_type *source = _GPU_DFIELD_BLK(f, f->type->copy_from[i], (_size));\
            int mem_size = (_size)*(f->type->copy_len[i])*sizeof(*(f->gpu_ptr));\
            CHECK_CUDA(cudaMemcpy(target, source, mem_size, cudaMemcpyDeviceToDevice)); \
        } \
    }

#define _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    void start_sendrecv_gpu_##_name(_field_type *f) \
    { \
        for (int i = 0; i < f->type->nbuffers_##_geom; ++i) \
        { \
            /* Destination Parameters */ \
            double *dest = (double*)_GPU_DFIELD_BLK(f, f->type->rbuf_start[i], (_size)); \
            int dest_proc = f->type->sbuf_to_proc[i];\
            int recv_size_in_dbl = (_size)*(f->type->rbuf_len[i])*sizeof(*(f->gpu_ptr))/sizeof(double);\
            \
            /* Origin Parameters */ \
            double *origin = (double*)_GPU_DFIELD_BLK(f, f->type->sbuf_start[i], (_size)); \
            int origin_proc = f->type->rbuf_from_proc[i];\
            int send_size_in_dbl = (_size)*(f->type->sbuf_len[i])*sizeof(*(f->gpu_ptr))/sizeof(double); \
            \
            /* Start to send */\
            CHECK_MPI(MPI_Isend(origin, send_size_in_dbl, MPI_DOUBLE, dest_proc, i, cart_comm, &(f->comm_req[2*i])));\
            \
            /* Start to receive */\
            CHECK_MPI(MPI_Irecv(dest, recv_size_in_dbl, MPI_DOUBLE, origin_proc, i, cart_comm, &(f->comm_req[2*i + 1])));\
        } \
    }

#define _DECLARE_COMPLETE_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    void complete_sendrecv_gpu_##_name(_field_type *f) \
    { \
        int nreq = 2 * f->type->nbuffers_##_geom;\
        if (nreq > 0) \
        { \
            MPI_Status status[nreq]; \
            CHECK_MPI(MPI_Waitall(nreq, f->comm_req, status)); \
        } \
    } 

#define _DECLARE_COMMS(_name, _field_type, _site_type, _size, _geom) \
    _DECLARE_SYNC(_name, _field_type, _site_type, _size, _geom) \
    _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    _DECLARE_COMPLETE_SENDRECV(_name, _field_type, _site_type, _size, _geom) 

/* Spinor fields */
_DECLARE_COMMS(spinor_field_f, spinor_field, suNf_spinor, 1, spinor);
_DECLARE_COMMS(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, spinor);
_DECLARE_COMMS(sfield, scalar_field, double, 1, spinor);

/* Gauge fields */
_DECLARE_COMMS(gfield, suNg_field, suNg, 4, gauge);
_DECLARE_COMMS(gfield_flt, suNg_field_flt, suNg_flt, 4, gauge);
_DECLARE_COMMS(gfield_f, suNf_field, suNf, 4, gauge);
_DECLARE_COMMS(gfield_f_flt, suNf_field_flt, suNf_flt, 4, gauge);
_DECLARE_COMMS(scalar_field, suNg_scalar_field, suNg_vector, 1, gauge);
_DECLARE_COMMS(avfield, suNg_av_field, suNg_algebra_vector, 4, gauge);
_DECLARE_COMMS(gtransf, suNg_field, suNg, 1, gauge);
_DECLARE_COMMS(clover_term, suNfc_field, suNfc, 4, gauge);
_DECLARE_COMMS(clover_force, suNf_field, suNf, 6, gauge);

#undef _DECLARE_COMMS
#undef _DECLARE_SYNC
#undef _DECLARE_START_SENDRECV
#undef _DECLARE_COMPLETE_SENDRECV

#endif