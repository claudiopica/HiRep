#include "global.h"
#include "logger.h"
#include "error.h"
#include "geometry.h"
#include "random.h"
#include "spinor_field.h"
#include "suN_types.h"
#include "utils.h"
#include <string.h>
#include "new_geometry.h"

// TODO: Single precision will not work, because then we need to cast to flt
// TODO: put gpu as last suffix
// TODO: fill buffers needs gpu suffix

#ifdef WITH_GPU
#ifdef WITH_MPI

#define random_double ranlxd
#define random_float ranlxs

void zeroes_double(double* dbl, int n) 
{
    for (int i = 0; i < n; ++i) 
    {
        dbl[i] = 0.0;
    }
}

void zeroes_float(float* flt, int n) 
{
    for (int i = 0; i < n; ++i) 
    {
        flt[i] = 0.0f;
    }
}

#define _DECLARE_SYNC_FIELD(_name, _type, _geom) \
    static void sync_field_to_buffer_gpu_##_name(geometry_descriptor *gd, \
                        _type *lattice,  \
                        void *sendbuf) \
    { \
        error(geometryBoxes==NULL, 1, __func__, "geometryBoxes are not initialized.\n"); \
        \
        /* Query first buffer box in the list */ \
        box_t *L = geometryBoxes->next; \
        \
        /* Iterate over all boxes*/ \
        /* The i-counter is necessary, because spinor-like and gauge-like fields have */ \
        /* different numbers of buffers */ \
        int i = 0; \
        while (L && i < gd->nbuffers_##_geom) \
        { \
            sync_box_to_buffer_gpu_##_name(gd, L->sendBox, lattice, sendbuf); \
            L=L->next; i++; \
        } \
    } 

#define _DECLARE_SYNC_FUNCTIONS(_name, _type, _size, _geom) \
    _DECLARE_SYNC_FIELD(_name, _type, _geom)

_DECLARE_SYNC_FUNCTIONS(spinor_field_f, suNf_spinor, 1, spinor);
_DECLARE_SYNC_FUNCTIONS(spinor_field_f_flt, suNf_spinor_flt, 1, spinor);
_DECLARE_SYNC_FUNCTIONS(sfield, double, 1, spinor);

_DECLARE_SYNC_FUNCTIONS(gfield, suNg, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_flt, suNg_flt, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_f, suNf, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_f_flt, suNf_flt, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(suNg_scalar_field, suNg_vector, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(avfield, suNg_algebra_vector, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gtransf, suNg, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_ldl, ldl_t, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_term, suNfc, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_force, suNf, 6, gauge);

#define _DECLARE_SYNC(_name, _field_type, _site_type, _size, _geom) \
    void sync_gpu_##_name(_field_type *f) \
    { \
        sync_field_to_buffer_gpu_##_name(f->type, f->gpu_ptr, f->sendbuf_gpu_ptr); \
    }

#define _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    void start_sendrecv_gpu_##_name(_field_type *f) \
    { \
        sync_gpu_##_name(f); \
        MPI_Status status[f->type->nbuffers_##_geom];\
        for (int i = 0; i < f->type->nbuffers_##_geom; ++i) \
        { \
            /* Destination Parameters */ \
            double *recv_buffer = (double*)(f->gpu_ptr + (_size)*f->type->rbuf_start[i]);\
            int recv_proc = f->type->rbuf_from_proc[i];\
            int recv_size_in_dbl = (_size)*(f->type->rbuf_len[i])*sizeof(*(f->gpu_ptr))/sizeof(double);\
            \
            /* Origin Parameters */ \
            double *send_buffer = (double*)(f->sendbuf_gpu_ptr + (_size)*f->type->sbuf_start[i]);\
            int send_proc = f->type->sbuf_to_proc[i];\
            int send_size_in_dbl = (_size)*(f->type->sbuf_len[i])*sizeof(*(f->gpu_ptr))/sizeof(double); \
            \
            /* Start to receive */\
            int mpiret; (void)mpiret; \
            mpiret = MPI_Irecv(recv_buffer, recv_size_in_dbl, MPI_DOUBLE, recv_proc, i, cart_comm, &(f->comm_req[2*i+1]));\
            if (mpiret != MPI_SUCCESS) { \
                char mesg[MPI_MAX_ERROR_STRING]; \
                int mesglen; \
                MPI_Error_string(mpiret, mesg, &mesglen); \
                lprintf("MPI", 0, "ERROR: %s\n", mesg); \
                error(1, 1, "global_sum_gpu " __FILE__, ": Cannot start recv comms.\n"); \
            } \
            \
            /* Start to send */\
            mpiret = MPI_Isend(send_buffer, send_size_in_dbl, MPI_DOUBLE, send_proc, i, cart_comm, &(f->comm_req[2*i]));\
            if (mpiret != MPI_SUCCESS) { \
                char mesg[MPI_MAX_ERROR_STRING]; \
                int mesglen; \
                MPI_Error_string(mpiret, mesg, &mesglen); \
                lprintf("MPI", 0, "ERROR: %s\n", mesg); \
                error(1, 1, "global_sum_gpu " __FILE__, ": Cannot start send comms.\n");\
            } \
        } \
    }

#define _DECLARE_FILL_BUFFERS(_name, _field_type, _prec_type, _size, _geom) \
    void fill_buffers_##_name(_field_type *f) \
    { \
        for (int i = 0; i < f->type->nbuffers_##_geom; ++i) \
        { \
            int rlen = (_size)*(f->type->rbuf_len[i])*sizeof(*(f->gpu_ptr))/sizeof(_prec_type); \
            _prec_type* buf = (_prec_type*)malloc(rlen*sizeof(*(f->gpu_ptr)));\
            random_##_prec_type(buf, rlen); \
            cudaMemcpy((_prec_type*)(f->gpu_ptr + (_size)*f->type->rbuf_start[i]), \
                                buf, rlen*sizeof(_prec_type), cudaMemcpyHostToDevice);\
        }\
    }\
    \
    void fill_buffers_with_zeroes_##_name(_field_type *f) \
    { \
        for (int i = 0; i < f->type->nbuffers_##_geom; ++i) \
        { \
            int rlen = (_size)*(f->type->rbuf_len[i])*sizeof(*(f->gpu_ptr))/sizeof(_prec_type); \
            _prec_type* buf = (_prec_type*)malloc(rlen*sizeof(*(f->gpu_ptr)));\
            zeroes_##_prec_type(buf, rlen); \
            cudaMemcpy((_prec_type*)(f->gpu_ptr + (_size)*f->type->rbuf_start[i]),  \
                                buf, rlen*sizeof(_prec_type), cudaMemcpyHostToDevice);\
        }\
    } 

#define _DECLARE_COMPLETE_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    void complete_sendrecv_gpu_##_name(_field_type *f) \
    { \
        int nreq = 2 * f->type->nbuffers_##_geom;\
        if (nreq > 0) \
        { \
            MPI_Status status[nreq]; \
            MPI_Waitall(nreq, f->comm_req, status); \
        } \
    } 

#define _DECLARE_COMMS(_name, _field_type, _site_type, _size, _geom, _prec_type) \
    _DECLARE_SYNC(_name, _field_type, _site_type, _size, _geom) \
    _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    _DECLARE_COMPLETE_SENDRECV(_name, _field_type, _site_type, _size, _geom)  \
    _DECLARE_FILL_BUFFERS(_name, _field_type, _prec_type, _size, _geom)

/* Spinor fields */
_DECLARE_COMMS(spinor_field_f, spinor_field, suNf_spinor, 1, spinor, double);
_DECLARE_COMMS(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, spinor, float);
_DECLARE_COMMS(sfield, scalar_field, double, 1, spinor, double);

/* Gauge fields */
_DECLARE_COMMS(gfield, suNg_field, suNg, 4, gauge, double);
_DECLARE_COMMS(gfield_flt, suNg_field_flt, suNg_flt, 4, gauge, float);
_DECLARE_COMMS(gfield_f, suNf_field, suNf, 4, gauge, double);
_DECLARE_COMMS(gfield_f_flt, suNf_field_flt, suNf_flt, 4, gauge, float);
_DECLARE_COMMS(suNg_scalar_field, suNg_scalar_field, suNg_vector, 1, gauge, double);
_DECLARE_COMMS(avfield, suNg_av_field, suNg_algebra_vector, 4, gauge, double);
_DECLARE_COMMS(gtransf, suNg_field, suNg, 1, gauge, double);
_DECLARE_COMMS(clover_ldl, ldl_field, ldl_t, 1, gauge, double);
_DECLARE_COMMS(clover_term, suNfc_field, suNfc, 4, gauge, double);
_DECLARE_COMMS(clover_force, suNf_field, suNf, 6, gauge, double);

#undef _DECLARE_COMMS
#undef _DECLARE_SYNC
#undef _DECLARE_START_SENDRECV
#undef _DECLARE_COMPLETE_SENDRECV
#undef random_double
#undef random_float
<<<<<<< HEAD
=======

#else //not WITH_MPI

#define _DECLARE_SYNC(_name, _field_type, _site_type, _size, _geom) \
    void sync_gpu_##_name(_field_type *f) {}

#define _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    void start_sendrecv_gpu_##_name(_field_type *f) {}

#define _DECLARE_COMPLETE_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    void complete_sendrecv_gpu_##_name(_field_type *f) {}

#define _DECLARE_FILL_BUFFERS(_name, _field_type, _prec_type, _size, _geom) \
    void fill_buffers_##_name(_field_type *f) {}

#define _DECLARE_COMMS(_name, _field_type, _site_type, _size, _geom, _prec_type) \
    _DECLARE_SYNC(_name, _field_type, _site_type, _size, _geom) \
    _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
    _DECLARE_COMPLETE_SENDRECV(_name, _field_type, _site_type, _size, _geom)  \
    _DECLARE_FILL_BUFFERS(_name, _field_type, _prec_type, _size, _geom)


/* Spinor fields */
_DECLARE_COMMS(spinor_field_f, spinor_field, suNf_spinor, 1, spinor, double);
_DECLARE_COMMS(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, spinor, float);
_DECLARE_COMMS(sfield, scalar_field, double, 1, spinor, double);

/* Gauge fields */
_DECLARE_COMMS(gfield, suNg_field, suNg, 4, gauge, double);
_DECLARE_COMMS(gfield_flt, suNg_field_flt, suNg_flt, 4, gauge, float);
_DECLARE_COMMS(gfield_f, suNf_field, suNf, 4, gauge, double);
_DECLARE_COMMS(gfield_f_flt, suNf_field_flt, suNf_flt, 4, gauge, float);
_DECLARE_COMMS(suNg_scalar_field, suNg_scalar_field, suNg_vector, 1, gauge, double);
_DECLARE_COMMS(avfield, suNg_av_field, suNg_algebra_vector, 4, gauge, double);
_DECLARE_COMMS(gtransf, suNg_field, suNg, 1, gauge, double);
_DECLARE_COMMS(clover_ldl, ldl_field, ldl_t, 1, gauge, double);
_DECLARE_COMMS(clover_term, suNfc_field, suNfc, 4, gauge, double);
_DECLARE_COMMS(clover_force, suNf_field, suNf, 6, gauge, double);

#endif

>>>>>>> remotes/upstream/HiRep-CUDA
#endif