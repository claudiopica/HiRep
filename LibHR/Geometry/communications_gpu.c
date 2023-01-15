/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include "io.h"
#include "random.h"

// TODO: put gpu as last suffix
// TODO: fill buffers needs gpu suffix
// TODO: MPI error management

#ifdef WITH_GPU
    
    #ifdef WITH_MPI

        #define random_double ranlxd
        #define random_float ranlxs

        static inline void zeroes_double(double* dbl, int n) {
            for (int i = 0; i < n; ++i) { dbl[i] = 0.0; }
        }

        static inline void zeroes_float(float* flt, int n) {
            for (int i = 0; i < n; ++i) { flt[i] = 0.0f; }
        }

        #define _DECLARE_SYNC_FIELD(_name, _field_type, _site_type, _size, _geom)                                                   \
            void sync_gpu_##_name(_field_type *f)                                                                                   \
            {                                                                                                                       \
                error(geometryBoxes==NULL, 1, __func__, "geometryBoxes are not initialized.\n");                                    \
                                                                                                                                    \
                box_t *L = geometryBoxes->next;                                                                                     \
                int i = 0;                                                                                                          \
                                                                                                                                    \
                int nbuffers = f->type->nbuffers_##_geom;                                                                           \
                if (f->type==&glattice) nbuffers /= 2;                                                                              \
                while (L && i < nbuffers)                                                                                           \
                {                                                                                                                   \
                    sync_box_to_buffer_gpu_##_name(f->type, L->sendBox, (void*)f->gpu_ptr, (void*)f->sendbuf_gpu_ptr);                            \
                    L=L->next; i++;                                                                                                 \
                }                                                                                                                   \
            }

        #define _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom, _ctype, _mpictype)                            \
            void start_sendrecv_gpu_##_name(_field_type *f)                                                                         \
            {                                                                                                                       \
                sync_gpu_##_name(f);                                                                                                \
                for (int i = 0; i < f->type->nbuffers_##_geom; ++i)                                                                 \
                {                                                                                                                   \
                    int mpi_chunks_per_site = sizeof(_site_type)/sizeof(_ctype);                                                    \
                    /* Data Destination Parameters */                                                                               \
                    _ctype *recv_buffer = (_ctype*)(_BUF_GPU_DFIELD_BLK(f, i, _size));                                              \
                    int recv_proc = f->type->rbuf_from_proc[i];                                                                     \
                    int number_of_sites = f->type->rbuf_len[i];                                                                     \
                    int recv_size_in_dbl = (_size)*number_of_sites*mpi_chunks_per_site;                                             \
                                                                                                                                    \
                    /* Data Origin Parameters */                                                                                    \
                    _ctype *send_buffer = (_ctype*)(_DFIELD_AT_PTR(f->sendbuf_gpu_ptr, f->type->sbuf_start[i], 0, 0, (_size)));     \
                    int send_proc = f->type->sbuf_to_proc[i];                                                                       \
                    number_of_sites = f->type->sbuf_len[i];                                                                         \
                    int send_size_in_dbl = (_size)*number_of_sites*mpi_chunks_per_site;                                             \
                                                                                                                                    \
                    /* Communications */                                                                                            \
                    MPI_Irecv(recv_buffer, recv_size_in_dbl, _mpictype, recv_proc, i, cart_comm, &(f->comm_req[2*i+1]));            \
                    MPI_Isend(send_buffer, send_size_in_dbl, _mpictype, send_proc, i, cart_comm, &(f->comm_req[2*i]));              \
                }                                                                                                                   \
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

        #define _DECLARE_FILL_BUFFERS(_name, _field_type, _site_type, _prec_type, _size, _geom)                                     \
            void fill_buffers_##_name(_field_type *f)                                                                               \
            {                                                                                                                       \
                for (int i = 0; i < f->type->nbuffers_##_geom; ++i)                                                                 \
                {                                                                                                                   \
                    int chunks_per_site = sizeof(_site_type)/sizeof(_prec_type);                                                    \
                    int number_of_sites = f->type->rbuf_len[i];                                                                     \
                    int buffer_length = (_size)*number_of_sites*chunks_per_site;                                                    \
                                                                                                                                    \
                    _prec_type* random_array = (_prec_type*)malloc(buffer_length*sizeof(_prec_type));                               \
                    _prec_type* buffer = (_prec_type*)(_BUF_GPU_FIELD_BLK(f,i));                                                     \
                                                                                                                                    \
                    random_##_prec_type(random_array, buffer_length);                                                               \
                    cudaMemcpy(buffer, random_array, buffer_length*sizeof(_prec_type), cudaMemcpyHostToDevice);                     \
                }                                                                                                                   \
            }                                                                                                                       \
                                                                                                                                    \
            void fill_buffers_with_zeroes_##_name(_field_type *f)                                                                   \
            {                                                                                                                       \
                for (int i = 0; i < f->type->nbuffers_##_geom; ++i)                                                                 \
                {                                                                                                                   \
                    int chunks_per_site = sizeof(_site_type)/sizeof(_prec_type);                                                    \
                    int number_of_sites = f->type->rbuf_len[i];                                                                     \
                    int buffer_length = (_size)*number_of_sites*chunks_per_site;                                                    \
                                                                                                                                    \
                                                                                                                                    \
                    _prec_type* zero_array = (_prec_type*)malloc(buffer_length*sizeof(_prec_type));                                 \
                    _prec_type* buffer = (_prec_type*)(_BUF_GPU_FIELD_BLK(f,i));                                                     \
                                                                                                                                    \
                    zeroes_##_prec_type(buffer, buffer_length);                                                                     \
                    cudaMemcpy(buffer, zero_array, buffer_length*sizeof(_prec_type), cudaMemcpyHostToDevice);                       \
                }                                                                                                                   \
            } 

    #else //not WITH_MPI

        #define _DECLARE_SYNC_FIELD(_name, _field_type, _site_type, _size, _geom) \
            void sync_gpu_##_name(_field_type *f) {}

        #define _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
            void start_sendrecv_gpu_##_name(_field_type *f) {}

        #define _DECLARE_COMPLETE_SENDRECV(_name, _field_type, _site_type, _size, _geom) \
            void complete_sendrecv_gpu_##_name(_field_type *f) {}

        #define _DECLARE_FILL_BUFFERS(_name, _field_type, _site_type, _prec_type, _size, _geom) \
            void fill_buffers_##_name(_field_type *f) {}

    #endif

    #define _DECLARE_COMMS(_name, _field_type, _site_type, _size, _geom, _prec_type, _mpictype) \
            _DECLARE_SYNC_FIELD(_name, _field_type, _site_type, size, _geom) \
            _DECLARE_START_SENDRECV(_name, _field_type, _site_type, _size, _geom, _prec_type, _mpictype) \
            _DECLARE_COMPLETE_SENDRECV(_name, _field_type, _site_type, _size, _geom)  \
            _DECLARE_FILL_BUFFERS(_name, _field_type, _site_type, _prec_type, _size, _geom)

    /* Spinor fields */
    _DECLARE_COMMS(spinor_field_f, spinor_field, suNf_spinor, 1, spinor, double, MPI_DOUBLE)
    _DECLARE_COMMS(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, spinor, float, MPI_FLOAT)
    _DECLARE_COMMS(sfield, scalar_field, double, 1, spinor, double, MPI_DOUBLE)

    /* Gauge fields */
    _DECLARE_COMMS(gfield, suNg_field, suNg, 4, gauge, double, MPI_DOUBLE)
    _DECLARE_COMMS(gfield_flt, suNg_field_flt, suNg_flt, 4, gauge, float, MPI_FLOAT)
    _DECLARE_COMMS(gfield_f, suNf_field, suNf, 4, gauge, double, MPI_DOUBLE)
    _DECLARE_COMMS(gfield_f_flt, suNf_field_flt, suNf_flt, 4, gauge, float, MPI_FLOAT)
    _DECLARE_COMMS(suNg_scalar_field, suNg_scalar_field, suNg_vector, 1, gauge, double, MPI_DOUBLE)
    _DECLARE_COMMS(avfield, suNg_av_field, suNg_algebra_vector, 4, gauge, double, MPI_DOUBLE)
    _DECLARE_COMMS(gtransf, suNg_field, suNg, 1, gauge, double, MPI_DOUBLE)
    _DECLARE_COMMS(clover_ldl, ldl_field, ldl_t, 1, gauge, double, MPI_DOUBLE)
    _DECLARE_COMMS(clover_term, suNfc_field, suNfc, 4, gauge, double, MPI_DOUBLE)
    _DECLARE_COMMS(clover_force, suNf_field, suNf, 6, gauge, double, MPI_DOUBLE)

#endif