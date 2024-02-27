/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file hr_sendrecv.c
 * @brief Abstracted sendrecv of buffers of given fields
 * 
*/

#include "geometry.h"
#include "memory.h"
#include "libhr_core.h"

#ifdef WITH_MPI

#define _GET_RECV_BUFFER(_buffer, _i, _field_dim, _type, _chars_per_site) \
    ((_buffer) + _chars_per_site * _field_dim * ((_type->rbuf_start[(_i)]) - (_type->master_shift)))

// For field_dim == 4 we need DFIELD_AT_PTR...
#define _GET_SEND_BUFFER(_buffer, _i, _field_dim, _type, _chars_per_site) \
    ((_buffer) + _chars_per_site * _field_dim * (_type->sbuf_start[(_i)]))

#define _BUFFER_FOR(_i, _nbuffers) for (int _i = 0; _i < _nbuffers; ++_i)

#define roundUp(_val, _mod) ((((_val - 1) / (_mod)) + 1) * (_mod))

void hr_sendrecv_complete(int nreq, MPI_Request *field_reqs) {
    if (nreq > 0) {
        MPI_Status status[nreq];
        CHECK_MPI(MPI_Waitall(nreq, field_reqs, status));
    }
}

void hr_sendrecv(void *sendbuffer, void *recvbuffer, geometry_descriptor *type, MPI_Datatype mpi_real_type, int field_dim,
                 int size_of_real, int mpi_chunks_per_site, int nbuffers, MPI_Request *field_reqs) {
    int chars_per_site = mpi_chunks_per_site * size_of_real / sizeof(char);

    _BUFFER_FOR(i, nbuffers) {
        char *recv_buffer = (_GET_RECV_BUFFER((char *)recvbuffer, i, field_dim, type, chars_per_site));
        int recv_proc = type->rbuf_from_proc[i];
        int number_of_sites = roundUp(type->rbuf_len[i], THREADSIZE) / 2;
        int recv_size_in_dbl = field_dim * number_of_sites * mpi_chunks_per_site;
        CHECK_MPI(MPI_Irecv(recv_buffer, recv_size_in_dbl, mpi_real_type, recv_proc, i, cart_comm, &(field_reqs[2 * i + 1])));
    }

    CHECK_MPI(MPI_Barrier(cart_comm));

    _BUFFER_FOR(i, nbuffers) {
        char *send_buffer = (_GET_SEND_BUFFER((char *)sendbuffer, i, field_dim, type, chars_per_site));
        int send_proc = type->sbuf_to_proc[i];
        int number_of_sites = roundUp(type->sbuf_len[i], THREADSIZE) / 2;
        int send_size_in_dbl = field_dim * number_of_sites * mpi_chunks_per_site;
        CHECK_MPI(MPI_Isend(send_buffer, send_size_in_dbl, mpi_real_type, send_proc, i, cart_comm, &(field_reqs[2 * i])));
    }
}

#endif
