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

#ifdef WITH_MPI

#define _GET_RECV_BUFFER(_buffer, _i, _field_dim, _type, _chars_per_site) \
    ((_buffer) + _chars_per_site * _field_dim * ((_type->rbuf_start[(_i)]) - (_type->master_shift)))

// For field_dim == 4 we need DFIELD_AT_PTR...
#define _GET_SEND_BUFFER(_buffer, _i, _field_dim, _type, _chars_per_site) \
    ((_buffer) + _chars_per_site * _field_dim * (_type->sbuf_start[(_i)]))

#define _BUFFER_FOR(_i, _nbuffers) for (int _i = 0; _i < _nbuffers; ++_i)

#define roundUp(_val, _mod) ((((_val - 1) / (_mod)) + 1) * (_mod))

#define min(a, b) (a > b) ? b : a

void hr_sendrecv_complete(int nreq, MPI_Request *field_reqs) {
    if (nreq > 0) {
        MPI_Status status[nreq];
        CHECK_MPI(MPI_Waitall(nreq, field_reqs, status));
    }
}

void hr_sendrecv(void *sendbuffer, void *recvbuffer, geometry_descriptor *type, MPI_Datatype mpi_real_type, int field_dim,
                 int size_of_real, int mpi_chunks_per_site, int nbuffers, MPI_Request *field_reqs) {
    int chars_per_site = mpi_chunks_per_site * size_of_real / sizeof(char);
    int tag_counter = 0;

    _BUFFER_FOR(i, nbuffers) {
        char *recv_buffer = (_GET_RECV_BUFFER((char *)recvbuffer, i, field_dim, type, chars_per_site));
        int recv_proc = type->rbuf_from_proc[i];
        size_t number_of_sites = roundUp(type->rbuf_len[i], THREADSIZE) / 2;
        size_t recv_size_in_dbl = field_dim * number_of_sites * mpi_chunks_per_site;

        const size_t n_messages = (recv_size_in_dbl - 1) / MAX_MSG_SIZE + 1;
        for (size_t msg_id = 0; msg_id < n_messages; msg_id++) {
            char *recv_msg = recv_buffer + msg_id * MAX_MSG_SIZE * size_of_real;
            size_t start_of_message = MAX_MSG_SIZE * msg_id;
            size_t end_of_message = min(MAX_MSG_SIZE * (msg_id + 1), recv_size_in_dbl);
            const int msg_size = (int)(end_of_message - start_of_message);
            CHECK_MPI(
                MPI_Irecv(recv_msg, msg_size, mpi_real_type, recv_proc, tag_counter, cart_comm, &(field_reqs[tag_counter])));
            tag_counter++;
        }
    }

    tag_counter = 0;
    MPI_Barrier(cart_comm);

    _BUFFER_FOR(i, nbuffers) {
        char *send_buffer = (_GET_SEND_BUFFER((char *)sendbuffer, i, field_dim, type, chars_per_site));
        int send_proc = type->sbuf_to_proc[i];
        size_t number_of_sites = roundUp(type->sbuf_len[i], THREADSIZE) / 2;
        size_t send_size_in_dbl = field_dim * number_of_sites * mpi_chunks_per_site;

        const size_t n_messages = (send_size_in_dbl - 1) / MAX_MSG_SIZE + 1;
        for (size_t msg_id = 0; msg_id < n_messages; msg_id++) {
            char *send_msg = send_buffer + msg_id * MAX_MSG_SIZE * size_of_real;
            size_t start_of_message = MAX_MSG_SIZE * msg_id;
            size_t end_of_message = min(MAX_MSG_SIZE * (msg_id + 1), send_size_in_dbl);
            const int msg_size = (int)(end_of_message - start_of_message);
            CHECK_MPI(MPI_Isend(send_msg, msg_size, mpi_real_type, send_proc, tag_counter, cart_comm,
                                &(field_reqs[2 * tag_counter])));
            tag_counter++;
        }
    }
}

#endif