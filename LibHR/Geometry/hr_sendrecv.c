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
#include "io.h"

// For reference https://github.com/NikashKumar/Multi-Threaded-Producer-Consumer-Work-Queue-Handling/blob/master/Consumer_Producer.c
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#define NO_OF_COMMS_THREADS 1
pthread_t thread[NO_OF_COMMS_THREADS];
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
pthread_barrier_t barrier;

#ifdef WITH_MPI

#define _GET_RECV_BUFFER(_buffer, _i, _field_dim, _type, _chars_per_site) \
    ((_buffer) + _chars_per_site * _field_dim * ((_type->rbuf_start[(_i)]) - (_type->master_shift)))

// For field_dim == 4 we need DFIELD_AT_PTR...
#define _GET_SEND_BUFFER(_buffer, _i, _field_dim, _type, _chars_per_site) \
    ((_buffer) + _chars_per_site * _field_dim * (_type->sbuf_start[(_i)]))

#define _BUFFER_FOR(_i, _nbuffers) for (int _i = 0; _i < _nbuffers; ++_i)

#define roundUp(_val, _mod) ((((_val - 1) / (_mod)) + 1) * (_mod))

static int finalize = 0;
static int ready = 0;
static int init = 0;

// small internal locking mechanism
static int lock = 0;

static void lock_comms() {
    lock = 1;
}

static void release_comms() {
    lock = 0;
}

static int is_locked() {
    return lock;
}

// Shared struct & static instance
typedef struct {
    char *send_buffer;
    char *recv_buffer;
    geometry_descriptor *type;
    MPI_Datatype mpi_real_type;
    int field_dim;
    int size_of_real;
    int mpi_chunks_per_site;
    int nbuffers;
    MPI_Request *field_reqs;
} comms_args;

static comms_args *c_args = NULL;

// Communication routine using shared struct
void communicate(comms_args *l) {
    // Initialize correct CUDA device context on consumer thread
    if (!init) {
        CHECK_CUDA(cudaSetDevice(LID));
        init = 1;
    }

    // Standard blocking communications
    int chars_per_site = l->mpi_chunks_per_site * l->size_of_real / sizeof(char);
    _BUFFER_FOR(i, l->nbuffers) {
        char *recv_buffer = (_GET_RECV_BUFFER((char *)l->recv_buffer, i, l->field_dim, l->type, chars_per_site));
        int recv_proc = l->type->rbuf_from_proc[i];
        int number_of_sites = roundUp(l->type->rbuf_len[i], THREADSIZE) / 2;
        int recv_size_in_dbl = l->field_dim * number_of_sites * l->mpi_chunks_per_site;

        char *send_buffer = (_GET_SEND_BUFFER((char *)l->send_buffer, i, l->field_dim, l->type, chars_per_site));
        int send_proc = l->type->sbuf_to_proc[i];
        number_of_sites = roundUp(l->type->sbuf_len[i], THREADSIZE) / 2;
        int send_size_in_dbl = l->field_dim * number_of_sites * l->mpi_chunks_per_site;

        CHECK_MPI(MPI_Sendrecv(send_buffer, send_size_in_dbl, MPI_DOUBLE, send_proc, i, recv_buffer, recv_size_in_dbl,
                               MPI_DOUBLE, recv_proc, i, cart_comm, MPI_STATUS_IGNORE));
    }
}

void *wait_for_signal(void *argv) {
    comms_args *l = (comms_args *)argv;
    while (1) {
        pthread_mutex_lock(&mutex);
        while (!ready) {
            pthread_cond_wait(&cond, &mutex);
        }
        if (!finalize) {
            // If the program isn't finalizing, the signal indicates communications
            communicate(l);
            ready = 0;
            pthread_mutex_unlock(&mutex);
            pthread_barrier_wait(&barrier);
        } else {
            // If the program is finalizing then break the while loop
            pthread_mutex_unlock(&mutex);
            break;
        }
    }
    return 0;
}

void signal() {
    pthread_mutex_lock(&mutex);
    ready = 1;
    pthread_cond_broadcast(&cond);
    pthread_mutex_unlock(&mutex);
}

int spawn_threads() {
    pthread_attr_t attributes;
    pthread_attr_init(&attributes);
    pthread_barrier_init(&barrier, NULL, 2);
    if (pthread_create(&thread[0], &attributes, wait_for_signal, c_args)) {
        error(1, 1, __func__, "Unable to spawn communication thread(s)\n");
        printf("Error creating pthreads\n");
    }
    return 0;
}

void hr_sendrecv_complete(int nreq, MPI_Request *field_reqs) {
    pthread_barrier_wait(&barrier);
    release_comms();
}

void init_hr_comms() {
    if (c_args == NULL) { c_args = (comms_args *)malloc(sizeof(comms_args)); }
    spawn_threads();
}

void finalize_hr_comms() {
    pthread_mutex_lock(&mutex);
    finalize = 1;
    pthread_mutex_unlock(&mutex);
    signal();
    pthread_join(thread[0], NULL);
}

void hr_sendrecv(void *sendbuffer, void *recvbuffer, geometry_descriptor *type, MPI_Datatype mpi_real_type, int field_dim,
                 int size_of_real, int mpi_chunks_per_site, int nbuffers, MPI_Request *field_reqs) {
    pthread_mutex_lock(&mutex);
    c_args->send_buffer = sendbuffer;
    c_args->recv_buffer = recvbuffer;
    c_args->type = type;
    c_args->mpi_real_type = mpi_real_type;
    c_args->field_dim = field_dim;
    c_args->size_of_real = size_of_real;
    c_args->mpi_chunks_per_site = mpi_chunks_per_site;
    c_args->nbuffers = nbuffers;
    c_args->field_reqs = field_reqs;

    if (!is_locked()) {
        lock_comms();
        pthread_mutex_unlock(&mutex);
        signal();
    } else {
        error(1, 1, __func__, "Simultaneous reduced communications attempted\n\n");
    }
}

#endif
