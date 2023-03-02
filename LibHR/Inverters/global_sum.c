/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *
 * All rights reserved.                                                      *
 \***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include <string.h>
#include "IO/logger.h"

#ifdef MPI_TIMING
struct timeval gfstart, gfend, gfetime, sfstart, sfend, sfetime;
int gf_control = 0, sf_control = 0;
#endif

void global_sum(double *d, int n) {
#ifdef WITH_MPI
    int mpiret;
    (void)mpiret; // Remove warning of variable set but not used
    double pres[n];

#ifdef MPI_TIMING
    struct timeval start, end, etime;
    gettimeofday(&start, 0);
#endif

    mpiret = MPI_Allreduce(d, pres, n, MPI_DOUBLE, MPI_SUM, GLB_COMM);

#ifdef MPI_TIMING
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MPI TIMING", 0, "global_sum " __FILE__ " %ld sec %ld usec\n", etime.tv_sec, etime.tv_usec);
#endif

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "global_sum " __FILE__, "Cannot perform global_sum");
    }
#endif
    while (n > 0) {
        --n;
        d[n] = pres[n];
    }
#else

    /* for non mpi do nothing */
    return;
#endif
}

void global_sum_int(int *d, int n) {
#ifdef WITH_MPI
    int mpiret;
    (void)mpiret; // Remove warning of variable set but not used
    int pres[n];

#ifdef MPI_TIMING
    struct timeval start, end, etime;
    gettimeofday(&start, 0);
#endif

    mpiret = MPI_Allreduce(d, pres, n, MPI_INT, MPI_SUM, GLB_COMM);

#ifdef MPI_TIMING
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MPI TIMING", 0, "global_sum_int " __FILE__ " %ld sec %ld usec\n", etime.tv_sec, etime.tv_usec);
#endif

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "global_sum_int " __FILE__, "Cannot perform global_sum");
    }
#endif
    while (n > 0) {
        --n;
        d[n] = pres[n];
    }
#else
    /* for non mpi do nothing */
    return;
#endif
}

void global_max(double *d, int n) {
#ifdef WITH_MPI
    int mpiret;
    (void)mpiret; // Remove warning of variable set but not used
    double pres[n];

#ifdef MPI_TIMING
    struct timeval start, end, etime;
    gettimeofday(&start, 0);
#endif

    mpiret = MPI_Allreduce(d, pres, n, MPI_DOUBLE, MPI_MAX, GLB_COMM);

#ifdef MPI_TIMING
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MPI TIMING", 0, "global_max " __FILE__ " %ld sec %ld usec\n", etime.tv_sec, etime.tv_usec);
#endif

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "global_max " __FILE__, "Cannot perform global_sum");
    }
#endif
    while (n > 0) {
        --n;
        d[n] = pres[n];
    }
#else
    /* for non mpi do nothing */
    return;
#endif
}

void global_min(double *d, int n) {
#ifdef WITH_MPI
    int mpiret;
    (void)mpiret; // Remove warning of variable set but not used
    double pres[n];

#ifdef MPI_TIMING
    struct timeval start, end, etime;
    gettimeofday(&start, 0);
#endif

    mpiret = MPI_Allreduce(d, pres, n, MPI_DOUBLE, MPI_MIN, GLB_COMM);

#ifdef MPI_TIMING
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MPI TIMING", 0, "global_max " __FILE__ " %ld sec %ld usec\n", etime.tv_sec, etime.tv_usec);
#endif

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "global_min " __FILE__, "Cannot perform global_sum");
    }
#endif
    while (n > 0) {
        --n;
        d[n] = pres[n];
    }
#else
    /* for non mpi do nothing */
    return;
#endif
}

void bcast(double *d, int n) {
#ifdef WITH_MPI
    int mpiret;
    (void)mpiret; // Remove warning of variable set but not used

#ifdef MPI_TIMING
    struct timeval start, end, etime;
    gettimeofday(&start, 0);
#endif

    mpiret = MPI_Bcast(d, n, MPI_DOUBLE, 0, GLB_COMM);

#ifdef MPI_TIMING
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MPI TIMING", 0, "bcast " __FILE__ " %ld sec %ld usec\n", etime.tv_sec, etime.tv_usec);
#endif

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "bcast " __FILE__, "Cannot perform global_sum");
    }
#endif

#else
    /* for non mpi do nothing */
    return;
#endif
}

void bcast_int(int *i, int n) {
#ifdef WITH_MPI
    int mpiret;
    (void)mpiret; // Remove warning of variable set but not used

#ifdef MPI_TIMING
    struct timeval start, end, etime;
    gettimeofday(&start, 0);
#endif

    mpiret = MPI_Bcast(i, n, MPI_INT, 0, GLB_COMM);

#ifdef MPI_TIMING
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MPI TIMING", 0, "bcast " __FILE__ " %ld sec %ld usec\n", etime.tv_sec, etime.tv_usec);
#endif

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "bcast " __FILE__, "Cannot perform global_sum");
    }
#endif

#else
    /* for non mpi do nothing */
    return;
#endif
}

void global_max_flt(float *d, int n) {
#ifdef WITH_MPI
    int mpiret;
    (void)mpiret; // Remove warning of variable set but not used
    float pres[n];

#ifdef MPI_TIMING
    struct timeval start, end, etime;
    gettimeofday(&start, 0);
#endif

    mpiret = MPI_Allreduce(d, pres, n, MPI_FLOAT, MPI_MAX, GLB_COMM);

#ifdef MPI_TIMING
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MPI TIMING", 0, "global_max_flt " __FILE__ " %ld sec %ld usec\n", etime.tv_sec, etime.tv_usec);
#endif

#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "global_max " __FILE__, "Cannot perform global_sum");
    }
#endif
    while (n > 0) {
        --n;
        d[n] = pres[n];
    }
#else
    /* for non mpi do nothing */
    return;
#endif
}
