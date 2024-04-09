/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File archive.c
 *
 * Write and read routines for archiving configurations
 *
 *******************************************************************************/
#include "libhr.h"
/*#include "Core/global.h"
#include "Geometry/communications.h"
#include "Utils/boundary_conditions.h"
#include "Utils/timing.h"
#include "Update/avr_plaquette.h"
#include "Random/ranlux.h"
#include "io.h"
#include "geometry.h"
//#include "inverters.h"
#include "Inverters/global_sum.h"*/

void write_gauge_field(char filename[]) {
#if NG == 2 && !defined(WITH_QUATERNIONS)
    write_gauge_field_su2q(filename);
#else
    write_gauge_field_matrix(filename);
#endif
}

void write_gauge_field_matrix(char filename[]) {
    FILE *fp = NULL;
    int g[4], p[4];
    double *buff = NULL;
    int pid = 0;
    int zsize, rz;
    double plaq;

#ifdef WITH_MPI
    /* MPI variables */
    MPI_Group wg, cg;
    MPI_Status st;
    int cid;
    int mpiret;
    (void)mpiret;
#endif

#ifndef ALLOCATE_REPR_GAUGE_FIELD
    complete_sendrecv_suNg_field(u_gauge);
    apply_BCs_on_represented_gauge_field(); //Save the link variables with periodic boundary conditions
#endif

    plaq = avr_plaquette(); /* to use as a checksum in the header */

#ifdef WITH_GPU
    copy_from_gpu_suNg_field(u_gauge);
#endif
    if (PID == 0) {
        int d[5] = { NG, GLB_T, GLB_X, GLB_Y, GLB_Z };
        error((fp = fopen(filename, "wb")) == NULL, 1, "write_gauge_field", "Failed to open file for writing");
        /* write NG and global size */
        error(fwrite_BE_int(d, (size_t)(5), fp) != (5), 1, "write_gauge_field", "Failed to write gauge field geometry");
        /* write average plaquette */
        error(fwrite_BE_double(&plaq, (size_t)(1), fp) != (1), 1, "write_gauge_field", "Failed to write gauge field plaquette");
    }

#ifdef WITH_MPI
    MPI_Comm_group(GLB_COMM, &wg);
    MPI_Comm_group(cart_comm, &cg);
#endif

    Timer clock;
    timer_set(&clock);

    zsize = GLB_Z / NP_Z;
    rz = GLB_Z - zsize * NP_Z;
    if (!four_fermion_active) {
        buff = malloc(sizeof(suNg) * 4 * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));
    } else {
        buff = malloc((sizeof(suNg) * 4 + 2 * sizeof(double)) * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));
    }

    g[3] = 0;
    for (g[0] = 0; g[0] < GLB_T; ++g[0]) { /* loop over T, X and Y direction */
        for (g[1] = 0; g[1] < GLB_X; ++g[1]) {
            for (g[2] = 0; g[2] < GLB_Y; ++g[2]) {
#ifdef WITH_MPI
                glb_to_proc(g, p); /* get the processor coordinate */
#endif
                for (p[3] = 0; p[3] < NP_Z; ++p[3]) { /* loop over processors in Z direction */
                    int bsize;
                    if (!four_fermion_active) {
                        bsize = sizeof(suNg) / sizeof(double) * 4 *
                                (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */
                    } else {
                        bsize = (sizeof(suNg) / sizeof(double) * 4 + 2) *
                                (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */
                    }

#ifdef WITH_MPI
                    MPI_Cart_rank(cart_comm, p, &cid);
                    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
                    if (pid == PID) { /* fill link buffer */
                        int lsite[4];
                        suNg *cm;

                        /* convert global to local coordinate */
                        origin_coord(lsite);
                        lsite[0] = g[0] - lsite[0];
                        lsite[1] = g[1] - lsite[1];
                        lsite[2] = g[2] - lsite[2];

                        /* fill buffer */
                        cm = (suNg *)buff;
                        for (lsite[3] = 0; lsite[3] < Z; ++lsite[3]) { /* loop on local Z */
                            int ix = ipt(lsite[0], lsite[1], lsite[2], lsite[3]);
                            suNg *pm = pu_gauge(ix, 0);
                            *(cm++) = *(pm++); /* copy 4 directions */
                            *(cm++) = *(pm++);
                            *(cm++) = *(pm++);
                            *(cm++) = *(pm);

                            if (four_fermion_active) {
                                double *s_buff = (double *)cm;
                                double *s = _FIELD_AT(ff_sigma, ix);
                                double *pl = _FIELD_AT(ff_pi, ix);
                                *(s_buff++) = *s;
                                *(s_buff++) = *pl;
                                cm = (suNg *)s_buff;
                            }
                        }
                    }
#ifdef WITH_MPI
                    MPI_Barrier(GLB_COMM);
                    if (pid != 0) { /* do send/receive only if the data is not on PID 0 */
                        /* send buffer */
                        if (pid == PID) {
#ifndef NDEBUG
                            error(Z != (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)), 1, "write_gauge_field",
                                  "Local lattice size mismatch!");
#endif
                            mpiret = MPI_Send(buff, bsize, MPI_DOUBLE, 0, 999, GLB_COMM);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                error(1, 1, "write_gauge_field " __FILE__, "Cannot send buffer");
                            }
#endif
                        }
                        /* receive buffer */
                        if (PID == 0) {
                            mpiret = MPI_Recv(buff, bsize, MPI_DOUBLE, pid, 999, GLB_COMM, &st);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                if (st.MPI_ERROR != MPI_SUCCESS) {
                                    MPI_Error_string(st.MPI_ERROR, mesg, &mesglen);
                                    lprintf("MPI", 0, "Source [%d] Tag [%] ERROR: %s\n", st.MPI_SOURCE, st.MPI_TAG, mesg);
                                }
                                error(1, 1, "write_gauge_field " __FILE__, "Cannot receive buffer");
                            }
#endif
                        }
                    }
#endif

                    /* write buffer to file */
                    if (PID == 0) {
                        error(fwrite_BE_double(buff, (size_t)(bsize), fp) != (bsize), 1, "write_gauge_field",
                              "Failed to write gauge field to file");
                    }

                } /* end loop over processors in Z direction */
            }
        }
    } /* end loop over T, X and Y direction */

    if (PID == 0) { fclose(fp); }
    free(buff);

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("IO", 0, "Configuration [%s] saved [%lf sec]\n", filename, elapsed_sec);

#ifdef WITH_MPI
    MPI_Barrier(GLB_COMM);
#endif

#ifndef ALLOCATE_REPR_GAUGE_FIELD
    complete_sendrecv_suNg_field(u_gauge);
    apply_BCs_on_represented_gauge_field(); //Restore the right boundary conditions
#endif
}

void read_gauge_field(char filename[]) {
#if (NG == 2) && !defined(WITH_QUATERNIONS)
    read_gauge_field_su2(filename);
#else
    read_gauge_field_matrix(filename);
#endif
    apply_BCs_on_fundamental_gauge_field();

#ifndef ALLOCATE_REPR_GAUGE_FIELD
    u_gauge_f = (suNf_field *)((void *)u_gauge);
    apply_BCs_on_represented_gauge_field();
#endif
}

void read_gauge_field_matrix(char filename[]) {
    FILE *fp = NULL;
    int g[4], p[4];
    double *buff = NULL;
    int pid = 0;
    int zsize, rz;
    double plaq, testplaq;

#ifdef WITH_MPI
    /* MPI variables */
    MPI_Group wg, cg;
    MPI_Status st;
    int cid;
    int mpiret;
    (void)mpiret;
#endif

    Timer clock;
    timer_set(&clock);

    if (PID == 0) {
        int d[5] = { 0 }; /* contains NG,GLB_T,GLB_X,GLB_Y,GLB_Z */
        error((fp = fopen(filename, "rb")) == NULL, 1, "read_gauge_field", "Failed to open file for reading");
        /* read NG and global size */
        error(fread_BE_int(d, (size_t)(5), fp) != (5), 1, "read_gauge_field", "Failed to read gauge field geometry");
        /* Check Gauge group and Lattice dimesions */
        if (NG != d[0]) {
            lprintf("ERROR", 0, "Read value of NG [%d] do not match this code [NG=%d].\nPlease recompile.\n", d[0], NG);
            error(1, 1, "read_gauge_field " __FILE__, "Gauge group mismatch");
        }
        if (GLB_T != d[1] || GLB_X != d[2] || GLB_Y != d[3] || GLB_Z != d[4]) {
            lprintf("ERROR", 0, "Read value of global lattice size (%d,%d,%d,%d) do not match input file (%d,%d,%d,%d).\n",
                    d[1], d[2], d[3], d[4], GLB_T, GLB_X, GLB_Y, GLB_Z);
            error(1, 1, "read_gauge_field " __FILE__, "Lattice dimension mismatch");
        }
        error(fread_BE_double(&plaq, (size_t)(1), fp) != (1), 1, "read_gauge_field", "Failed to read gauge field plaquette");
    }

#ifdef WITH_MPI
    MPI_Comm_group(GLB_COMM, &wg);
    MPI_Comm_group(cart_comm, &cg);
#endif

    zsize = GLB_Z / NP_Z;
    rz = GLB_Z - zsize * NP_Z;
    if (!four_fermion_active) {
        buff = malloc(sizeof(suNg) * 4 * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));
    } else {
        buff = malloc((sizeof(suNg) * 4 + 2 * sizeof(double)) * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));
    }

    g[3] = 0;
    for (g[0] = 0; g[0] < GLB_T; ++g[0]) { /* loop over T, X and Y direction */
        for (g[1] = 0; g[1] < GLB_X; ++g[1]) {
            for (g[2] = 0; g[2] < GLB_Y; ++g[2]) {
#ifdef WITH_MPI
                glb_to_proc(g, p); /* get the processor coordinate */
#endif
                for (p[3] = 0; p[3] < NP_Z; ++p[3]) { /* loop over processors in Z direction */
                    int bsize;
                    if (!four_fermion_active) {
                        bsize = sizeof(suNg) / sizeof(double) * 4 *
                                (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */
                    } else {
                        bsize = (sizeof(suNg) / sizeof(double) * 4 + 2) *
                                (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */
                    }

#ifdef WITH_MPI
                    MPI_Cart_rank(cart_comm, p, &cid);
                    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
                    /* read buffer from file */
                    if (PID == 0) {
                        error(fread_BE_double(buff, (size_t)(bsize), fp) != (bsize), 1, "read_gauge_field",
                              "Failed to read gauge field from file");
                    }

#ifdef WITH_MPI
                    /* MPI_Barrier(GLB_COMM); */
                    if (pid != 0) { /* do send/receive only if the data is not on PID 0 */
                        /* send buffer */
                        if (PID == 0) {
                            mpiret = MPI_Send(buff, bsize, MPI_DOUBLE, pid, 998, GLB_COMM);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                error(1, 1, "read_gauge_field " __FILE__, "Cannot send buffer");
                            }
#endif
                        }
                        /* receive buffer */
                        if (PID == pid) {
#ifndef NDEBUG
                            error(Z != (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)), 1, "read_gauge_field",
                                  "Local lattice size mismatch!");
#endif
                            mpiret = MPI_Recv(buff, bsize, MPI_DOUBLE, 0, 998, GLB_COMM, &st);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                if (st.MPI_ERROR != MPI_SUCCESS) {
                                    MPI_Error_string(st.MPI_ERROR, mesg, &mesglen);
                                    lprintf("MPI", 0, "Source [%d] Tag [%] ERROR: %s\n", st.MPI_SOURCE, st.MPI_TAG, mesg);
                                }
                                error(1, 1, "read_gauge_field " __FILE__, "Cannot receive buffer");
                            }
#endif
                        }
                    }
#endif

                    if (pid == PID) { /* copy buffer into memory */
                        int lsite[4];
                        suNg *cm;

                        /* convert global to local coordinate */
                        origin_coord(lsite);
                        lsite[0] = g[0] - lsite[0];
                        lsite[1] = g[1] - lsite[1];
                        lsite[2] = g[2] - lsite[2];

                        /* copy buffer in place */
                        cm = (suNg *)buff;
                        for (lsite[3] = 0; lsite[3] < Z; ++lsite[3]) { /* loop on local Z */
                            int ix = ipt(lsite[0], lsite[1], lsite[2], lsite[3]);
                            suNg *pm = pu_gauge(ix, 0);
                            *(pm++) = *(cm++); /* copy 4 directions */
                            *(pm++) = *(cm++);
                            *(pm++) = *(cm++);
                            *(pm) = *(cm++);

                            if (four_fermion_active) {
                                double *s_buff = (double *)cm;
                                double *s = _FIELD_AT(ff_sigma, ix);
                                double *pl = _FIELD_AT(ff_pi, ix);
                                *(s) = *(s_buff++);
                                *(pl) = *(s_buff++);
                                cm = (suNg *)s_buff;
                            }
                        }
                    }

                } /* end loop over processors in Z direction */
            }
        }
    } /* end loop over T, X and Y direction */

    if (PID == 0) { fclose(fp); }
    free(buff);

#ifdef WITH_GPU
    copy_to_gpu_suNg_field(u_gauge);
#endif

    /* start sendrecv of global gauge field */
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);

    /* check average plaquette */
    testplaq = avr_plaquette();
    if (PID == 0) {
        if (fabs(testplaq - plaq) > 1.e-12) {
            lprintf("WARNING", 0, "Stored plaquette value [%e] do not match the configuration! [diff=%e]\n", plaq,
                    fabs(testplaq - plaq));
            //error(1, 1, "read_gauge_field " __FILE__, "Plaquette value mismatch");
        }
    }

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("IO", 0, "Configuration [%s] read [%lf sec] Plaquette=%e\n", filename, elapsed_sec, testplaq);

#ifdef WITH_MPI
    MPI_Barrier(GLB_COMM);
#endif
}

void write_ranlxd_state(char filename[]) {
    FILE *fp = NULL;
    const int nproc = NP_T * NP_X * NP_Y * NP_Z;
    int p[4];
    int *buff = NULL;
    int pid = 0;
    int rsize = 0;

#ifdef WITH_MPI
    /* MPI variables */
    MPI_Group wg, cg;
    MPI_Status st;
    int cid;
    int mpiret;
    (void)mpiret; /*No warning for unused variable. */
#endif

    rsize = rlxd_size();

    if (PID == 0) {
        int d[2] = { nproc, rsize };
        error((fp = fopen(filename, "wb")) == NULL, 1, "write_ranlxd_state", "Failed to open file for writing");
        /* write nproc and rsize */
        error(fwrite_BE_int(d, (size_t)(2), fp) != (2), 1, "write_ranlxd_state", "Failed to write header");
    }

#ifdef WITH_MPI
    MPI_Comm_group(GLB_COMM, &wg);
    MPI_Comm_group(cart_comm, &cg);
#endif

    Timer clock;
    timer_set(&clock);

    buff = malloc(sizeof(*buff) * rsize);

    for (p[0] = 0; p[0] < NP_T; ++p[0]) { /* loop over processor grid */
        for (p[1] = 0; p[1] < NP_X; ++p[1]) {
            for (p[2] = 0; p[2] < NP_Y; ++p[2]) {
                for (p[3] = 0; p[3] < NP_Z; ++p[3]) {
#ifdef WITH_MPI
                    MPI_Cart_rank(cart_comm, p, &cid);
                    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
                    if (pid == PID) { /* get state */
                        rlxd_get(buff);
                    }
#ifdef WITH_MPI
                    MPI_Barrier(GLB_COMM);
                    if (pid != 0) { /* do send/receive only if the data is not on PID 0 */
                        /* send buffer */
                        if (pid == PID) {
                            mpiret = MPI_Send(buff, rsize, MPI_INT, 0, 991, GLB_COMM);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                error(1, 1, "write_ranlxd_state " __FILE__, "Cannot send buffer");
                            }
#endif
                        }
                        /* receive buffer */
                        if (PID == 0) {
                            mpiret = MPI_Recv(buff, rsize, MPI_INT, pid, 991, GLB_COMM, &st);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                if (st.MPI_ERROR != MPI_SUCCESS) {
                                    MPI_Error_string(st.MPI_ERROR, mesg, &mesglen);
                                    lprintf("MPI", 0, "Source [%d] Tag [%] ERROR: %s\n", st.MPI_SOURCE, st.MPI_TAG, mesg);
                                }
                                error(1, 1, "write_ranlxd_state " __FILE__, "Cannot receive buffer");
                            }
#endif
                        }
                    }
#endif

                    /* write buffer to file */
                    if (PID == 0) {
                        error(fwrite_BE_int(buff, (size_t)(rsize), fp) != (rsize), 1, "write_ranlxd_state",
                              "Failed to write gauge field to file");
                    }

                } /* end loop over processors in Z direction */
            }
        }
    } /* end loop over T, X and Y direction */

    if (PID == 0) { fclose(fp); }
    free(buff);

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("IO", 0, "Ranlxd state [%s] saved [%lf sec]\n", filename, elapsed_sec);

#ifdef WITH_MPI
    MPI_Barrier(GLB_COMM);
#endif
}

void read_ranlxd_state(char filename[]) {
    FILE *fp = NULL;
    const int nproc = NP_T * NP_X * NP_Y * NP_Z;
    int hproc = 0; /* number of states in the file */
    int cproc = 0;
    int p[4];
    int *buff = NULL;
    int pid = 0;
    int rsize = 0;

#ifdef WITH_MPI
    /* MPI variables */
    MPI_Group wg, cg;
    MPI_Status st;
    int cid;
    int mpiret;
    (void)mpiret; /*No warning for unused variable. */
#endif

    Timer clock;
    timer_set(&clock);

    rsize = rlxd_size();

    if (PID == 0) {
        int d[2] = { 0 };
        error((fp = fopen(filename, "rb")) == NULL, 1, "read_ranlxd_state", "Failed to open file for reading");
        /* read number of states in the file */
        error(fread_BE_int(d, (size_t)(2), fp) != (2), 1, "read_ranlxd_state", "Failed to read header");
        /* give an error if hproc != nproc */
        hproc = d[0];
        if (hproc != nproc) {
            lprintf("IO", 10, "the number of ranlxd states read [%d] doesn't match the number of processes [%d].\n", hproc,
                    nproc);
            error(1, 1, "read_ranlxd_state " __FILE__, "number of ranlxd states mismatch");
        }
        /* check if ranlxd size is the same as in the header */
        if (rsize != d[1]) {
            lprintf("IO", 0, "Read value of ranlxd size [%d] do not match this code [%d].\n", d[1], rsize);
            error(1, 1, "read_ranlxd_state " __FILE__, "ranlxd size mismatch");
        }
    }

    bcast_int(&hproc, 1);

#ifdef WITH_MPI
    MPI_Comm_group(GLB_COMM, &wg);
    MPI_Comm_group(cart_comm, &cg);
#endif

    buff = malloc(sizeof(*buff) * rsize);

    for (p[0] = 0; p[0] < NP_T; ++p[0]) { /* loop over processor grid */
        for (p[1] = 0; p[1] < NP_X; ++p[1]) {
            for (p[2] = 0; p[2] < NP_Y; ++p[2]) {
                for (p[3] = 0; p[3] < NP_Z; ++p[3]) {
                    /* check if there are more states available */
                    if (hproc > cproc) {
                        cproc++;

#ifdef WITH_MPI
                        MPI_Cart_rank(cart_comm, p, &cid);
                        MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
                        /* read buffer from file */
                        if (PID == 0) {
                            error(fread_BE_int(buff, (size_t)(rsize), fp) != (rsize), 1, "read_ranlxd_state",
                                  "Failed to read ranlxd state from file");
                        }

#ifdef WITH_MPI
                        /* MPI_Barrier(GLB_COMM); */
                        if (pid != 0) { /* do send/receive only if the data is not on PID 0 */
                            /* send buffer */
                            if (PID == 0) {
                                mpiret = MPI_Send(buff, rsize, MPI_INT, pid, 990, GLB_COMM);
#ifndef NDEBUG
                                if (mpiret != MPI_SUCCESS) {
                                    char mesg[MPI_MAX_ERROR_STRING];
                                    int mesglen;
                                    MPI_Error_string(mpiret, mesg, &mesglen);
                                    lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                    error(1, 1, "read_ranlxd_state " __FILE__, "Cannot send buffer");
                                }
#endif
                            }
                            /* receive buffer */
                            if (PID == pid) {
                                mpiret = MPI_Recv(buff, rsize, MPI_INT, 0, 990, GLB_COMM, &st);
#ifndef NDEBUG
                                if (mpiret != MPI_SUCCESS) {
                                    char mesg[MPI_MAX_ERROR_STRING];
                                    int mesglen;
                                    MPI_Error_string(mpiret, mesg, &mesglen);
                                    lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                    if (st.MPI_ERROR != MPI_SUCCESS) {
                                        MPI_Error_string(st.MPI_ERROR, mesg, &mesglen);
                                        lprintf("MPI", 0, "Source [%d] Tag [%] ERROR: %s\n", st.MPI_SOURCE, st.MPI_TAG, mesg);
                                    }
                                    error(1, 1, "read_ranlxd_state " __FILE__, "Cannot receive buffer");
                                }
#endif
                            }
                        }
#endif

                        if (pid == PID) { /* reset state */
                            rlxd_reset(buff);
                        }
                    }

                } /* end loop over processors in Z direction */
            }
        }
    } /* end loop over T, X and Y direction */

    if (PID == 0) { fclose(fp); }
    free(buff);

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("IO", 0, "Ranlxd state [%s] read [%lf sec]\n", filename, elapsed_sec);

#ifdef WITH_MPI
    MPI_Barrier(GLB_COMM);
#endif
}

void write_scalar_field(char filename[]) {
    FILE *fp = NULL;
    int g[4], p[4];
    double *buff = NULL;
    int pid = 0;
    int zsize, rz;

#ifdef WITH_MPI
    /* MPI variables */
    MPI_Group wg, cg;
    MPI_Status st;
    int cid;
    int mpiret;
    (void)mpiret;
#endif

#ifndef ALLOCATE_REPR_GAUGE_FIELD
    complete_sendrecv_suNg_scalar_field(u_scalar);
//  apply_BCs_on_represented_gauge_field(); //Save the link variables with periodic boundary conditions
#endif

    if (PID == 0) {
        int d[5] = { NG, GLB_T, GLB_X, GLB_Y, GLB_Z };
        error((fp = fopen(filename, "wb")) == NULL, 1, "write_gauge_field", "Failed to open file for writing");
        /* write NG and global size */
        error(fwrite_LE_int(d, (size_t)(5), fp) != (5), 1, "write_gauge_field", "Failed to write gauge field geometry");
    }

#ifdef WITH_MPI
    MPI_Comm_group(GLB_COMM, &wg);
    MPI_Comm_group(cart_comm, &cg);
#endif

    Timer clock;
    timer_set(&clock);

    zsize = GLB_Z / NP_Z;
    rz = GLB_Z - zsize * NP_Z;
    buff = malloc(sizeof(suNg_vector) * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));

    g[3] = 0;
    for (g[0] = 0; g[0] < GLB_T; ++g[0]) { /* loop over T, X and Y direction */
        for (g[1] = 0; g[1] < GLB_X; ++g[1]) {
            for (g[2] = 0; g[2] < GLB_Y; ++g[2]) {
#ifdef WITH_MPI
                glb_to_proc(g, p); /* get the processor coordinate */
#endif
                for (p[3] = 0; p[3] < NP_Z; ++p[3]) { /* loop over processors in Z direction */
                    int bsize;
                    bsize = sizeof(suNg_vector) / sizeof(double) *
                            (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */

#ifdef WITH_MPI
                    MPI_Cart_rank(cart_comm, p, &cid);
                    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
                    if (pid == PID) { /* fill link buffer */
                        int lsite[4];
                        suNg_vector *cm;

                        /* convert global to local coordinate */
                        origin_coord(lsite);
                        lsite[0] = g[0] - lsite[0];
                        lsite[1] = g[1] - lsite[1];
                        lsite[2] = g[2] - lsite[2];

                        /* fill buffer */
                        cm = (suNg_vector *)buff;
                        for (lsite[3] = 0; lsite[3] < Z; ++lsite[3]) { /* loop on local Z */
                            int ix = ipt(lsite[0], lsite[1], lsite[2], lsite[3]);
                            suNg_vector *pm = pu_scalar(ix);
                            *(cm++) = *(pm);
                        }
                    }
#ifdef WITH_MPI
                    MPI_Barrier(GLB_COMM);
                    if (pid != 0) { /* do send/receive only if the data is not on PID 0 */
                        /* send buffer */
                        if (pid == PID) {
#ifndef NDEBUG
                            error(Z != (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)), 1, "write_scalar_field",
                                  "Local lattice size mismatch!");
#endif
                            mpiret = MPI_Send(buff, bsize, MPI_DOUBLE, 0, 999, GLB_COMM);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                error(1, 1, "write_scalar_field " __FILE__, "Cannot send buffer");
                            }
#endif
                        }
                        /* receive buffer */
                        if (PID == 0) {
                            mpiret = MPI_Recv(buff, bsize, MPI_DOUBLE, pid, 999, GLB_COMM, &st);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                if (st.MPI_ERROR != MPI_SUCCESS) {
                                    MPI_Error_string(st.MPI_ERROR, mesg, &mesglen);
                                    lprintf("MPI", 0, "Source [%d] Tag [%] ERROR: %s\n", st.MPI_SOURCE, st.MPI_TAG, mesg);
                                }
                                error(1, 1, "write_scalar_field " __FILE__, "Cannot receive buffer");
                            }
#endif
                        }
                    }
#endif

                    /* write buffer to file */
                    if (PID == 0) {
                        // Little endian output CHANGE LATER
                        error(fwrite_LE_double(buff, (size_t)(bsize), fp) != (bsize), 1, "write_scalar_field",
                              "Failed to write scalar field to file");
                    }

                } /* end loop over processors in Z direction */
            }
        }
    } /* end loop over T, X and Y direction */

    if (PID == 0) { fclose(fp); }
    free(buff);

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("IO", 0, "Configuration [%s] saved [%lf sec]\n", filename, elapsed_sec);

#ifdef WITH_MPI
    MPI_Barrier(GLB_COMM);
#endif

#ifndef ALLOCATE_REPR_GAUGE_FIELD
    complete_sendrecv_suNg_scalar_field(u_scalar);
#endif
}

void read_scalar_field(char filename[]) {
    FILE *fp = NULL;
    int g[4], p[4];
    double *buff = NULL;
    int pid = 0;
    int zsize, rz;

#ifdef WITH_MPI
    /* MPI variables */
    MPI_Group wg, cg;
    MPI_Status st;
    int cid;
    int mpiret;
    (void)mpiret;
#endif

    Timer clock;
    timer_set(&clock);

    if (PID == 0) {
        int d[5] = { 0 }; /* contains NG,GLB_T,GLB_X,GLB_Y,GLB_Z */
        error((fp = fopen(filename, "rb")) == NULL, 1, "read_scalar_field", "Failed to open file for reading");
        /* read NG and global size */
        error(fread_LE_int(d, (size_t)(5), fp) != (5), 1, "read_scalar_field", "Failed to read scalar field geometry");
        /* Check Gauge group and Lattice dimesions */
        if (NG != d[0]) {
            lprintf("ERROR", 0, "Read value of NG [%d] do not match this code [NG=%d].\nPlease recompile.\n", d[0], NG);
            error(1, 1, "read_scalar_field " __FILE__, "Gauge group mismatch");
        }
        if (GLB_T != d[1] || GLB_X != d[2] || GLB_Y != d[3] || GLB_Z != d[4]) {
            lprintf("ERROR", 0, "Read value of global lattice size (%d,%d,%d,%d) do not match input file (%d,%d,%d,%d).\n",
                    d[1], d[2], d[3], d[4], GLB_T, GLB_X, GLB_Y, GLB_Z);
            error(1, 1, "read_scalar_field " __FILE__, "Gauge group mismatch");
        }
    }

#ifdef WITH_MPI
    MPI_Comm_group(GLB_COMM, &wg);
    MPI_Comm_group(cart_comm, &cg);
#endif

    zsize = GLB_Z / NP_Z;
    rz = GLB_Z - zsize * NP_Z;
    buff = malloc(sizeof(suNg_vector) * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));

    g[3] = 0;
    for (g[0] = 0; g[0] < GLB_T; ++g[0]) { /* loop over T, X and Y direction */
        for (g[1] = 0; g[1] < GLB_X; ++g[1]) {
            for (g[2] = 0; g[2] < GLB_Y; ++g[2]) {
#ifdef WITH_MPI
                glb_to_proc(g, p); /* get the processor coordinate */
#endif
                for (p[3] = 0; p[3] < NP_Z; ++p[3]) { /* loop over processors in Z direction */
                    int bsize;
                    bsize = sizeof(suNg_vector) / sizeof(double) *
                            (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */

#ifdef WITH_MPI
                    MPI_Cart_rank(cart_comm, p, &cid);
                    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
                    /* read buffer from file */
                    if (PID == 0) {
                        error(fread_LE_double(buff, (size_t)(bsize), fp) != (bsize), 1, "read_scalar_field",
                              "Failed to read scalar field from file");
                    }

#ifdef WITH_MPI
                    /* MPI_Barrier(GLB_COMM); */
                    if (pid != 0) { /* do send/receive only if the data is not on PID 0 */
                        /* send buffer */
                        if (PID == 0) {
                            mpiret = MPI_Send(buff, bsize, MPI_DOUBLE, pid, 998, GLB_COMM);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                error(1, 1, "read_scalar_field " __FILE__, "Cannot send buffer");
                            }
#endif
                        }
                        /* receive buffer */
                        if (PID == pid) {
#ifndef NDEBUG
                            error(Z != (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)), 1, "read_scalar_field",
                                  "Local lattice size mismatch!");
#endif
                            mpiret = MPI_Recv(buff, bsize, MPI_DOUBLE, 0, 998, GLB_COMM, &st);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                if (st.MPI_ERROR != MPI_SUCCESS) {
                                    MPI_Error_string(st.MPI_ERROR, mesg, &mesglen);
                                    lprintf("MPI", 0, "Source [%d] Tag [%] ERROR: %s\n", st.MPI_SOURCE, st.MPI_TAG, mesg);
                                }
                                error(1, 1, "read_scalar_field " __FILE__, "Cannot receive buffer");
                            }
#endif
                        }
                    }
#endif

                    if (pid == PID) { /* copy buffer into memory */
                        int lsite[4];
                        suNg_vector *cm;

                        /* convert global to local coordinate */
                        origin_coord(lsite);
                        lsite[0] = g[0] - lsite[0];
                        lsite[1] = g[1] - lsite[1];
                        lsite[2] = g[2] - lsite[2];

                        /* copy buffer in place */
                        cm = (suNg_vector *)buff;
                        for (lsite[3] = 0; lsite[3] < Z; ++lsite[3]) { /* loop on local Z */
                            int ix = ipt(lsite[0], lsite[1], lsite[2], lsite[3]);
                            suNg_vector *pm = pu_scalar(ix);
                            *(pm) = *(cm++);
                        }
                    }

                } /* end loop over processors in Z direction */
            }
        }
    } /* end loop over T, X and Y direction */

    if (PID == 0) { fclose(fp); }
    free(buff);

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("IO", 0, "Configuration [%s] read [%lf sec]\n", filename, elapsed_sec);

    /* start sendrecv of global scalar field */
    start_sendrecv_suNg_scalar_field(u_scalar);
    complete_sendrecv_suNg_scalar_field(u_scalar);

#ifdef WITH_MPI
    MPI_Barrier(GLB_COMM);
#endif
}
