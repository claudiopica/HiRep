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

#if NG == 2 && !defined(WITH_QUATERNIONS)

#ifdef NDEBUG
#define MPIRET(type)
#else
#define MPIRET(type) type =
#endif

/* 
 * The conversion is based on the following convention
 * for su(2) matrices:
 * | c0 c1 |   | q0+i*q3 -q2+i*q1 |
 * | c2 c3 | = | q2+i*q1  q0-i*q3 |
 */

//#define _FAST_CONVERSION
static void translate_su2_quat(double *q, suNg *m) {
#ifdef _FAST_CONVERSION
    q[0] = creal((*m).c[0]);
    q[1] = cimag((*m).c[1]);
    q[2] = -creal((*m).c[1]);
  q[3]=cimag((*m).c[0]));
#else
    q[0] = (creal((*m).c[0]) + creal((*m).c[3])) / 2.;
    q[1] = (cimag((*m).c[1]) + cimag((*m).c[2])) / 2.;
    q[2] = (creal((*m).c[2]) - creal((*m).c[1])) / 2.;
    q[3] = (cimag((*m).c[0]) - cimag((*m).c[3])) / 2.;
#endif
}
static void translate_quat_su2(suNg *m, double *q) {
    (*m).c[0] = q[0] + I * q[3];
    (*m).c[1] = -q[2] + I * q[1];
    (*m).c[2] = q[2] + I * q[1];
    (*m).c[3] = q[0] - I * q[3];
}

void write_gauge_field_su2q(char filename[]) {
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
#ifndef NDEBUG
    int mpiret;
#endif
#endif

#ifndef ALLOCATE_REPR_GAUGE_FIELD
    complete_sendrecv_suNg_field(u_gauge);
    apply_BCs_on_represented_gauge_field(); //Save the link variables with periodic boundary conditions
#endif

#ifdef WITH_GPU
    copy_from_gpu_suNg_field(u_gauge);
#endif

    error((NG != 2), 1, "write_gauge_field_su2q", "This function cannot be called if NG!=2");

    plaq = avr_plaquette(); /* to use as a checksum in the header */
    if (PID == 0) {
        int d[5] = { NG, GLB_T, GLB_X, GLB_Y, GLB_Z };
        error((fp = fopen(filename, "wb")) == NULL, 1, "write_gauge_field_su2q", "Failed to open file for writing");
        /* write NG and global size */
        error(fwrite_BE_int(d, (size_t)(5), fp) != (5), 1, "write_gauge_field_su2q", "Failed to write gauge field geometry");
        /* write average plaquette */
        error(fwrite_BE_double(&plaq, (size_t)(1), fp) != (1), 1, "write_gauge_field_su2q",
              "Failed to write gauge field plaquette");
    }

#ifdef WITH_MPI
    MPI_Comm_group(GLB_COMM, &wg);
    MPI_Comm_group(cart_comm, &cg);
#endif

    Timer clock;
    timer_set(&clock);

    zsize = GLB_Z / NP_Z;
    rz = GLB_Z - zsize * NP_Z;
    if (four_fermion_active) {
        buff = malloc(sizeof(double) * (4 * 4 + 2) * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));
    } else {
        buff = malloc(4 * sizeof(double) * 4 * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));
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
                    if (four_fermion_active) {
                        bsize = (2 + 4 * 4) * (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */
                    } else {
                        bsize = 4 * 4 * (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0));
                    }
#ifdef WITH_MPI
                    MPI_Cart_rank(cart_comm, p, &cid);
                    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
                    if (pid == PID) { /* fill link buffer */
                        int lsite[4];
                        double *cm;

                        /* convert global to local coordinate */
                        origin_coord(lsite);
                        lsite[0] = g[0] - lsite[0];
                        lsite[1] = g[1] - lsite[1];
                        lsite[2] = g[2] - lsite[2];

                        /* fill buffer */
                        cm = buff;
                        for (lsite[3] = 0; lsite[3] < Z; ++lsite[3]) { /* loop on local Z */
                            int ix = ipt(lsite[0], lsite[1], lsite[2], lsite[3]);
                            suNg *pm = pu_gauge(ix, 0);
                            /* copy 4 directions */
                            translate_su2_quat(cm, pm);
                            cm += 4;
                            translate_su2_quat(cm, pm + 1);
                            cm += 4;
                            translate_su2_quat(cm, pm + 2);
                            cm += 4;
                            translate_su2_quat(cm, pm + 3);
                            cm += 4;

                            if (four_fermion_active) {
                                double *s_buff = (double *)cm;
                                double *s = _FIELD_AT(ff_sigma, ix);
                                double *pl = _FIELD_AT(ff_pi, ix);
                                *(s_buff++) = *s;
                                *(s_buff++) = *pl;
                                cm = s_buff;
                            }
                        }
                    }
#ifdef WITH_MPI
                    MPI_Barrier(GLB_COMM);
                    if (pid != 0) { /* do send/receive only if the data is not on PID 0 */
                        /* send buffer */
                        if (pid == PID) {
#ifndef NDEBUG
                            error(Z != (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)), 1, "write_gauge_field_su2q",
                                  "Local lattice size mismatch!");
#endif
                            MPIRET(mpiret)
                            MPI_Send(buff, bsize, MPI_DOUBLE, 0, 999, GLB_COMM);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                error(1, 1, "write_gauge_field_su2q " __FILE__, "Cannot send buffer");
                            }
#endif
                        }
                        /* receive buffer */
                        if (PID == 0) {
                            MPIRET(mpiret)
                            MPI_Recv(buff, bsize, MPI_DOUBLE, pid, 999, GLB_COMM, &st);
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
                                error(1, 1, "write_gauge_field_su2q " __FILE__, "Cannot receive buffer");
                            }
#endif
                        }
                    }
#endif

                    /* write buffer to file */
                    if (PID == 0) {
                        error(fwrite_BE_double(buff, (size_t)(bsize), fp) != (bsize), 1, "write_gauge_field_su2q",
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
    apply_BCs_on_represented_gauge_field(); //Save the link variables with periodic boundary conditions
#endif
}

void read_gauge_field_su2q(char filename[]) {
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
#ifndef NDEBUG
    int mpiret;
#endif
#endif

    error((NG != 2), 1, "read_gauge_field_su2q", "This function cannot be called if NG!=2");

    Timer clock;
    timer_set(&clock);

    if (PID == 0) {
        int d[5] = { 0 }; /* contains NG,GLB_T,GLB_X,GLB_Y,GLB_Z */
        error((fp = fopen(filename, "rb")) == NULL, 1, "read_gauge_field_su2q", "Failed to open file for reading");
        /* read NG and global size */
        error(fread_BE_int(d, (size_t)(5), fp) != (5), 1, "read_gauge_field_su2q", "Failed to read gauge field geometry");
        /* Check Gauge group and Lattice dimesions */
        if (NG != d[0]) {
            lprintf("ERROR", 0, "Read value of NG [%d] do not match this code [NG=%d].\nPlease recompile.\n", d[0], NG);
            error(1, 1, "read_gauge_field_su2q " __FILE__, "Gauge group mismatch");
        }
        if (GLB_T != d[1] || GLB_X != d[2] || GLB_Y != d[3] || GLB_Z != d[4]) {
            lprintf("ERROR", 0, "Read value of global lattice size (%d,%d,%d,%d) do not match input file (%d,%d,%d,%d).\n",
                    d[1], d[2], d[3], d[4], GLB_T, GLB_X, GLB_Y, GLB_Z);
            error(1, 1, "read_gauge_field_su2q " __FILE__, "Gauge group mismatch");
        }
        error(fread_BE_double(&plaq, (size_t)(1), fp) != (1), 1, "read_gauge_field_su2q",
              "Failed to read gauge field plaquette");
    }

#ifdef WITH_MPI
    MPI_Comm_group(GLB_COMM, &wg);
    MPI_Comm_group(cart_comm, &cg);
#endif

    zsize = GLB_Z / NP_Z;
    rz = GLB_Z - zsize * NP_Z;
    if (four_fermion_active) {
        buff = malloc(sizeof(double) * (4 * 4 + 2) * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));
    } else {
        buff = malloc(4 * 4 * sizeof(double) * (GLB_Z / NP_Z + ((rz > 0) ? 1 : 0)));
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
                    if (four_fermion_active) {
                        bsize = (2 + 4 * 4) * (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)); /* buffer size in doubles */
                    } else {
                        bsize = 4 * 4 * (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0));
                    }
#ifdef WITH_MPI
                    MPI_Cart_rank(cart_comm, p, &cid);
                    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
                    /* read buffer from file */
                    if (PID == 0) {
                        error(fread_BE_double(buff, (size_t)(bsize), fp) != (bsize), 1, "read_gauge_field_su2q",
                              "Failed to read gauge field from file");
                    }

#ifdef WITH_MPI
                    /* MPI_Barrier(GLB_COMM); */
                    if (pid != 0) { /* do send/receive only if the data is not on PID 0 */
                        /* send buffer */
                        if (PID == 0) {
                            MPIRET(mpiret)
                            MPI_Send(buff, bsize, MPI_DOUBLE, pid, 999, GLB_COMM);
#ifndef NDEBUG
                            if (mpiret != MPI_SUCCESS) {
                                char mesg[MPI_MAX_ERROR_STRING];
                                int mesglen;
                                MPI_Error_string(mpiret, mesg, &mesglen);
                                lprintf("MPI", 0, "ERROR: %s\n", mesg);
                                error(1, 1, "read_gauge_field_su2q " __FILE__, "Cannot send buffer");
                            }
#endif
                        }
                        /* receive buffer */
                        if (PID == pid) {
#ifndef NDEBUG
                            error(Z != (GLB_Z / NP_Z + ((p[3] < rz) ? 1 : 0)), 1, "read_gauge_field",
                                  "Local lattice size mismatch!");
#endif
                            MPIRET(mpiret)
                            MPI_Recv(buff, bsize, MPI_DOUBLE, 0, 999, GLB_COMM, &st);
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
                                error(1, 1, "read_gauge_field_su2q " __FILE__, "Cannot receive buffer");
                            }
#endif
                        }
                    }
#endif

                    if (pid == PID) { /* copy buffer into memory */
                        int lsite[4];
                        double *cm;

                        /* convert global to local coordinate */
                        origin_coord(lsite);
                        lsite[0] = g[0] - lsite[0];
                        lsite[1] = g[1] - lsite[1];
                        lsite[2] = g[2] - lsite[2];

                        /* copy buffer in place */
                        cm = buff;
                        for (lsite[3] = 0; lsite[3] < Z; ++lsite[3]) { /* loop on local Z */
                            int ix = ipt(lsite[0], lsite[1], lsite[2], lsite[3]);
                            suNg *pm = pu_gauge(ix, 0);
                            /* copy 4 directions */
                            translate_quat_su2(pm, cm);
                            cm += 4;
                            translate_quat_su2(pm + 1, cm);
                            cm += 4;
                            translate_quat_su2(pm + 2, cm);
                            cm += 4;
                            translate_quat_su2(pm + 3, cm);
                            cm += 4;

                            if (four_fermion_active) {
                                double *s_buff = (double *)cm;
                                double *s = _FIELD_AT(ff_sigma, ix);
                                double *pl = _FIELD_AT(ff_pi, ix);
                                *(s) = *(s_buff++);
                                *(pl) = *(s_buff++);
                                cm = s_buff;
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
            lprintf("ERROR", 0, "Stored plaquette value [%e] do not match the configuration! [diff=%e]\n", plaq,
                    fabs(testplaq - plaq));
            error(1, 1, "read_gauge_field_su2q " __FILE__, "Plaquette value mismatch");
        }
    }

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("IO", 0, "Configuration [%s] read [%lf sec] Plaquette=%e\n", filename, elapsed_sec, testplaq);
}

void read_gauge_field_su2(char filename[]) {
    FILE *fp = NULL;
    double plaq;
    int quaternions = 0;

    error((NG != 2), 1, "read_gauge_field_su2", "This function cannot be called if NG!=2");

    if (PID == 0) {
        int file_size;
        int d[5] = { 0 }; /* contains NG,GLB_T,GLB_X,GLB_Y,GLB_Z */
        error((fp = fopen(filename, "rb")) == NULL, 1, "read_gauge_field_su2", "Failed to open file for reading");
        /* Get the file size */
        fseek(fp, 0, SEEK_END);
        file_size = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        /* read NG and global size */
        error(fread_BE_int(d, (size_t)(5), fp) != (5), 1, "read_gauge_field_su2", "Failed to read gauge field geometry");
        /* Check Gauge group and Lattice dimesions */
        if (NG != d[0]) {
            lprintf("ERROR", 0, "Read value of NG [%d] do not match this code [NG=%d].\nPlease recompile.\n", d[0], NG);
            error(1, 1, "read_gauge_field_su2 " __FILE__, "Gauge group mismatch");
        }
        if (GLB_T != d[1] || GLB_X != d[2] || GLB_Y != d[3] || GLB_Z != d[4]) {
            lprintf("ERROR", 0, "Read value of global lattice size (%d,%d,%d,%d) do not match input file (%d,%d,%d,%d).\n",
                    d[1], d[2], d[3], d[4], GLB_T, GLB_X, GLB_Y, GLB_Z);
            error(1, 1, "read_gauge_field_su2 " __FILE__, "Gauge group mismatch");
        }
        error(fread_BE_double(&plaq, (size_t)(1), fp) != (1), 1, "read_gauge_field_su2",
              "Failed to read gauge field plaquette");

        double buff[8];
        error(fread_BE_double(buff, 8, fp) != 8, 1, "read_gauge_field_su2", "Failed to read gauge field from file");

        int quaternion_buffer_size = 8 * 4 * GLB_T * GLB_X * GLB_Y * GLB_Z * sizeof(double) + 5 * sizeof(int) + sizeof(double);
        if (four_fermion_active) { quaternion_buffer_size += 2 * GLB_T * GLB_X * GLB_Y * GLB_Z * sizeof(double); }

        if (fabs(buff[0] * buff[4] + buff[1] * buff[5] + buff[2] * buff[6] + buff[3] * buff[7]) < 1e-10 &&
            fabs(buff[0] * buff[5] - buff[4] * buff[1] + buff[2] * buff[7] - buff[6] * buff[3]) < 1e-10 &&
            quaternion_buffer_size == file_size) {
            quaternions = (1 == 0);
            lprintf("IO", 0, "SU2 matrix representation\n");
        } else {
            quaternions = (1 == 1);
            lprintf("IO", 0, "SU2 quaternion representation\n");
        }

        fclose(fp);
    }
#ifdef WITH_MPI
    MPI_Bcast(&quaternions, 1, MPI_INT, 0, GLB_COMM);
#endif
    if (quaternions) {
        read_gauge_field_su2q(filename);
        apply_BCs_on_fundamental_gauge_field();

#ifndef ALLOCATE_REPR_GAUGE_FIELD
        u_gauge_f = (suNf_field *)((void *)u_gauge);
        apply_BCs_on_represented_gauge_field();
#endif
    } else {
        read_gauge_field_matrix(filename);
    }
}

#endif /* NG==2 */
