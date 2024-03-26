/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica, Sofie Martins                           *
 * All rights reserved.                                                      *
 \***************************************************************************/

/*******************************************************************************
 *
 * File process_init.c
 *
 * Inizialization of geometry structures
 *
 *******************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include "memory.h"
#include "io.h"
#include "utils.h"
#include "random.h"
#include "update.h"
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>

/* setup_process
 * Assign a unique RID, PID to each process and setup
 * communications as necessary
 *
 * return codes:
 * 0 => success
 *
 * OUTPUT:
 * MPI: GLB_COMM, RID, PID, WORLD_SIZE
 */

static char input_filename[256] = "input_file";
static char output_filename[256] = "out";
static char error_filename[256] = "err";

char *get_input_filename() {
    return input_filename;
}
char *get_output_filename() {
    return output_filename;
}
char *get_error_filename() {
    return error_filename;
}

static int setup_replicas();
static void setup_random();

static int setup_level = 0;

static void read_cmdline(int argc, char **argv) {
    int option, ai = 0;

    while ((option = getopt(argc, argv, "i:o:mh")) != -1) { // get option from the getopt() method
        switch (option) {
        // For option i, r, l, print that these are options
        case 'i':
            ai = 1;
            strcpy(input_filename, optarg);
            break;
        case 'o':
            strcpy(output_filename, optarg);
            break;
        case 'm':
            print_compiling_info();
            exit(0);
        case 'h':
            lprintf("read_cmdline", 0, "Standard cmdline:\n\t-i <input file>\n\t-o <log file>\n\t-m (compilation information)");
            exit(0);
        }
    }
    if (ai != 1) { lprintf("PROCESS INIT", 1, "Using default input file [%s]\n", input_filename); }
}

void setup_gauge_fields() {
    if (setup_level != 1) {
        error(0, 1, "SETUP_GAUGE_FIELDS", "setup_process has not yet been called\n");
    } else {
        setup_level = 2;
    }

    u_gauge = alloc_suNg_field(&glattice);

#ifndef REPR_FUNDAMENTAL
    u_gauge_f = alloc_suNf_field(&glattice);
#endif

#ifndef ALLOCATE_REPR_GAUGE_FIELD
    u_gauge_f = (suNf_field *)(u_gauge);
#endif

#ifdef DPHI_FLT
    u_gauge_f_flt = alloc_suNf_field_flt(&glattice);
#endif

#ifdef WITH_SMEARING
    init_smearing(1.0, 1.0);
#endif

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    clover_init(1.0);
#endif

    reset_wrk_pointers();
}

int setup_process(int *argc, char ***argv) {
    if (setup_level != 0) {
        printf("Error: the setup_process should be the first initialization function call\n");
        exit(1);
    } else {
        setup_level = 1;
    }

    register_sighandlers();

    read_cmdline(*argc, *argv);

#ifdef WITH_MPI
    /* INIT MPI*/
    int mpiret;
    int required = MPI_THREAD_SERIALIZED;

    int provided;
    mpiret = MPI_Init_thread(argc, argv, required, &provided);
    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "setup_process " __FILE__, "MPI inizialization failed");
    }
    error(provided < MPI_THREAD_SERIALIZED, 1, "setup_process " __FILE__,
          "MPI inizialization failed, The threading support level is lesser than that demanded.\n");

    init_hr_comms();

    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_PID);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_WORLD_SIZE);
    PID = MPI_PID;
    WORLD_SIZE = MPI_WORLD_SIZE;
    RID = 0;
    MPI_Comm shared_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shared_comm);
    MPI_Comm_rank(shared_comm, &LID);
#else
    RID = 0;
    PID = 0;
    LID = 0;
    WORLD_SIZE = 1;
#endif

    /* read global variables file */
    read_input(glb_var.read, input_filename);

    setup_replicas();

    /* logger setup */
    read_input(logger_var.read, input_filename);
    logger_set_input(&logger_var);

#ifndef LOG_ALLPIDS
    if (PID != 0) {
        logger_disable();
    } /* disable logger for MPI processes != 0 */
    else
#endif
    {
        FILE *stderrp;
        char sbuf[270];
        sprintf(sbuf, ">>%s_%d", output_filename, PID);
        logger_stdout(sbuf);
        sprintf(sbuf, "%s_%d", error_filename, PID);
        stderrp = freopen(sbuf, "w", stderr);
        error(stderrp == NULL, 1, "setup_process [process_init.c]", "Cannot redirect the stderr");
    }

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
        { lprintf("OMP", 0, "Number of Threads requested = %i\n", omp_get_num_threads()); }
    }
#endif

#ifdef WITH_GPU
    read_input(gpu_var.read, input_filename);
    init_gpu(gpu_var);
#endif

    lprintf("SYSTEM", 0, "Gauge group: SU(%d)\n", NG);
    lprintf("SYSTEM", 0, "Fermion representation: dim = %d\n", NF);
    lprintf("SYSTEM", 0, "[RepID: %d][world_size: %d]\n[MPI_ID: %d][MPI_size: %d]\n", RID, WORLD_SIZE, MPI_PID, MPI_WORLD_SIZE);

    print_compiling_info_short();

    //  lprintf("MAIN",0,"Logger lelvel: %d\n",logger_getlevel(0));

    /* setup lattice geometry */
    if (geometry_init() == 1) {
        finalize_process();
        return 0;
    }

#ifdef AVX2_HIREP
    lprintf("SYSTEM", 0, "AVX2 support enabled\n");
#endif

#ifdef SIMD_VECTOR_HIREP
    lprintf("SYSTEM", 0, "SIMD VECTORIZATION support enabled\n");
#endif

    setup_random();

    return 0;
}

static void setup_random() {
    read_input(rlx_var.read, get_input_filename());

    if (strcmp(rlx_var.rlxd_start, "continue") == 0 && rlx_var.rlxd_state[0] != '\0') {
        /*load saved state*/
        lprintf("MAIN", 0, "Loading rlxd state from file [%s]\n", rlx_var.rlxd_state);
        read_ranlxd_state(rlx_var.rlxd_state);
        return;
    }

    if (rlx_var.rlxd_level == 0) {
        lprintf("SYSTEM", 0, "No rlx_level defined, skipping init of the random number generator.\n");
    } else {
        lprintf("SETUP_RANDOM", 0, "RLXD [%d,%d]\n", rlx_var.rlxd_level, rlx_var.rlxd_seed + MPI_PID);
        rlxd_init(rlx_var.rlxd_level, rlx_var.rlxd_seed);
    }
}

/* this function is intended to clean up before process ending
 *
 */
void finalize_process() {
    free_ghmc();

    free_BCs();

    free_wrk_space();

    /* free memory */
    free_suNg_field(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
    free_suNf_field(u_gauge_f);
#endif
    if (u_scalar != NULL) { free_suNg_scalar_field(u_scalar); }

    if (u_gauge_f_flt != NULL) { free_suNf_field_flt(u_gauge_f_flt); }

    // #ifndef WITH_NEW_GEOMETRY
    //   free_geometry_mpi_eo();
    // #endif

#ifdef WITH_MPI
    finalize_hr_comms();
    /* MPI variables */
    int init;
    MPI_Initialized(&init);
    if (init) { MPI_Finalize(); }
#endif
    lprintf("SYSTEM", 0, "Process finalized.\n");
}

/* setup_replicas
 * Split MPI_COMM_WORLD into replica communicators GLB_COMM
 *
 * return codes:
 * 0 => success
 *
 * AFFECTS THE GLOBAL VARIABLES: GLB_COMM, RID, PID, WORLD_SIZE
 */
static int setup_replicas() {
#ifdef WITH_MPI
    if (N_REP > 1) {
        int mpiret;

        MPI_Initialized(&mpiret);
        if (!mpiret) {
            lprintf("MPI", 0, "ERROR: MPI has not been initialized!!!\n");
            error(1, 1, "setup_replicas " __FILE__, "Cannot create replicas");
        }

        /* set up replicas */
        char sbuf[64];
        if ((MPI_WORLD_SIZE % N_REP) != 0) {
            error(1, 1, "setup_replicas " __FILE__, "MPI_WORLD_SIZE is not a multiple of the number of replicas!");
        }

        RID = MPI_PID / (MPI_WORLD_SIZE / N_REP); /* Replica ID */
        mpiret = MPI_Comm_split(MPI_COMM_WORLD, RID, 0, &GLB_COMM);
        if (mpiret != MPI_SUCCESS) {
            char mesg[MPI_MAX_ERROR_STRING];
            int mesglen;
            MPI_Error_string(mpiret, mesg, &mesglen);
            lprintf("MPI", 0, "ERROR: %s\n", mesg);
            error(1, 1, "setup_replicas " __FILE__, "Inizialization of replicas failed");
        }
        MPI_Comm_rank(GLB_COMM, &PID);
        MPI_Comm_size(GLB_COMM, &WORLD_SIZE);

        /* chdir to replica dir */
        sprintf(sbuf, "Rep_%d", RID);
        mpiret = chdir(sbuf);
    }
#endif // ifdef WITH_MPI

    return 0;
}
