/*************************************************************************** \
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
 *
 * File update.c
 *
 * Update programs
 *
 *******************************************************************************/

#define PROJECT_INTERVAL 10

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "io.h"
#include "utils.h"
#include <math.h>

static int *dyn_gauge = NULL;

void project_gauge_field(void) {
#ifdef WITH_GPU
    exec_project();
#else
    _MASTER_FOR(&glattice, ix) {
        project_to_suNg(pu_gauge(ix, 0));
        project_to_suNg(pu_gauge(ix, 1));
        project_to_suNg(pu_gauge(ix, 2));
        project_to_suNg(pu_gauge(ix, 3));
    }
#endif

    start_sendrecv_suNg_field(u_gauge);
}

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
static void g_up_Dirichlet_BCs() {
    int ix, iy, iz, index;

    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    index = ipt(T - 1, ix, iy, iz);
                    dyn_gauge[index * 4] = dyn_gauge[index * 4 + 1] = dyn_gauge[index * 4 + 2] = dyn_gauge[index * 4 + 3] = 0;
                }
            }
        }
    }
}
#endif

#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED) || defined(BC_T_MIXED)
static void g_dn_Dirichlet_BCs() {
    int ix, iy, iz, index;

    if (COORD[0] == 0) {
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    index = ipt(0, ix, iy, iz);
                    dyn_gauge[index * 4] = dyn_gauge[index * 4 + 1] = dyn_gauge[index * 4 + 2] = dyn_gauge[index * 4 + 3] = 0;
                }
            }
        }
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    index = ipt(1, ix, iy, iz);
                    dyn_gauge[index * 4 + 1] = dyn_gauge[index * 4 + 2] = dyn_gauge[index * 4 + 3] = 0;
                }
            }
        }
    }
}
#endif

#if defined(BC_T_OPEN) || defined(BC_T_MIXED)
static void g_up_open_BCs() {
    int ix, iy, iz, index;

    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    index = ipt(T - 1, ix, iy, iz);
                    dyn_gauge[index * 4] = 0;
                }
            }
        }
    }
}
#endif

static void free_hb_boundary() {
    if (dyn_gauge != NULL) {
        free(dyn_gauge);
        dyn_gauge = NULL;
    }
}

static void init_hb_boundary() {
    dyn_gauge = malloc(sizeof(*dyn_gauge) * glattice.gsize_gauge * 4);
    atexit(&free_hb_boundary); //register cleanup function at exit

    for (int i = 0; i < glattice.gsize_gauge * 4; i++) {
        dyn_gauge[i] = 1;
    }
#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
    g_up_Dirichlet_BCs();
    g_dn_Dirichlet_BCs();
#endif
#ifdef BC_T_MIXED
    g_up_open_BCs();
    g_dn_Dirichlet_BCs();
#endif
#ifdef BC_T_OPEN
    g_up_open_BCs();
#endif
}

#ifdef WITH_NEW_GEOMETRY
#define FULL_MASK ((1u << 8) - 1)
static void update_all(double *beta, int type) {
    static int count = PROJECT_INTERVAL;

    if (count >= PROJECT_INTERVAL) {
        project_gauge_field();
        count = 0;
    }
    ++count;

    _OMP_PRAGMA(_omp_parallel) {
        suNg v;

        for (int mu = 0; mu < 4; mu++) {
            char mask = FULL_MASK;
#ifdef WITH_MPI
            _OMP_PRAGMA(master) {
                start_sendrecv_suNg_field(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            _OMP_PRAGMA(_omp_for)
            for (int j = glat_even.master_start[0]; j <= glat_even.master_end[0]; j++) {
                if (dyn_gauge[j * 4 + mu] != 0 && (imask[j] == mask)) {
                    staples(j, mu, &v);
                    cabmar(*beta, pu_gauge(j, mu), &v, type);
                }
            }
#ifdef WITH_MPI
            _OMP_PRAGMA(master) {
                complete_sendrecv_suNg_field(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            _OMP_PRAGMA(_omp_for)
            for (int j = glat_even.master_start[0]; j <= glat_even.master_end[0]; j++) {
                if (dyn_gauge[j * 4 + mu] != 0 && !(imask[j] == mask)) {
                    staples(j, mu, &v);
                    cabmar(*beta, pu_gauge(j, mu), &v, type);
                }
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            char mask = FULL_MASK;
#ifdef WITH_MPI
            _OMP_PRAGMA(master) {
                start_sendrecv_suNg_field(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            _OMP_PRAGMA(_omp_for)
            for (int j = glat_odd.master_start[0]; j <= glat_odd.master_end[0]; j++) {
                if (dyn_gauge[j * 4 + mu] != 0 && (imask[j] == mask)) {
                    staples(j, mu, &v);
                    cabmar(*beta, pu_gauge(j, mu), &v, type);
                }
            }
#ifdef WITH_MPI
            _OMP_PRAGMA(master) {
                complete_sendrecv_suNg_field(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            _OMP_PRAGMA(_omp_for)
            for (int j = glat_odd.master_start[0]; j <= glat_odd.master_end[0]; j++) {
                if (dyn_gauge[j * 4 + mu] != 0 && !(imask[j] == mask)) {
                    staples(j, mu, &v);
                    cabmar(*beta, pu_gauge(j, mu), &v, type);
                }
            }
        }
    }
}
#else
static void update_all(double *beta, int type) {
    static int count = PROJECT_INTERVAL;

    if (count >= PROJECT_INTERVAL) {
        project_gauge_field();
        count = 0;
    }
    ++count;

    _OMP_PRAGMA(_omp_parallel) {
        suNg v;

        for (int mu = 0; mu < 4; mu++) {
#ifdef WITH_MPI
            _OMP_PRAGMA(master) {
                start_sendrecv_suNg_field(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            _OMP_PRAGMA(_omp_for)
            for (int j = glat_even.master_start[0]; j <= glat_even.master_end[0]; j++) {
                if (dyn_gauge[j * 4 + mu] != 0) {
                    staples(j, mu, &v);
                    cabmar(*beta, pu_gauge(j, mu), &v, type);
                }
            }
#ifdef WITH_MPI
            _OMP_PRAGMA(master) {
                complete_sendrecv_suNg_field(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            for (int i = 1; i < glat_even.local_master_pieces; i++) {
                _OMP_PRAGMA(_omp_for)
                for (int j = glat_even.master_start[i]; j <= glat_even.master_end[i]; j++) {
                    if (dyn_gauge[j * 4 + mu] != 0) {
                        staples(j, mu, &v);
                        cabmar(*beta, pu_gauge(j, mu), &v, type);
                    }
                }
            }
        }

        for (int mu = 0; mu < 4; mu++) {
#ifdef WITH_MPI
            _OMP_PRAGMA(master) {
                start_sendrecv_suNg_field(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            _OMP_PRAGMA(_omp_for)
            for (int j = glat_odd.master_start[0]; j <= glat_odd.master_end[0]; j++) {
                if (dyn_gauge[j * 4 + mu] != 0) {
                    staples(j, mu, &v);
                    cabmar(*beta, pu_gauge(j, mu), &v, type);
                }
            }
#ifdef WITH_MPI
            _OMP_PRAGMA(master) {
                complete_sendrecv_suNg_field(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            for (int i = 1; i < glat_odd.local_master_pieces; i++) {
                _OMP_PRAGMA(_omp_for)
                for (int j = glat_odd.master_start[i]; j <= glat_odd.master_end[i]; j++) {
                    if (dyn_gauge[j * 4 + mu] != 0) {
                        staples(j, mu, &v);
                        cabmar(*beta, pu_gauge(j, mu), &v, type);
                    }
                }
            }
        }
    }
}
#endif

void update(double *beta, int nhb, int nor) {
    if (dyn_gauge == NULL) { init_hb_boundary(); }

    for (int n = 0; n < nhb; n++) {
        update_all(beta, 0);
    }

    for (int n = 0; n < nor; n++) {
        update_all(beta, 1);
    }

    start_sendrecv_suNg_field(u_gauge);
}
