/*************************************************************************** \
 * Copyright (c)                                  *   
 * All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
 *
 * File update_hb_multilevel.c
 *
 * Update programs
 *
 *******************************************************************************/

#define PROJECT_INTERVAL 10

#include "suN.h"
#include "utils.h"
#include "global.h"
#include "update.h"
#include "communications.h"
#include "logger.h"
#include "glueballs.h"

#define PI 3.141592653589793238462643383279502884197
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int *dyn_gauge = NULL;
static int max_mh_level;

#if defined(BASIC_SF) || defined(ROTATED_SF)
static void g_up_Dirichlet_BCs()
{
    int ix, iy, iz, index, lev;

    if (COORD[0] == NP_T - 1)
    {
        for (ix = 0; ix < X; ++ix)
            for (iy = 0; iy < Y; ++iy)
                for (iz = 0; iz < Z; ++iz)
                {
                    index = ipt(T - 1, ix, iy, iz);
                    for (lev = 0; lev < max_mh_level; lev++)
                    {
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4] = 0;
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 1] = 0;
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 2] = 0;
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 3] = 0;
                    }
                }
    }
}
#endif

#if defined(BASIC_SF) || defined(ROTATED_SF) || defined(BC_T_MIXED)
static void g_dn_Dirichlet_BCs()
{
    int ix, iy, iz, index, lev;

    if (COORD[0] == 0)
    {
        for (ix = 0; ix < X; ++ix)
            for (iy = 0; iy < Y; ++iy)
                for (iz = 0; iz < Z; ++iz)
                {
                    index = ipt(0, ix, iy, iz);
                    for (lev = 0; lev < max_mh_level; lev++)
                    {
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4] = 0;
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 1] = 0;
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 2] = 0;
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 3] = 0;
                    }
                }
        for (ix = 0; ix < X; ++ix)
            for (iy = 0; iy < Y; ++iy)
                for (iz = 0; iz < Z; ++iz)
                {
                    index = ipt(1, ix, iy, iz);
                    for (lev = 0; lev < max_mh_level; lev++)
                    {
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 1] = 0;
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 2] = 0;
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 3] = 0;
                    }
                }
    }
}
#endif

#if defined(BC_T_OPEN) || defined(BC_T_MIXED)
static void g_up_open_BCs()
{
    int ix, iy, iz, index, lev;

    if (COORD[0] == NP_T - 1)
    {
        for (ix = 0; ix < X; ++ix)
            for (iy = 0; iy < Y; ++iy)
                for (iz = 0; iz < Z; ++iz)
                {
                    index = ipt(T - 1, ix, iy, iz);
                    for (lev = 0; lev < max_mh_level; lev++)
                    {
                        dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4] = 0;
                    }
                }
    }
}
#endif

static void free_hb_boundary()
{
    if (dyn_gauge != NULL)
    {
        free(dyn_gauge);
        dyn_gauge = NULL;
    }
}

static void init_hb_multihit_boundary()
{
    dyn_gauge = malloc(sizeof(*dyn_gauge) * glattice.gsize_gauge * 4 * (max_mh_level));
    atexit(&free_hb_boundary); //register cleanup function at exit

    for (int i = 0; i < glattice.gsize_gauge * 4 * (max_mh_level); i++)
        dyn_gauge[i] = 1;

#if defined(BASIC_SF) || defined(ROTATED_SF)
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

    int ix, iy, iz, index, lev, it;
    for (lev = 0; lev < max_mh_level; lev++)
        for (it = 0; it < T; ++it)
        {

            if ((it + zerocoord[0] + 1) % (GLB_T / (1 << (lev + 1))) == 0 && it + zerocoord[0] + 1 != GLB_T)
            {
                for (ix = 0; ix < X; ++ix)
                    for (iy = 0; iy < Y; ++iy)
                        for (iz = 0; iz < Z; ++iz)
                        {
                            index = ipt(it, ix, iy, iz);
                            dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 1] = 0;
                            dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 2] = 0;
                            dyn_gauge[lev * (glattice.gsize_gauge * 4) + index * 4 + 3] = 0;
                        }
            }
        }
}

static void update_mh_all(int lev, double *beta, int type)
{
    static int count = PROJECT_INTERVAL;
    static int *loc_dyn;
    loc_dyn = dyn_gauge + lev * (glattice.gsize_gauge * 4);

    if (count >= PROJECT_INTERVAL)
    {

        _MASTER_FOR(&glattice, ix)
        {
            for (int mu = 0; mu < 4; mu++)
            {
                if (loc_dyn[ix * 4 + mu] != 0)
                    project_to_suNg(pu_gauge(ix, mu));
            }
            count = 0;
        }
    }
    ++count;

    _OMP_PRAGMA(_omp_parallel)
    {
        suNg v;

        for (int mu = 0; mu < 4; mu++)
        {
#ifdef WITH_MPI
            _OMP_PRAGMA(master)
            {
                start_gf_sendrecv(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            _OMP_PRAGMA(_omp_for)
            for (int j = glat_even.master_start[0]; j <= glat_even.master_end[0]; j++)
            {
                if (loc_dyn[j * 4 + mu] != 0)
                {
                    staples(j, mu, &v);
                    cabmar(*beta, pu_gauge(j, mu), &v, type);
                }
            }
#ifdef WITH_MPI
            _OMP_PRAGMA(master)
            {
                complete_gf_sendrecv(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            for (int i = 1; i < glat_even.local_master_pieces; i++)
            {
                _OMP_PRAGMA(_omp_for)
                for (int j = glat_even.master_start[i]; j <= glat_even.master_end[i]; j++)
                {
                    if (loc_dyn[j * 4 + mu] != 0)
                    {
                        staples(j, mu, &v);
                        cabmar(*beta, pu_gauge(j, mu), &v, type);
                    }
                }
            }
        }

        for (int mu = 0; mu < 4; mu++)
        {
#ifdef WITH_MPI
            _OMP_PRAGMA(master)
            {
                start_gf_sendrecv(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            _OMP_PRAGMA(_omp_for)
            for (int j = glat_odd.master_start[0]; j <= glat_odd.master_end[0]; j++)
            {
                if (loc_dyn[j * 4 + mu] != 0)
                {
                    staples(j, mu, &v);
                    cabmar(*beta, pu_gauge(j, mu), &v, type);
                }
            }
#ifdef WITH_MPI
            _OMP_PRAGMA(master)
            {
                complete_gf_sendrecv(u_gauge);
            }
            _OMP_PRAGMA(barrier)
#endif
            for (int i = 1; i < glat_odd.local_master_pieces; i++)
            {
                _OMP_PRAGMA(_omp_for)
                for (int j = glat_odd.master_start[i]; j <= glat_odd.master_end[i]; j++)
                {
                    if (loc_dyn[j * 4 + mu] != 0)
                    {
                        staples(j, mu, &v);
                        cabmar(*beta, pu_gauge(j, mu), &v, type);
                    }
                }
            }
        }
    }
}

static void update_mh(int lev, double *beta, int nhb, int nor)
{

    for (int n = 0; n < nhb; n++)
    {
        update_mh_all(lev, beta, 0);
    }

    for (int n = 0; n < nor; n++)
    {
        update_mh_all(lev, beta, 1);
    }

    start_gf_sendrecv(u_gauge);
}

void set_max_mh_level(int lev)
{
    max_mh_level = lev;
}

void update_hb_multilevel_gb_measure(int lev, double *beta, int nhb, int nor, int *ml_up, int *ml_skip, int nblockingstart, int nblockingend, double *smear_val, cor_list *lcor)
{
    int i, j;
    static double complex *one_point_gb;
    static long double norm = 1.0;
    struct timeval start, end, etime;
    int nblocking = nblockingend - nblockingstart +1;

    if (lev == 0)
    {
        if (dyn_gauge == NULL)
        {
            init_hb_multihit_boundary();
            one_point_gb = malloc(sizeof(double complex) * total_n_glue_op * nblocking * n_active_slices);
            for (i = 0; i < max_mh_level; i++)
                norm *= ml_up[i];
            norm *= GLB_VOL3 * NG;
        }
        gettimeofday(&start, 0);

        memset(one_point_gb, 0, sizeof(double complex) * total_n_glue_op * nblocking * n_active_slices);
    }

    if (lev < max_mh_level - 1)
    {
        for (i = 0; i < ml_up[lev]; i++)
        {
            for (j = 0; j < ml_skip[lev]; j++)
                update_mh(lev, beta, nhb, nor);

            update_hb_multilevel_gb_measure(lev + 1, beta, nhb, nor, ml_up, ml_skip, nblockingstart, nblockingend, smear_val, lcor);
        }
    }
    else
    {
        for (i = 0; i < ml_up[lev]; i++)
        {
            for (j = 0; j < ml_skip[lev]; j++)
                update_mh(lev, beta, nhb, nor);

            measure_1pt_glueballs(nblockingstart, nblockingend, smear_val, one_point_gb);
        }
    }

    if (lev == 0)
    {
        gettimeofday(&end, 0);

        timeval_subtract(&etime, &end, &start);
        lprintf("HB MULTILEVEL", 0, "Update and 1pt measure done [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

        for (i = 0; i < nblocking * n_active_slices * total_n_glue_op; i++)
            one_point_gb[i] /= norm;

        evaluate_1pt_functions(lcor, nblocking, one_point_gb);
        gettimeofday(&start, 0);
        timeval_subtract(&etime, &start, &end);

        lprintf("HB MULTILEVEL", 0, "1pt writing done [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);
    }
}
