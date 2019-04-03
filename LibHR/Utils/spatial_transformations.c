/***************************************************************************\
* Copyright (c)                                                             *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "global.h"
#include "utils.h"
#include "suN.h"
#include "memory.h"
#include "communications.h"
#include "logger.h"
#include "hr_complex.h"
#include "random.h"

#define gauge_field_pointer(u, ix, mu) ((u->ptr) + coord_to_index(ix, mu))

static int wrk_idx1;
static int wrk_idx2;
static suNg_field *wrk_g_1 = NULL;
static suNg_field *wrk_g_2 = NULL;
static int *wrk_iup_1;
static int *wrk_idn_1;
static int *wrk_iup_2;
static int *wrk_idn_2;
static suNg_field *tmp_g;

static int *tmp_iup;
static int *tmp_idn;

static int blocking_level;

static int *active_slices_list = NULL;
static int n_active_slices;

void initialize_spatial_active_slices(int *tlist)
{
    if (NP_X != 1 || NP_Y != 1 || NP_Z != 1)
        error(1, 1, "initialize_spatial_active_slices [spatial_transformations.c]",
              "Error the spatial transformation module can be used only with non parallel spatial boundary conditions");

    if (active_slices_list == NULL)
    {
        wrk_idx2 = reserve_wrk_space_with_pointers(&wrk_g_2, &wrk_iup_2, &wrk_idn_2);

        tmp_g = wrk_g_2;
        tmp_iup = wrk_iup_2;
        tmp_idn = wrk_idn_2;

        int t = 0;

        if (tlist != NULL)
        {
            n_active_slices = 0;
            for (t = 0; t < GLB_T; t++)
            {
                if (tlist[t] != 0 && t >= zerocoord[0] && t < zerocoord[0] + T)
                    n_active_slices++;
            }
            active_slices_list = malloc(n_active_slices * sizeof(int));
            n_active_slices = 0;

            for (t = 0; t < GLB_T; t++)
            {
                if (tlist[t] != 0 && t >= zerocoord[0] && t < zerocoord[0] + T)
                {
                    active_slices_list[n_active_slices] = t - zerocoord[0];
                    n_active_slices++;
                }
            }
            lprintf("SPATIAL TRASFORMATION", 0, "Partial Spatial lattice trasformations enabled (%d t-slices per proc)\n", n_active_slices);
        }
        else
        {
            lprintf("SPATIAL TRASFORMATION", 0, "Full Spatial lattice trasformations enabled\n");
            active_slices_list = malloc(T * sizeof(int));
            n_active_slices = T;

            for (t = 0; t < T; t++)
                active_slices_list[t] = t;
        }
    }
    else
        lprintf("SPATIAL TRASFORMATION", 0, "Already initialized\n");
}

void free_spatial_active_slices()
{
    if (active_slices_list != NULL)
    {
        free(active_slices_list);
    }
}

static void slice_spatial_blocking(int t, suNg_field *g_in, int *iup_in, int *idn_in, suNg_field *g_out, int *iup_out, int *idn_out)
{

    int ix, iy, iz, mu, mid, midn;
    suNg *u, *w, *v;
    for (ix = 0; ix < X; ix++)
        for (iy = 0; iy < Y; iy++)
            for (iz = 0; iz < Z; iz++)
            {
                mid = ipt(t, ix, iy, iz);
                for (mu = 1; mu < 4; mu++)
                {
                    u = gauge_field_pointer(g_out, mid, mu);
                    v = gauge_field_pointer(g_in, mid, mu);
                    midn = iup_in[4 * mid + mu];
                    iup_out[(mid)*4 + (mu)] = iup_in[midn * 4 + mu];
                    w = gauge_field_pointer(g_in, midn, mu);

                    _suNg_times_suNg(*u, *v, *w);

                    midn = idn_in[4 * mid + mu];
                    idn_out[4 * mid + mu] = idn_in[4 * midn + mu];
                }
            }
}

void spatial_blocking_wrkspace(eval_spat_block eval, unsigned int level)
{
    int t, idx_tmp;
    static int entry = 1;
    int last_level = 0;
    if (blocking_level == level)
        return;

    if (entry)
    {
        wrk_idx1 = reserve_wrk_space_with_pointers(&wrk_g_1, &wrk_iup_1, &wrk_idn_1);
        entry = 0;
        last_level = 1;
    }
    if (eval || blocking_level > level)
    {
        tmp_g = u_gauge;
        tmp_iup = iup;
        tmp_idn = idn;
        blocking_level = 0;
    }

    if (level != blocking_level + 1)
        spatial_blocking_wrkspace(CONT_SBLK, level - 1);

    for (t = 0; t < n_active_slices; t++)
        slice_spatial_blocking(active_slices_list[t], tmp_g, tmp_iup, tmp_idn, wrk_g_1, wrk_iup_1, wrk_idn_1);

    tmp_g = wrk_g_1;
    tmp_iup = wrk_iup_1;
    tmp_idn = wrk_idn_1;

    wrk_g_1 = wrk_g_2;
    wrk_iup_1 = wrk_iup_2;
    wrk_idn_1 = wrk_idn_2;

    wrk_g_2 = tmp_g;
    wrk_iup_2 = tmp_iup;
    wrk_idn_2 = tmp_idn;

    idx_tmp = wrk_idx1;
    wrk_idx1 = wrk_idx2;
    wrk_idx2 = idx_tmp;

    blocking_level++;

    if (last_level)
    {
        release_wrk_space(wrk_idx1);
        set_wrk_space(wrk_idx2);
        entry = 1;
    }
}

#undef gauge_field_pointer

void assign_spatial_rotated_wrkspace(int *map, int idx_wrkspace)
{

    error(NP_X != 1 || NP_Y != 1 || NP_Z != 1, 1, "generate_spatial_rotated_wrkspace [spatial_transformations.c]",
          "Error the spatial rotation module can be used on with periodic spatial boundary conditions");

    error(GLB_X != GLB_Y || GLB_X != GLB_Z, 1, "generate_spatial_rotated_wrkspace [spatial_transformations.c]",
          "Error the spatial rotation module can be used on with cubic spatial lattices");

    error(map[0] != 0, 1, "generate_spatial_rotated_wrkspace [spatial_transformations.c]",
          "Error the map transformation must send time into time (map[0]=0)");

    int *_idn, *_iup;
    suNg_field *_g;

    set_wrk_space_and_pointers(idx_wrkspace, &_g, &_iup, &_idn);
    memcpy(_iup, iup, 4 * glattice.gsize_gauge * sizeof(int));
    memcpy(_idn, idn, 4 * glattice.gsize_gauge * sizeof(int));

    int n0, n1, n2, n3;
    int nsteps = X;
    int ip, ia;
    int tmp;
    int mu;

    for (n0 = 0; n0 < T; n0++)
    {
        ip = ipt(n0, 0, 0, 0);
        ia = ip;
        for (n1 = 0; n1 < nsteps; n1++)
        {
            for (n2 = 0; n2 < nsteps; n2++)
            {
                for (n3 = 0; n3 < nsteps; n3++)
                {
                    *pu_gauge_wrk(ia, 0) = *pu_gauge(ip, 0);

                    for (mu = 1; mu < 4; mu++)
                    {
                        if (map[mu] > 0)
                            *pu_gauge_wrk(ia, map[mu]) = *pu_gauge(ip, mu);
                        else
                        {
                            tmp = idn_wrk(ia, -map[mu]);
                            _suNg_dagger(*pu_gauge_wrk(tmp, -map[mu]), *pu_gauge(ip, mu));
                        }
                    }

                    if (map[3] > 0)
                        ia = iup_wrk(ia, map[3]);
                    else
                        ia = idn_wrk(ia, -map[3]);
                    ip = iup_wrk(ip, 3);
                }
                if (map[2] > 0)
                    ia = iup_wrk(ia, map[2]);
                else
                    ia = idn_wrk(ia, -map[2]);
                ip = iup_wrk(ip, 2);
            }
            if (map[1] > 0)
                ia = iup_wrk(ia, map[1]);
            else
                ia = idn_wrk(ia, -map[1]);
            ip = iup_wrk(ip, 1);
        }
    }
}
/*
void spatial_APE_smear_wrkspace(double *smear_val)
{
    suNg staple, tr1, tr2;

    int ix, iy, iz, it, t;
    int mu, nu, mid, midpmu, midpnu, midmnu, midpmumnu;
    suNg *v;
    static suNg *storage = NULL;

    if (storage == NULL)
        storage = malloc(sizeof(suNg) * 3 * VOL3);

    for (t = 0; t < n_active_slices; t++)
    {
        it = active_slices_list[t];
        for (ix = 0; ix < X; ix++)
            for (iy = 0; iy < Y; iy++)
                for (iz = 0; iz < Z; iz++)
                {
                    mid = ipt(it, ix, iy, iz);

                    v = storage + 3 * (ix + X * (iy + Y * (iz)));

                    for (mu = 1; mu < 4; mu++)
                    {

                        _suNg_zero(*v);

                        midpmu = iup_wrk(mid, mu);

                        for (int i = 0; i < 2-1; i++)
                        {
                            nu = (mu + i) % 3 + 1;
                            midpnu = iup_wrk(mid, nu);
                            midmnu = idn_wrk(mid, nu);
                            midpmumnu = idn_wrk(midpmu, nu);

                            //Up Staple
                            _suNg_times_suNg(tr1, *pu_gauge_wrk(mid, nu), *pu_gauge_wrk(midpnu, mu));
                            _suNg_times_suNg_dagger(staple,tr1,*pu_gauge_wrk(midpmu, nu));

                            _suNg_add_assign(*v, staple);

                            //Down Staple
                            _suNg_times_suNg(tr2, *pu_gauge_wrk(midmnu, mu), *pu_gauge_wrk(midpmumnu, nu));
                            _suNg_dagger_times_suNg(staple, *pu_gauge_wrk(midmnu, nu), tr2);

                            _suNg_add_assign(*v, staple);
                        }
                        _suNg_mul_assign(*v, *smear_val);
                        v++;
                    }
                }

        for (ix = 0; ix < X; ix++)
            for (iy = 0; iy < Y; iy++)
                for (iz = 0; iz < Z; iz++)
                {
                    mid = ipt(it, ix, iy, iz);
                    v = storage + 3 * (ix + X * (iy + Y * (iz)));
                    for (mu = 1; mu < 4; mu++)
                    {
                        
                        _suNg_add_assign(*pu_gauge_wrk(mid, mu), *v);
                        project_to_suNg(pu_gauge_wrk(mid, mu));
                        v++;
                    }
                }
    }
}*/