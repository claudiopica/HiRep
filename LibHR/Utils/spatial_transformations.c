/***************************************************************************\
* Copyright (c)                                                             *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*
void initialize_spatial_active_slices(int *tlist)
    Initialize the list of "t" slices on which to perform the spatial transofrmations.
    if list is NULL then the code assume that all the t slices are active

void free_spatial_active_slices()
    free the list of active slices

int spatial_blocking_wrkspace(eval_spat_block eval, unsigned int level)
    Performs a spatial blocking on a new workspace up to and absolute level of blocking "level".
    The blocking can be incremental (i.e. starting from a previously blocked conf) or start from a new conf.
    The routine returns the id of the new workspace and set the pointer to the new workspace automatically.

void single_level_spatial_blocking_wrkspace(wrk_in)
    Performs a spatial blocking of the in workspace on a new workspace of one blocking "level".
    The routine returns the id of the new workspace and set the pointer to the new workspace automatically.

void assign_spatial_rotated_wrkspace(int *map,int idx_wrkspace)
    assigns to the workspace "idx_wrkspace" a spatially rotated copy of u_gauge, the rotation is defined through the map {t',x',y',z'}= Permutation{t,x,y,z}
    with the constraint t=t'.

int spatial_APE_smear_wrkspace(double *smear_val,int wrkspace_in)
    evaluate the APE smeared field of workspace "wrkspace_in" into a new wrkspace with a smear value of smear_val.
    The routine returns the id of the new workspace and set the pointer to the new workspace automatically.
*/

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

void initialize_spatial_active_slices(int *tlist)
{
    if (NP_X != 1 || NP_Y != 1 || NP_Z != 1)
        error(1, 1, "initialize_spatial_active_slices [spatial_transformations.c]",
              "Error the spatial transformation module can be used only with non parallel spatial boundary conditions");

    if (active_slices_list == NULL)
    {

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

        glbT_to_active_slices = malloc(GLB_T * sizeof(int));
        for (t = 0; t < GLB_T; t++)
            glbT_to_active_slices[t] = -1;

        for (t = 0; t < n_active_slices; t++)
            glbT_to_active_slices[active_slices_list[t]+zerocoord[0]] = t;
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

int single_level_spatial_blocking_wrkspace(int wrk_in)
{
    suNg_field *g_in, *g_out;

    int *iup_in, *idn_in, *iup_out, *idn_out;
    int wrk_out;
    if (wrk_in == -1)
    {
        g_in = u_gauge;
        iup_in = iup;
        idn_in = idn;
    }
    else
        set_wrk_space_and_pointers(wrk_in, &g_in, &iup_in, &idn_in);

    wrk_out = reserve_wrk_space_with_pointers(&g_out, &iup_out, &idn_out);
    int t;
    for (t = 0; t < n_active_slices; t++)
        slice_spatial_blocking(active_slices_list[t], g_in, iup_in, idn_in, g_out, iup_out, idn_out);

    return wrk_out;
}

int spatial_blocking_wrkspace(eval_spat_block eval, unsigned int level)
{
    int idx_tmp, t;

    static int entry = 1;
    static int init = 1;
    int last_level = 0;
    if (blocking_level == level)
        return wrk_idx2;

    if (entry)
    {
        wrk_idx1 = reserve_wrk_space_with_pointers(&wrk_g_1, &wrk_iup_1, &wrk_idn_1);
        entry = 0;
        last_level = 1;
        if (init)
        {
            wrk_idx2 = reserve_wrk_space_with_pointers(&wrk_g_2, &wrk_iup_2, &wrk_idn_2);

            tmp_g = wrk_g_2;
            tmp_iup = wrk_iup_2;
            tmp_idn = wrk_idn_2;
            init = 0;
        }
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
    return wrk_idx2;
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

int spatial_APE_smear_wrkspace(double *smear_val, int wrkspace_in)
{
    suNg staple, tr1, tr2;

    int ix, iy, iz, it, t;
    int mu, nu, mid, midpmu, midpnu, midmnu, midpmumnu;
    suNg v;
    suNg_field *gout;
    int *_iup, *_idn;
    suNg *vout;
    int wrkspace_out;

    double sm1 = 1. - *smear_val;

    double sm2 = *smear_val / (4.0);

    wrkspace_out = reserve_wrk_space_with_pointers(&gout, &_iup, &_idn);

    if (wrkspace_in == -1)
        reset_wrk_pointers();
    else
        set_wrk_space(wrkspace_in);

    for (t = 0; t < n_active_slices; t++)
    {
        it = active_slices_list[t];
        for (ix = 0; ix < X; ix++)
            for (iy = 0; iy < Y; iy++)
                for (iz = 0; iz < Z; iz++)
                {
                    mid = ipt(it, ix, iy, iz);

                    for (mu = 1; mu < 4; mu++)
                    {

                        _suNg_zero(v);
                        vout = (gout->ptr) + coord_to_index(mid, mu);
                        midpmu = iup_wrk(mid, mu);

                        for (int i = 0; i < 2; i++)
                        {
                            nu = (mu + i) % 3 + 1;
                            midpnu = iup_wrk(mid, nu);
                            midmnu = idn_wrk(mid, nu);
                            midpmumnu = idn_wrk(midpmu, nu);

                            //Up Staple
                            _suNg_times_suNg(tr1, *pu_gauge_wrk(mid, nu), *pu_gauge_wrk(midpnu, mu));
                            _suNg_times_suNg_dagger(staple, tr1, *pu_gauge_wrk(midpmu, nu));

                            _suNg_add_assign(v, staple);

                            //Down Staple
                            _suNg_times_suNg(tr2, *pu_gauge_wrk(midmnu, mu), *pu_gauge_wrk(midpmumnu, nu));
                            _suNg_dagger_times_suNg(staple, *pu_gauge_wrk(midmnu, nu), tr2);

                            _suNg_add_assign(v, staple);
                        }

                        _suNg_mul_add(*vout, sm1, *pu_gauge_wrk(mid, mu), sm2, v);

                        covariant_project_to_suNg(vout);
                    }
                }
    }
    set_wrk_space(wrkspace_out);
    return wrkspace_out;
}