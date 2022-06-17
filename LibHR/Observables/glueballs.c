/***************************************************************************\
* Copyright (c) 2019, Antonio Rago                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#include "utils.h"
#include "glueballs.h"
#include "global.h"

void measure_1pt_glueballs(int nblockingstart, int nblockingend, double *smear_val, double complex *gb_storage)
{
    int i, nt;
    int wrk1, wrk2 = -1;
    double complex *point_gb;

    point_gb = gb_storage;

    for (i = 0; i < nblockingstart; i++)
    {

        wrk1 = spatial_APE_smear_wrkspace(smear_val, wrk2);

        release_wrk_space(wrk2);

        wrk2 = single_level_spatial_blocking_wrkspace(wrk1);

        release_wrk_space(wrk1);
    }

    for (i = 0; i < nblockingend - nblockingstart; i++)
    {
        point_gb = gb_storage + i * total_n_glue_op;
        for (nt = 0; nt < n_active_slices; nt++)
        {
            eval_all_glueball_ops(active_slices_list[nt], point_gb);
            point_gb += total_n_glue_op * (nblockingend - nblockingstart + 1);
        }

        wrk1 = spatial_APE_smear_wrkspace(smear_val, wrk2);

        release_wrk_space(wrk2);

        wrk2 = single_level_spatial_blocking_wrkspace(wrk1);

        release_wrk_space(wrk1);
    }

    point_gb = gb_storage + (nblockingend - nblockingstart) * total_n_glue_op;

    for (nt = 0; nt < n_active_slices; nt++)
    {
        eval_all_glueball_ops(active_slices_list[nt], point_gb);
        point_gb += total_n_glue_op * (nblockingend - nblockingstart + 1);
    }

    release_wrk_space(wrk2);
}

void measure_1pt_torellons(double *smear_val, double complex *tor_storage)
{
    int nt;
    int wrk1;
    double complex *point_tor;

    point_tor = tor_storage;

    wrk1 = spatial_APE_smear_wrkspace(smear_val, -1);

    for (nt = 0; nt < n_active_slices; nt++)
    {
        eval_all_torellon_ops(active_slices_list[nt], point_tor);
        point_tor += total_n_tor_op;
    }

    release_wrk_space(wrk1);
}

wilson_lines *polyleg(int ix, int d)
{

    int in = ix, in1 = 0;
    suNg *w1, *w2, *w3, s1, s2;

    w2 = &s1;
    w3 = &s2;
    static wilson_lines *poly = NULL;
    static int *LL;

    if (poly == NULL)
    {
        poly = (wilson_lines *)malloc(3 * sizeof(wilson_lines));
        poly[0].p = (suNg *)malloc(npoly_dist * sizeof(suNg));
        poly[1].p = (suNg *)malloc(npoly_dist * sizeof(suNg));
        poly[2].p = (suNg *)malloc(npoly_dist * sizeof(suNg));
        poly[0].ix = poly[1].ix = poly[2].ix = -1;
        LL = (int *)malloc(4 * sizeof(int));
        LL[0] = -T;
        LL[1] = X;
        LL[2] = Y;
        LL[3] = Z;
    }

    if (ix == poly[d - 1].ix)
        return poly + d - 1;

    if (idn_wrk(ix, d) != poly[d - 1].ix)
    {
        for (int i = 0; i < npoly_dist; i++)
            in = iup_wrk(in, d);

        in1 = in;
        _suNg_unit(*w2);

        for (int i = npoly_dist; i < LL[d] - 1; i++)
        {

            w1 = pu_gauge_wrk(in, d);
            _suNg_times_suNg(*w3, *w2, *w1);
            w1 = w2;
            w2 = w3;
            w3 = w1;
            in = iup_wrk(in, d);
        }

        if (npoly_dist < LL[d])
        {
            w1 = pu_gauge_wrk(in, d);
            _suNg_times_suNg(poly[d - 1].p[npoly_dist - 1], *w2, *w1);
        }
        else
        {
            _suNg_unit(poly[d - 1].p[npoly_dist - 1]);
        }

        in = in1;
        for (int i = npoly_dist - 2; i >= 0; i--)
        {
            in = idn_wrk(in, d);
            w1 = pu_gauge_wrk(in, d);
            _suNg_times_suNg(poly[d - 1].p[i], *w1, poly[d - 1].p[i + 1]);
        }

        in = idn_wrk(in, d);
        w1 = pu_gauge_wrk(in, d);
        _suNg_times_suNg(*w3, *w1, poly[d - 1].p[0]);
        _suNg_trace(poly[d - 1].tr, *w3);
    }
    else
    {
        in1 = idn_wrk(ix, d);
        in = ix;

        w1 = pu_gauge_wrk(in1, d);
        w2 = pu_gauge_wrk(in, d);

        _suNg_times_suNg(s2, poly[d - 1].p[0], *w1);
        _suNg_dagger_times_suNg(poly[d - 1].p[0], *w2, s2);

        for (int i = 1; i < npoly_dist; i++)
        {
            in = iup_wrk(in, d);
            w2 = pu_gauge_wrk(in, d);
            _suNg_dagger_times_suNg(poly[d - 1].p[i], *w2, poly[d - 1].p[i - 1]);
        }
    }
    poly[d - 1].ix = ix;

    return poly + d - 1;
}
