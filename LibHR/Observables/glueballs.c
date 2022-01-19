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
