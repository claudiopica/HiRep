/*******************************************************************************
*
* Check of the generate rotated workspace routine
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "utils.h"
#include "update.h"
#include "observables.h"
#include "random.h"
#include "logger.h"
#include "communications.h"
#include "representation.h"
#include "glueballs.h"

#include "setup.h"

static double complex spat_avr_0pp_wrk()
{
    static double complex pa, tmp;
    suNg_field *_u = u_gauge_wrk();
    start_gf_sendrecv(_u);

    _OMP_PRAGMA(single)
    {
        pa = tmp = 0.;
    }

    _PIECE_FOR(&glattice, ixp)
    {
        if (ixp == glattice.inner_master_pieces)
        {
            _OMP_PRAGMA(master)
            /* wait for gauge field to be transfered */
            complete_gf_sendrecv(_u);
            _OMP_PRAGMA(barrier)
        }
        _SITE_FOR_SUM(&glattice, ixp, ix, pa)
        {
            cplaq_wrk(&tmp, ix, 1, 2);
            pa += tmp;
            cplaq_wrk(&tmp, ix, 2, 1);
            pa += tmp;
            cplaq_wrk(&tmp, ix, 1, 3);
            pa += tmp;
            cplaq_wrk(&tmp, ix, 3, 1);
            pa += tmp;
            cplaq_wrk(&tmp, ix, 3, 2);
            pa += tmp;
            cplaq_wrk(&tmp, ix, 2, 3);
            pa += tmp;
        }
    }

    global_sum((double *)(&pa), 2);

#ifdef BC_T_OPEN
    pa /= 6.0 * NG * GLB_VOLUME * (GLB_T - 1) / GLB_T;
#else
    pa /= 6.0 * NG * GLB_VOLUME;
#endif

    return pa;
}

int main(int argc, char *argv[])
{

    int return_value = 0;
    int idx_wrk;
    double complex plaq[2], test;
    double dop;
    double complex res;
    int mu, j;
    int **space_rotations;

    setup_process(&argc, &argv);

    setup_gauge_fields();

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_gf_sendrecv(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");

    plaq[0] = spat_avr_0pp_wrk();

    lprintf("MAIN", 0, "Requesting, one workspace gauge field\n");
    idx_wrk=reserve_wrk_space();
    lprintf("MAIN", 0, "done.\n\n");

    space_rotations = direct_spatial_rotations();

    lprintf("MAIN", 0, "Resetting and initializing (rotating) workspace gauge field once for each of the 48 cubic rotations\n\n");
    for (j = 0; j < 48; j++)
    {
        lprintf("MAIN", 0, "Rotation %d  %d->%d   %d->%d   %d->%d   %d->%d\n", j, 0, space_rotations[j][0], 1, space_rotations[j][1], 2, space_rotations[j][2], 3, space_rotations[j][3]);

        _MASTER_FOR(&glattice, i)
        {
            for (mu = 0; mu < 4; mu++)
            {
                _suNg_zero(*pu_gauge_wrk(i, mu));
            }
        }
        assign_spatial_rotated_wrkspace(space_rotations[j], idx_wrk);

        lprintf("MAIN", 0, "Checking if the workspace field is assigned to a SU(N) element on all the links\n");
        test = 0;
        _MASTER_FOR(&glattice, i)
        {
            for (mu = 0; mu < 4; mu++)
            {
                det_suNg(&res, pu_gauge_wrk(i, mu));
                test += (res - 1.0) / (4 * GLB_VOLUME);
            }
        }
        global_sum((double *)(&test), 2);
        dop = sqrt(_complex_prod_re(test, test));
        lprintf("MAIN", 0, "Maximal normalized difference = %.4e \n", dop);
        lprintf("MAIN", 0, "(should be smaller 1*10^(-15) or so)\n");
        if (dop > 10.e-14)
            return_value++;

        plaq[1] = spat_avr_0pp_wrk();

        test = plaq[0] - plaq[1];
        dop = sqrt(_complex_prod_re(test, test));
        lprintf("MAIN", 0, "Evaluating the difference of a 0++ op on the rotated and original configurations\n");

        lprintf("MAIN", 0, "Maximal normalized difference = %.4e \n", dop);
        lprintf("MAIN", 0, "(should be smaller 1*10^(-15) or so)\n");
        if (dop > 10.e-14)
            return_value++;
        lprintf("MAIN", 0, "done.\n\n");
    }
    finalize_process();
    return return_value;
}
