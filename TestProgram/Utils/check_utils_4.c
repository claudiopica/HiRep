/*******************************************************************************
*
* Check of the gauge invariance of the spatial APE smearing and of the covariant gauge projection.
*
*******************************************************************************/

#include "libhr.h"

static gtransf *g;

static void random_g(void) {
    _MASTER_FOR(&glattice, ix) {
        random_suNg(_FIELD_AT(g, ix));
    }

    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);
}

static void transform_u(void) {
    _MASTER_FOR(&glattice, ix) {
        suNg v;
        for (int mu = 0; mu < 4; mu++) {
            int iy = iup(ix, mu);
            suNg *u = pu_gauge(ix, mu);
            _suNg_times_suNg_dagger(v, *u, *_FIELD_AT(g, iy));
            _suNg_times_suNg(*u, *_FIELD_AT(g, ix), v);
        }
    }

    start_sendrecv_suNg_field(u_gauge);
    represent_gauge_field();
    // smear_gauge_field();
}

static hr_complex spat_avr_0pp_wrk() {
    static hr_complex pa, tmp;
    suNg_field *_u = u_gauge_wrk();
    start_sendrecv_suNg_field(_u);

    _OMP_PRAGMA(single) {
        pa = tmp = 0.;
    }

    _PIECE_FOR(&glattice, ixp) {
        if (ixp == glattice.inner_master_pieces) {
            _OMP_PRAGMA(master)
            /* wait for gauge field to be transfered */
            complete_sendrecv_suNg_field(_u);
            _OMP_PRAGMA(barrier)
        }
        _SITE_FOR_SUM(&glattice, ixp, ix, pa) {
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

int main(int argc, char *argv[]) {
    int return_value = 0;
    hr_complex Zeropp[2];
    double dop;

    double smear_par = 0.8;

    setup_process(&argc, &argv);

    setup_gauge_fields();
    initialize_spatial_active_slices(NULL);

    /* allocate additional memory */
    g = alloc_gtransf(&glattice);

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Requesting, one workspace... ");
    reserve_wrk_space();
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Smearing the gauge field of the workspace... ");
    spatial_APE_smear_wrkspace(&smear_par, -1);
    lprintf("MAIN", 0, "done.\n\n");

    Zeropp[0] = spat_avr_0pp_wrk();

    lprintf("MAIN", 0, "Generating and applying a random gauge transf... ");
    random_g();
    transform_u();

    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Smearing the gauge field of the workspace... ");
    spatial_APE_smear_wrkspace(&smear_par, -1);
    lprintf("MAIN", 0, "done.\n\n");

    Zeropp[1] = spat_avr_0pp_wrk();

    Zeropp[0] -= Zeropp[1];

    dop = sqrt(_complex_prod_re(Zeropp[0], Zeropp[0]));

    lprintf("MAIN", 0, "Checking gauge invariance of the spatial 0pp on smeared configurations.\n ");
    lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (dop > 10.e-14) { return_value++; }

    finalize_process();
    return return_value;
}
