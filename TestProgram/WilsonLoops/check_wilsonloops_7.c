/*******************************************************************************
*  NOCOMPILE= !NG==2
*
* Check gauge covariance of the HYP smearing
*
*******************************************************************************/

#include "libhr.h"

static void random_g(suNg_field *g) {
    //   _DECLARE_INT_ITERATOR(ix);

    _MASTER_FOR(g->type, ix) {
        random_suNg(_FIELD_AT(g, ix));
    }
}

static void transform_u(suNg_field *out, suNg_field *in, suNg_field *g) {
    //_DECLARE_INT_ITERATOR(ix);
    int iy, mu;
    suNg v;

    _MASTER_FOR(&glattice, ix) {
        for (mu = 0; mu < 4; mu++) {
            iy = iup(ix, mu);
            _suNg_times_suNg_dagger(v, *_4FIELD_AT(in, ix, mu), *_FIELD_AT(g, iy));
            _suNg_times_suNg(*_4FIELD_AT(out, ix, mu), *_FIELD_AT(g, ix), v);
        }
    }
}

int main(int argc, char *argv[]) {
    int return_value = 0;

    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    setup_gauge_fields();

    suNg_field *u[5];
    u[0] = alloc_gfield(&glattice);
    u[1] = alloc_gfield(&glattice);
    u[2] = alloc_gfield(&glattice);
    u[3] = alloc_gfield(&glattice);
    u[4] = alloc_gfield(&glattice);

    suNg_field *g = alloc_gtransf(&glattice);

    random_u(u[0]);
    start_sendrecv_gfield(u[0]);
    complete_sendrecv_gfield(u[0]);

    random_g(g);
    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);

    double HYP_weight[3] = { 1., 1., 1. };

    transform_u(u[1], u[0], g);
    start_sendrecv_gfield(u[1]);
    complete_sendrecv_gfield(u[1]);
    HYP_smearing(u[2], u[1], HYP_weight);

    HYP_smearing(u[3], u[0], HYP_weight);
    start_sendrecv_gfield(u[3]);
    complete_sendrecv_gfield(u[3]);
    transform_u(u[4], u[3], g);

    double err = 0.;
    double dtmp;
    int t, x, y, z, i;
    for (t = 0; t < T; t++) {
        for (x = 0; x < X; x++) {
            for (y = 0; y < Y; y++) {
                for (z = 0; z < Z; z++) {
                    i = ipt(t, x, y, z);
                    for (int mu = 0; mu < 4; mu++) {
                        _suNg_sub_assign(*_4FIELD_AT(u[4], i, mu), *_4FIELD_AT(u[2], i, mu));
                        _suNg_sqnorm(dtmp, *_4FIELD_AT(u[4], i, mu));
                        if (dtmp > err) { err = dtmp; }
                    }
                }
            }
        }
    }

    lprintf("MAIN", 0, "Checking gauge covariance of the HYP smearing\n");
    lprintf("MAIN", 0, "\n");

#ifdef WITH_MPI
    int mpiret;

    dtmp = err / (2 * NG * NG);
    mpiret = MPI_Allreduce(&dtmp, &err, 1, MPI_DOUBLE, MPI_MAX, GLB_COMM);

    if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret, mesg, &mesglen);
        lprintf("MPI", 0, "ERROR: %s\n", mesg);
        error(1, 1, "main [check_wilsonloops_7.c]", "Cannot compute global maximum");
    }
#endif

    lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", err);
    lprintf("MAIN", 0, "(should be around 1*10^(-14) or so)\n");
    lprintf("MAIN", 0, "\n");
    if (err < 1.e-14) {
        lprintf("MAIN", 0, "check_wilsonloops_7 ... OK\n");
    } else {
        return_value += 1;
        lprintf("MAIN", 0, "check_wilsonloops_7 ... FAILED\n");
    }

    WL_free();

    free_gfield(u[0]);
    free_gfield(u[1]);
    free_gfield(u[2]);
    free_gfield(u[3]);
    free_gfield(u[4]);
    free_gtransf(g);

    finalize_process();
    return return_value;
}
