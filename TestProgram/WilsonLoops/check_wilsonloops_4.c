/*******************************************************************************
*
* Check Wilson loops(t=0) = tr(id)
*
*******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int return_value = 0;

    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    setup_gauge_fields();

    suNg_field *h_gauge = alloc_suNg_field(&glattice);
    suNg *poly = amalloc(sizeof(suNg) * X * Y * Z, ALIGN);

    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);

    WL_initialize();

    int c[3] = { 2, 3, 0 };
    WL_load_path(c, 1);

    WL_Hamiltonian_gauge(h_gauge, u_gauge);

    WL_broadcast_polyakov(poly, h_gauge);

    double **WL;
    WL = amalloc(sizeof(double *), ALIGN);
    WL[0] = amalloc(sizeof(double) * GLB_T, ALIGN);

    double err = 0.;
    for (int p = 0; p < WL_path[0].nperms; p++) {
        int sign[3] = { 1, 1, 1 };

        WL_correlators(WL, h_gauge, poly, WL_path[0].nsteps, WL_path[0].path, WL_path[0].length, WL_path[0].perm[p], sign);

        double dtmp = (WL[0][0] - NG) * (WL[0][0] - NG);
        if (dtmp > err) { err = dtmp; }
    }

    lprintf("MAIN", 0, "Checking that Wilson loops(t=0) = tr(id)\n");
    lprintf("MAIN", 0, "\n");

    lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", err);
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
    lprintf("MAIN", 0, "\n");
    if (err < 1.e-15) {
        lprintf("MAIN", 0, "check_wilsonloops_4 ... OK\n");
    } else {
        return_value += 1;
        lprintf("MAIN", 0, "check_wilsonloops_4 ... FAILED\n");
    }
    WL_free();

    afree(WL[0]);
    afree(WL);
    free_suNg_field(u_gauge);
    free_suNg_field(h_gauge);
    afree(poly);

    finalize_process();
    return return_value;
}
