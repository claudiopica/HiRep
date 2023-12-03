/******************************************************************************
*
* NOCOMPILE= UPDATE_EO
* NOCOMPILE= !DPHI_FLT
*
* Test of modules
*
******************************************************************************/

#include "libhr.h"

int nhb, nor, nit, nth, nms, level, seed;
double beta;

int main(int argc, char *argv[]) {
    int return_value = 0;
    double tau;
    spinor_field *s1, *s2, *s3;

    g5QMR_fltacc_par par;
    mshift_par mpar;

    int cgiters;

    logger_map("DEBUG", "debug");

    setup_process(&argc, &argv);

    setup_gauge_fields();
    u_gauge_f_flt = alloc_suNf_field_flt(&glattice);

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    lprintf("MAIN", 0, "done.\n");

    random_u(u_gauge);

    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);

    represent_gauge_field();
    assign_ud2u_f();

    s1 = alloc_spinor_field(3, &glattice);
    s2 = s1 + 1;
    s3 = s2 + 1;

    gaussian_spinor_field(s1);

    /* TEST g5QMR_M */

    par.max_iter = 0;
    par.err2 = 1e-14;
    par.max_iter_flt = 0;
    par.err2_flt = 1e-6;

    set_dirac_mass(0.);

    lprintf("QMR TEST", 0, "\n");
    lprintf("QMR TEST", 0, "Testing g5QMR with single-precision acceleration\n");
    lprintf("QMR TEST", 0, "------------------------\n");

    cgiters = g5QMR_fltacc(&par, &D, &D_flt, s1, s3);
    lprintf("QMR TEST", 0, "Converged in %d iterations\n", cgiters);

    D(s2, s3);
    sub_assign_spinor_field(s2, s1);
    tau = sqnorm_spinor_field(s2) / sqnorm_spinor_field(s1);
    lprintf("QMR TEST", 0, "Res = %e\n", tau);

    if (tau > par.err2) { return_value += 1; }

    mpar.max_iter = 0;
    mpar.err2 = 1e-14;
    mpar.n = 1;
    double shift = 0;
    mpar.shift = &shift;

    lprintf("QMR TEST", 0, "\n");
    lprintf("QMR TEST", 0, "Testing g5QMR multishift\n");
    lprintf("QMR TEST", 0, "------------------------\n");

    zero_spinor_field(s3);
    cgiters = g5QMR_mshift(&mpar, &D, s1, s3);
    lprintf("QMR TEST", 0, "Converged in %d iterations\n", cgiters);

    D(s2, s3);
    sub_assign_spinor_field(s2, s1);
    tau = sqnorm_spinor_field(s2) / sqnorm_spinor_field(s1);
    lprintf("QMR TEST", 0, "Res = %e\n", tau);

    if (tau > par.err2) { return_value += 1; }

    free_spinor_field(s1);
    finalize_process();

    return return_value;
}
