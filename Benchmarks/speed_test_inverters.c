/*******************************************************************************
*
* Speed test of available inverters on existing configuration
*
*******************************************************************************/

#include "libhr.h"

int nhb, nor, nit, nth, nms, level, seed;
double beta;
spinor_field *tmp;
spinor_field_flt *tmp_flt;
static double hmass = 0.859296;

int main(int argc, char *argv[]) {
    mshift_par par;
    Timer clock;
    int cgiters;
    double elapsed;
    int n_dirac;
    double tau;

    spinor_field *s1, *s2, *tmp;
    spinor_field *res;

    setup_process(&argc, &argv);
    setup_random_gauge_fields();
    //setup_gauge_fields();

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    double csw_input = 1.905900;
    set_csw(&csw_input);
#endif

    //read_gauge_field(test_cnfg_filename);
    //represent_gauge_field();
    set_dirac_mass(-hmass);

    lprintf("SANITY CHECK", 0, "Initial plaquette: %1.8e\n", avr_plaquette());

    par.n = 1;
    par.shift = (double *)malloc(sizeof(double) * (par.n));
    par.err2 = 1.e-16;
    par.max_iter = 0;
#ifndef UPDATE_EO
    res = alloc_spinor_field(par.n + 3, &glattice);
#else
    res = alloc_spinor_field(par.n + 3, &glat_even);
#endif
    s1 = res + par.n;
    s2 = s1 + 1;
    tmp = s2 + 1;
    par.shift[0] = 0;

#ifndef UPDATE_EO
    spinor_field *check = alloc_spinor_field(1, &glattice);
#else
    spinor_field *check = alloc_spinor_field(1, &glat_even);
#endif
    gaussian_spinor_field(s1);
    copy(check, s1); // Save initial field for later
    zero(s2);

#ifndef UPDATE_EO
    geometry_descriptor *geometry = &glattice;
#else
    geometry_descriptor *geometry = &glat_even;
#endif

    timer_set(&clock);

    /*************  CG MSHIFT  ***************/
    lprintf("CG_MSHIFT TEST", 0, "---------------------\n");
    lprintf("CG_MSHIFT TEST", 0, "Testing CG multi-shift\n");
    lprintf("CG_MSHIFT TEST", 0, "---------------------\n");
    copy(s1, check);
    zero(s2);
    zero(res);

    n_dirac = getMVM();
    timer_lap(&clock);
    cgiters = cg_mshift(&par, &H2, s1, res);
    elapsed = timer_lap(&clock) * 1.e-3;
    n_dirac = getMVM();

    lprintf("CG_MSHIFT TEST", 0, "CG Mshift inversion took %lf msec\n", elapsed);
    lprintf("CG_MSHIFT TEST", 0, "Inverter needed %d hopping term applications.\n", n_dirac);
    lprintf("CG_MSHIFT TEST", 0, "Last CG mshift converged in = %d iterations\n", cgiters);

    /* check precision reached */
    H2(s2, res);
    sub_assign(s2, s1);
    tau = sqnorm(s2) / sqnorm(s1);
    lprintf("CG_MSHIFT TEST", 0, "test cg_mshift = %e (given relative inverter precision: %e)\n", tau, par.err2);
    lprintf("CG_MSHIFT TEST", 0, "Done\n");

    /*************  g5QMR  ***************/
    lprintf("g5QMR TEST", 0, "---------------------\n");
    lprintf("g5QMR TEST", 0, "Testrun g5QMR\n");
    lprintf("g5QMR TEST", 0, "---------------------\n");
    copy(s1, check);
    zero(s2);
    zero(res);

    /* Invert and print results */
    n_dirac = getMVM();
    timer_lap(&clock);
    cgiters = g5QMR_mshift(&par, &H2, s1, res);
    elapsed = timer_lap(&clock) * 1.e-3;
    n_dirac = getMVM();

    lprintf("g5QMR TEST", 0, "g5QMR needed %lf msec\n", elapsed);
    lprintf("g5QMR TEST", 0, "g5QMR converged in = %d steps\n", cgiters);
    lprintf("g5QMR TEST", 0, "Inverter needed %d hopping term applications.\n", n_dirac);

    /* check precision reached */
    H2(s2, res);
    sub_assign(s2, s1);
    tau = sqnorm(s2) / sqnorm(s1);
    lprintf("g5QMR TEST", 0, "test g5QMR = %e (given relative inverter precision: %e)\n", tau, par.err2);
    lprintf("g5QMR TEST", 0, "Done\n");

    /*************  BiCGstab  ***************/
    lprintf("BICGSTAB TEST", 0, "---------------------\n");
    lprintf("BICGSTAB TEST", 0, "Testing BiCGstab\n");
    lprintf("BICGSTAB TEST", 0, "---------------------\n");
    copy(s1, check);
    zero(s2);
    zero(res);

    n_dirac = getMVM();
    timer_lap(&clock);
    g5(tmp, s1);
    cgiters = BiCGstab(&par, &D, tmp, res);
    g5(tmp, res);
    cgiters += BiCGstab(&par, &D, tmp, res);
    elapsed = timer_lap(&clock) * 1.e-3;
    n_dirac = getMVM();

    lprintf("BICGSTAB TEST", 0, "BiCGstab inversion took %lf msec\n", elapsed);
    lprintf("BICGSTAB TEST", 0, "Inverter needed %d hopping term applications.\n", n_dirac);
    lprintf("BICGSTAB TEST", 0, "BiCGstab converged in = %d iterations\n", cgiters);

    H2(s2, res);
    mul_add_assign(s2, -par.shift[0], res);
    sub_assign(s2, s1);
    tau = sqnorm(s2) / sqnorm(s1);
    lprintf("BICGSTAB TEST", 0, "test BiCGstab = %e (given relative inverter precision: %e)\n", tau, par.err2);
    lprintf("BICGSTAB TEST", 0, "Done\n");

#if defined(DPHI_FLT) && defined(WITH_MPI)
    /*************  SAP  ***************/
    lprintf("SAP TEST", 0, "---------------------\n");
    lprintf("SAP TEST", 0, "Testing SAP only\n");
    lprintf("SAP TEST", 0, "---------------------\n");
    copy(s1, check);
    zero(s2);
    zero(res);

    n_dirac = getMVM() + getMVM_flt();

    timer_lap(&clock);
    int ncy = 5;
    int nit = 4;
    cgiters = SAP_prec(nit, ncy, &BiCGstab, &par, &D, s1, res);
    //cgiters += SAP_prec(nit, ncy, &BiCGstab, &par, &H, tmp, res);
    elapsed = timer_lap(&clock) * 1.e-3;
    n_dirac = getMVM() + getMVM_flt();

    lprintf("SAP TEST", 0, "SAP inversion took %lf msec\n", elapsed);
    lprintf("SAP TEST", 0, "Inverter needed %d hopping term applications.\n", n_dirac);
    lprintf("SAP TEST", 0, "SAP applied with nmr=%d ncy=%d\n", nit, ncy);

    D(s2, res);
    mul_add_assign(s2, -par.shift[0], res);
    sub_assign(s2, s1);
    tau = sqnorm(s2) / sqnorm(s1);
    lprintf("SAP TEST", 0, "SAP needed %lf msec\n", elapsed);
    lprintf("SAP TEST", 0, "test SAP = %e (given relative inverter precision: %e)\n", tau, par.err2);
    lprintf("SAP TEST", 0, "Done\n");

    /*************  GCR  ***************/
    lprintf("GCR TEST", 0, "---------------------\n");
    lprintf("GCR TEST", 0, "Testing GCR\n");
    lprintf("GCR TEST", 0, "---------------------\n");
    copy(s1, check);
    zero(s2);
    zero(res);

    /* Invert and print results */
    n_dirac = getMVM() + getMVM_flt();

    timer_lap(&clock);
    g5(tmp, s1);
    cgiters = gcr(&par, &D, tmp, res);
    g5(tmp, res);
    cgiters += gcr(&par, &D, tmp, res);
    elapsed = timer_lap(&clock) * 1.e-3;
    n_dirac = getMVM() + getMVM_flt();

    lprintf("GCR TEST", 0, "GCR mshift needed %lf msec\n", elapsed);
    lprintf("GCR TEST", 0, "Converged in %d iterations\n", cgiters);
    lprintf("GCR TEST", 0, "Inverter needed %d hopping term applications.\n", n_dirac);
    lprintf("GCR TEST", 0, "Outfield sqnorm: %0.15e\n", sqnorm(res));

    /* check precision reached */
    H2(s2, res);
    mul_add_assign(s2, -par.shift[0], res);
    sub_assign(s2, s1);
    tau = sqnorm(s2) / sqnorm(s1);
    lprintf("SAP TEST", 0, "test SAP = %e (given relative inverter precision: %e)\n", tau, par.err2);

    /*************  SAP+GCR  ***************/
    lprintf("SAP+GCR TEST", 0, "---------------------\n");
    lprintf("SAP+GCR TEST", 0, "Testing SAP+GCR\n");
    lprintf("SAP+GCR TEST", 0, "---------------------\n");
    copy(s1, check);
    zero(s2);
    zero(res);

    n_dirac = getMVM() + getMVM_flt();

    timer_lap(&clock);
    g5(tmp, s1);
    cgiters = sapgcr(&par, &D, tmp, res);
    g5(tmp, res);
    cgiters += sapgcr(&par, &D, tmp, res);
    elapsed = timer_lap(&clock) * 1.e-3;
    n_dirac = getMVM() + getMVM_flt();

    lprintf("SAP+GCR TEST", 0, "SAP+GCR mshift needed %lf msec\n", elapsed);
    lprintf("SAP+GCR TEST", 0, "Converged in %d iterations\n", cgiters);
    lprintf("SAP+GCR TEST", 0, "Inverter needed %d hopping term applications.\n", n_dirac);
    lprintf("SAP+GCR TEST", 0, "Outfield sqnorm: %0.15e\n", sqnorm(res));

    /* check precision reached */
    H2(s2, res);
    sub_assign(s2, s1);
    tau = sqnorm(s2) / sqnorm(s1);
    lprintf("SAP+GCR TEST", 0, "test SAP = %e (given relative inverter precision: %e)\n", tau, par.err2);
#endif

    /*************  Finalize  ***************/

    free_field(res);
    free(par.shift);
    finalize_process();

    return 0;
}
