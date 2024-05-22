/******************************************************************************
*
* Test of Schwarz Alternating Procedure
*
* NOCOMPILE= !WITH_GPU
* NOCOMPILE= !WITH_MPI
* NOCOMPILE= !DPHI_FLT
*
******************************************************************************/

#include "libhr.h"

int nhb, nor, nit, nth, nms, level, seed;
double beta;
spinor_field *tmp;
void M(spinor_field *out, spinor_field *in) {
    H(out, in);
}

int main(int argc, char *argv[]) {
    setup_process(&argc, &argv);
    setup_gauge_fields();

    double tau;
    spinor_field *s1, *s2;
    spinor_field *res;

    mshift_par par;
    int cgiters;
    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);
    represent_gauge_field();

    par.n = 1;
    par.shift = (double *)malloc(sizeof(double) * (par.n));
    par.err2 = 1.e-28;
    par.max_iter = 0;
#ifndef UPDATE_EO
    tmp = alloc_spinor_field(4, &glattice);
#else
    tmp = alloc_spinor_field(4, &glat_even);
#endif
    res = tmp + 1;
    s1 = res + 1;
    s2 = s1 + 1;
    par.shift[0] = 0.0;
    gaussian_spinor_field(s1);

    lprintf("SAP TEST", 0, "Testing SAP\n");
    lprintf("SAP TEST", 0, "---------------------\n");
    cgiters = 0;
    zero_spinor_field(res);

    cgiters = SAP_prec(10, 5, &BiCGstab, &par, &D, s1, res);

    lprintf("SAP TEST", 0, "Converged in %d iterations\n", cgiters);
    D(s2, res);
    printf("sqnorm res %0.15e\n", sqnorm(s1));
    sub_assign_spinor_field(s2, s1);
    tau = sqnorm_spinor_field(s2) / sqnorm_spinor_field(s1);
    lprintf("SAP TEST", 0, "test SAP = %e (req. %e)\n", tau, par.err2);

    lprintf("SAP TEST", 0, "Testing GCR\n");
    lprintf("SAP TEST", 0, "---------------------\n");
    cgiters = 0;
    zero_spinor_field(res);
    g5(tmp, s1);
    cgiters = gcr(&par, &D, tmp, res);
    g5(tmp, res);
    cgiters += gcr(&par, &D, tmp, res);

    lprintf("GCR TEST", 0, "Converged in %d iterations\n", cgiters);
    H2(s2, res);
    sub_assign_spinor_field(s2, s1);
    tau = sqnorm_spinor_field(s2) / sqnorm_spinor_field(s1);
    lprintf("GCR TEST", 0, "test GCR = %e (req. %e)\n", tau, par.err2);

    lprintf("SAP TEST", 0, "Testing SAP+GCR\n");
    lprintf("SAP TEST", 0, "---------------------\n");
    cgiters = 0;
    zero_spinor_field(res);

    g5(tmp, s1);
    cgiters = sapgcr(&par, &D, tmp, res);
    g5(tmp, res);
    cgiters += sapgcr(&par, &D, tmp, res);

    lprintf("SAP+GCR TEST", 0, "Converged in %d iterations\n", cgiters);
    H2(s2, res);
    sub_assign_spinor_field(s2, s1);
    tau = sqnorm_spinor_field(s2) / sqnorm_spinor_field(s1);
    lprintf("SAP+GCR TEST", 0, "test SAP+GCR = %e (req. %e)\n", tau, par.err2);

    free_spinor_field(res);
    free(par.shift);

    finalize_process();
    return 0;
}
