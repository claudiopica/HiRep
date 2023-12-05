/******************************************************************************
*
* Test of linear algebra against openMP
*
* NOCOMPILE= WITH_GPU
*
******************************************************************************/

#include "libhr.h"

static double spinors_max_difference(spinor_field *r0, spinor_field *r1) {
    double red = 0.;
    _MASTER_FOR(&glattice, ix) {
        suNf_spinor *v1 = _FIELD_AT(r0, ix);
        suNf_spinor *v2 = _FIELD_AT(r1, ix);
        for (int j = 0; j < 4; j++) {
            for (int n = 0; n < NF; n++) {
                if (_complex_prod_re(v1->c[j].c[n] - v2->c[j].c[n], v1->c[j].c[n] - v2->c[j].c[n]) > red) {
                    red = _complex_prod_re(v1->c[j].c[n] - v2->c[j].c[n], v1->c[j].c[n] - v2->c[j].c[n]);
                }
            }
        }
    }
    global_max(&red, 1);
    bcast(&red, 1);

    return red;
}
int main(int argc, char *argv[]) {
    /* setup process id and communications */
    int retval = 0;
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
#ifdef _OPENMP
    int max_nthreads = omp_get_max_threads();
#endif

    double red, res[2], rnd[2];
    hr_complex cres[2];

    spinor_field *s0, *s1, *s2, *s3;

    s0 = alloc_spinor_field(4, &glattice);
    s1 = s0 + 1;
    s2 = s1 + 1;
    s3 = s2 + 1;

#ifdef _OPENMP
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(max_nthreads);
#endif

    ranlxd(rnd, 2);

    gaussian_spinor_field(s1);
    gaussian_spinor_field(s0);

    copy_spinor_field(s2, s0);
    copy_spinor_field(s3, s1);

    add_assign_spinor_field(s3, s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    add_assign_spinor_field(s2, s1);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for add_assign_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    sub_spinor_field(s2, s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    sub_spinor_field(s3, s0, s1);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for sub_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    mul_spinor_field(s2, rnd[0], s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    mul_spinor_field(s3, rnd[0], s0);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for mul_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    copy_spinor_field(s2, s0);
    copy_spinor_field(s3, s0);
    mul_add_assign_spinor_field(s2, rnd[0], s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    mul_add_assign_spinor_field(s3, rnd[0], s0);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for mul_add_assign_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    copy_spinor_field(s2, s0);
    copy_spinor_field(s3, s0);
    mulc_add_assign_spinor_field(s2, rnd[0] + I * rnd[1], s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    mulc_add_assign_spinor_field(s3, rnd[0] + I * rnd[1], s0);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0,
            "Max point difference for mulc_add_assign_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    g5_spinor_field(s3, s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    g5_spinor_field(s2, s0);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for g5_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    lc_spinor_field(s2, rnd[0], s0, rnd[1], s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    lc_spinor_field(s3, rnd[0], s0, rnd[1], s1);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for lc_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = prod_re_spinor_field(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = prod_re_spinor_field(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for prod_re_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = prod_im_spinor_field(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = prod_im_spinor_field(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for prod_im_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    cres[0] = prod_spinor_field(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    cres[1] = prod_spinor_field(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for prod_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            cabs(cres[1] - cres[0]) / cabs(cres[1]));
    if (cabs(cres[1] - cres[0]) / cabs(cres[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = sqnorm_spinor_field(s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = sqnorm_spinor_field(s0);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for sqnorm_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = g5_prod_re_spinor_field(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = g5_prod_re_spinor_field(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for g5_prod_re_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = g5_prod_im_spinor_field(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = g5_prod_im_spinor_field(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for g5_prod_im_spinor_field_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif
    finalize_process();

    return retval;
}
