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

    spinor_field_copy_f(s2, s0);
    spinor_field_copy_f(s3, s1);

    spinor_field_add_assign_f(s3, s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    spinor_field_add_assign_f(s2, s1);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for spinor_field_add_assign_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    spinor_field_sub_f(s2, s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    spinor_field_sub_f(s3, s0, s1);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for spinor_field_sub_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    spinor_field_mul_f(s2, rnd[0], s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    spinor_field_mul_f(s3, rnd[0], s0);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for spinor_field_mul_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    spinor_field_copy_f(s2, s0);
    spinor_field_copy_f(s3, s0);
    spinor_field_mul_add_assign_f(s2, rnd[0], s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    spinor_field_mul_add_assign_f(s3, rnd[0], s0);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for spinor_field_mul_add_assign_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    spinor_field_copy_f(s2, s0);
    spinor_field_copy_f(s3, s0);
    spinor_field_mulc_add_assign_f(s2, rnd[0] + I * rnd[1], s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    spinor_field_mulc_add_assign_f(s3, rnd[0] + I * rnd[1], s0);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0,
            "Max point difference for spinor_field_mulc_add_assign_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    spinor_field_g5_f(s3, s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    spinor_field_g5_f(s2, s0);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for spinor_field_g5_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    spinor_field_lc_f(s2, rnd[0], s0, rnd[1], s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    spinor_field_lc_f(s3, rnd[0], s0, rnd[1], s1);
    red = spinors_max_difference(s3, s2);
    lprintf("MAIN", 0, "Max point difference for spinor_field_lc_f %.14e.\n(should be around 1*10^(-15) or so)\n\n", red);
    if (red > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = spinor_field_prod_re_f(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = spinor_field_prod_re_f(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for spinor_field_prod_re_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = spinor_field_prod_im_f(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = spinor_field_prod_im_f(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for spinor_field_prod_im_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    cres[0] = spinor_field_prod_f(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    cres[1] = spinor_field_prod_f(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for spinor_field_prod_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            cabs(cres[1] - cres[0]) / cabs(cres[1]));
    if (cabs(cres[1] - cres[0]) / cabs(cres[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = spinor_field_sqnorm_f(s0);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = spinor_field_sqnorm_f(s0);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for spinor_field_sqnorm_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = spinor_field_g5_prod_re_f(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = spinor_field_g5_prod_re_f(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for spinor_field_g5_prod_re_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif

    res[0] = spinor_field_g5_prod_im_f(s0, s1);
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    res[1] = spinor_field_g5_prod_im_f(s0, s1);
    bcast(res, 2);
    lprintf("MAIN", 0, "Max difference for spinor_field_g5_prod_im_f %.14e.\n(should be around 1*10^(-15) or so)\n\n",
            fabs(res[1] - res[0]) / fabs(res[1]));
    if (fabs(res[1] - res[0]) / fabs(res[1]) > 1.e-14) { retval++; }
#ifdef _OPENMP
    omp_set_num_threads(max_nthreads);
#endif
    finalize_process();

    return retval;
}
