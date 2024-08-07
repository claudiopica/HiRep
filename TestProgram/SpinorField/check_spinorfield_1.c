
/******************************************************************************
 *
 * File check11.c
 *
 * Consistency checks on the programs in the module linalg
 *
 * Author: luigi del debbio <luigi.del.debbio@ed.ac.uk>
 *
 ******************************************************************************/

#include "libhr.h"

#define MAX_ROTATE 50

static hr_complex v[25];
static double EPSILON = 1.e-12;
static spinor_field *ppk[5];

static int initr = 0;

static suNf_spinor *psi;

static void alloc_ws_rotate(void) {
    psi = calloc(MAX_ROTATE, sizeof(suNf_spinor));

    error((psi == NULL), 1, "alloc_ws_rotate [linalg.c]", "Unable to allocate workspace");

    initr = 1;
}

static void rotate_ptr(int n, spinor_field *pkk[], hr_complex vl[]) {
    if (initr == 0) { alloc_ws_rotate(); }
#ifdef WITH_GPU
    for (int i = 0; i < n; i++) {
        copy_from_gpu(pkk[i]);
    }
#endif

    error((n < 1) || (n > MAX_ROTATE), 1, "rotate [eva.c]", "Parameter n is out of range");

    for (int i = 0; i < n; i++) {
        error((*pkk)->type != (*(pkk + i))->type, 1, "not available", "Spinors don't match!");
    }

#undef _OMP_PRAGMA //This doesn't works with multiple threads
#define _OMP_PRAGMA(s)

    _MASTER_FOR(pkk[0]->type, ix) {
        for (int k = 0; k < n; k++) {
            suNf_spinor *pk = &(psi[k]);
            suNf_spinor *pj = _FIELD_AT(pkk[0], ix);
            hr_complex *z = &vl[k];

            _spinor_mulc_f(*pk, *z, *pj);

            for (int j = 1; j < n; j++) {
                pj = _FIELD_AT(pkk[j], ix);
                z += n;

                _spinor_mulc_add_assign_f(*pk, *z, *pj);
            }
        }

        for (int k = 0; k < n; k++) {
            *_FIELD_AT(pkk[k], ix) = psi[k];
        }
    }

#ifdef WITH_GPU
    for (int i = 0; i < n; i++) {
        copy_to_gpu(pkk[i]);
    }
#endif
}

static void project(spinor_field *pk, spinor_field *pl) {
    hr_complex sp;

    sp = -prod_spinor_field(pl, pk);

    mulc_add_assign_spinor_field(pk, sp, pl);
}

static double normalize(spinor_field *ps) {
    double r, ri;

    r = sqnorm_spinor_field(ps);
    r = sqrt(r);
    error(r < EPSILON, 1, "normalize [eva.c]", "vector has vanishing norm");

    ri = 1.0 / r;
    mul_spinor_field(ps, ri, ps);

    return (double)(r);
}

static hr_complex sp(spinor_field *pk, spinor_field *pl) {
    hr_complex z = 0.0;

    _TWO_SPINORS_FOR_SUM(pk, pl, x, y) {
        for (int i = 0; i < (4 * NF); i++) {
            hr_complex *rpk = (hr_complex *)_SPINOR_PTR(pk) + i;
            hr_complex *rpl = (hr_complex *)_SPINOR_PTR(pl) + i;
            z += conj(*rpk) * (*rpl);
            /* x+=(double)((*rpk).re*(*rpl).re+(*rpk).im*(*rpl).im); */
            /* y+=(double)((*rpk).re*(*rpl).im-(*rpk).im*(*rpl).re); */
            //rpk+=1; //?? why these increment
            //rpl+=1; //??
        }
    }

#ifdef WITH_MPI
    global_sum((double *)&z, 2);
#endif

    return z;
}

int main(int argc, char *argv[]) {
    int i, j;
    double r;
    double rd, zsqd;
    double d, dmax;
    hr_complex w;
    hr_complex zd, wd;
    spinor_field *ws;
    spinor_field *pk, *pl;
    spinor_field *tmp;
    int return_value = 0;

    logger_map("DEBUG", "debug");

    /* setup process id and communications */
    setup_process(&argc, &argv);

    lprintf("LA TEST", 0, "Consistency of the programs in the module linalg\n");
    lprintf("LA TEST", 0, "------------------------------------------------\n");

    tmp = alloc_spinor_field(1, &glattice);
    ws = alloc_spinor_field(10, &glattice);

    for (i = 0; i < 10; i++) {
        gaussian_spinor_field(&ws[i]);
    }

    dmax = 0.0;

    for (i = 0; i < 10; i++) {
        pk = &ws[i];
        pl = &ws[9 - i];
        w = sp(pk, pl);

        zd = prod_spinor_field(pk, pl);
        rd = sqnorm_spinor_field(pk) * sqnorm_spinor_field(pl);
        d = (zd - w) * conj(zd - w);
        /* d=((zd.re-(double)w.re)*(zd.re-(double)w.re)+ */
        /*    (zd.im-(double)w.im)*(zd.im-(double)w.im)); */
        d = sqrt(d / rd);
        if (d > dmax) { dmax = d; }

        rd = prod_re_spinor_field(pk, pl);
        d = fabs(creal(zd) / rd - 1.0);
        if (d > dmax) { dmax = d; }

        zd = prod_spinor_field(pk, pk);
        rd = sqnorm_spinor_field(pk);

        d = fabs(cimag(zd) / rd);
        if (d > dmax) { dmax = d; }

        d = fabs(creal(zd) / rd - 1.0f);
        if (d > dmax) { dmax = d; }
    }
    lprintf("LA TEST", 0, "Check of spinor_field_prod, spinor_field_prod_re\n");
    lprintf("LA TEST", 0, "and spinor_field_sqnorm: %.2e\n\n", dmax);
    if (dmax > 1e-14) {
        lprintf("LA TEST", 0, "Test failed ?\n");
        return_value += 1;
    }

    dmax = 0.0;
    zd = 0.345 - I * 0.876;
    zsqd = zd * conj(zd);

    for (i = 0; i < 9; i++) {
        pk = &ws[i];
        pl = &ws[i + 1];

        wd = prod_spinor_field(pk, pl);
        rd = sqnorm_spinor_field(pk) + zsqd * sqnorm_spinor_field(pl) + 2.0 * (creal(zd * wd));

        mulc_add_assign_spinor_field(pk, zd, pl);

        d = fabs(rd / sqnorm_spinor_field(pk) - 1.0);
        if (d > dmax) { dmax = d; }
    }
    lprintf("LA TEST", 0, "Consistency of spinor_prod, norm_square\n");
    lprintf("LA TEST", 0, "and mulc_spinor_add: %.2e\n\n", dmax);
    if (dmax > 1e-14) {
        lprintf("LA TEST", 0, "Test failed ?\n");
        return_value += 1;
    }
    for (i = 0; i < 10; i++) {
        gaussian_spinor_field(&ws[i]);
    }

    dmax = 0.0;

    for (i = 0; i < 10; i++) {
        pk = &ws[i];

        if (i > 0) {
            pl = &ws[i - 1];
            project(pk, pl);
            zd = prod_spinor_field(pk, pl);

            d = (fabs(creal(zd)) + fabs(cimag(zd))) / sqrt(sqnorm_spinor_field(pk));

            if (d > dmax) { dmax = d; }
        }

        normalize(pk);
        rd = sqnorm_spinor_field(pk);

        d = fabs(rd - 1.0f);
        if (d > dmax) { dmax = d; }
    }

    lprintf("LA TEST", 0, "Consistency of spinor_prod, norm_square,\n");
    lprintf("LA TEST", 0, "normalize and project: %.2e\n\n", dmax);
    if (dmax > 1e-14) {
        lprintf("LA TEST", 0, "Test failed ?\n");
        return_value += 1;
    }

    for (i = 0; i < 5; i++) {
        pk = &ws[i];
        pl = &ws[i + 5];

        gaussian_spinor_field(pk);
        copy_spinor_field(pl, pk);

        for (j = 0; j < 5; j++) {
            v[5 * i + j] =
                0.1234f * (double)(i ^ 2) - 0.8976f * (double)(j) + I * (0.2231f * (double)(i) + 0.9922f * (double)(j ^ 2));
        }

        ppk[i] = pl;
    }

    rotate_ptr(5, ppk, v);
    dmax = 0.0;

    for (i = 5; i < 10; i++) {
        pk = &ws[i];

        for (j = 0; j < 5; j++) {
            zd = -v[5 * j + (i - 5)];

            pl = &ws[j];
            mulc_add_assign_spinor_field(pk, zd, pl);
        }

        rd = sqnorm_spinor_field(pk);

        d = fabs(rd);
        if (d > dmax) { dmax = d; }
    }

    dmax /= sqnorm_spinor_field(&ws[0]);
    dmax = sqrt(dmax);

    lprintf("LA TEST", 0, "Consistency of mulc_spinor_add\n");
    lprintf("LA TEST", 0, "and rotate: %.2e\n\n", dmax);
    if (dmax > 1e-14) {
        lprintf("LA TEST", 0, "Test failed ?\n");
        return_value += 1;
    }

    dmax = 0.0;

    for (i = 0; i < 5; i++) {
        pk = &ws[i];
        pl = &ws[9 - i];
        gaussian_spinor_field(pk);
        copy_spinor_field(pl, pk);
        g5_spinor_field(tmp, pk);
        g5_spinor_field(pk, tmp);

        zd = -1.0;

        mulc_add_assign_spinor_field(pl, zd, pk);
        r = sqnorm_spinor_field(pl) / sqnorm_spinor_field(pk);
        d = sqrt(r);
        if (d > dmax) { dmax = d; }

        gaussian_spinor_field(pl);
        zd = prod_spinor_field(pk, pl);
        g5_spinor_field(pk, pk);
        g5_spinor_field(pl, pl);
        wd = prod_spinor_field(pk, pl);

        d = (fabs(creal(zd - wd)) + fabs(cimag(zd - wd))) / (fabs(creal(zd)) + fabs(cimag(zd)));
        if (d > dmax) { dmax = d; }
    }

    lprintf("LA TEST", 0, "Check of g5_spinor_field_f: %.2e\n\n", dmax);
    if (dmax > 1e-30) {
        lprintf("LA TEST", 0, "Test failed ?\n");
        return_value += 1;
    }

    dmax = 0.0;

    for (i = 0; i < 5; i++) {
        pk = &ws[i];
        pl = &ws[9 - i];
        gaussian_spinor_field(pk);
        copy_spinor_field(pl, pk);
        d = -2.5;
        lc1_spinor_field(d, pk, pl);

        zd = 1.5;
        mulc_add_assign_spinor_field(pk, zd, pl);
        d = sqnorm_spinor_field(pk) / sqnorm_spinor_field(pl);

        if (d > dmax) { dmax = d; }
    }

    lprintf("LA TEST", 0, "Check of lc1: %.2e\n\n", dmax);
    if (dmax > 1e-31) {
        lprintf("LA TEST", 0, "Test failed ?\n");
        return_value += 1;
    }
    dmax = 0.0;

    for (i = 0; i < 5; i++) {
        pk = &ws[i];
        pl = &ws[9 - i];
        gaussian_spinor_field(pk);
        copy_spinor_field(pl, pk);
        d = 1.0;
        r = 2.5;
        lc2_spinor_field(d, r, pk, pl);

        zd = -3.5;
        mulc_add_assign_spinor_field(pk, zd, pl);
        d = sqnorm_spinor_field(pk) / sqnorm_spinor_field(pl);

        if (d > dmax) { dmax = d; }
    }

    lprintf("LA TEST", 0, "Check of lc2: %.2e\n\n", dmax);
    if (dmax > 1e-31) {
        lprintf("LA TEST", 0, "Test failed ?\n");
        return_value += 1;
    }
    dmax = 0.0;

    for (i = 0; i < 5; i++) {
        pk = &ws[i];
        pl = &ws[9 - i];
        gaussian_spinor_field(pk);
        copy_spinor_field(pl, pk);
        d = 3.5;
        r = -1.5;
        lc3_spinor_field(d, r, pk, pl, pk);

        zd = -1.0;

        mulc_add_assign_spinor_field(pk, zd, pl);
        d = sqnorm_spinor_field(pk) / sqnorm_spinor_field(pl);

        if (d > dmax) { dmax = d; }
    }

    lprintf("LA TEST", 0, "Check of lc3: %.2e\n\n", dmax);
    if (dmax > 1e-31) {
        lprintf("LA TEST", 0, "Test failed ?\n");
        return_value += 1;
    }
    finalize_process();

    return return_value;
}
