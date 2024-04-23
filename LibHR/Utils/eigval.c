/***************************************************************************\
* Copyright (c) 2008-2024, Claudio Pica, Martin Hansen                      *
* All rights reserved.                                                      *
\***************************************************************************/

#include "utils.h"
#include "libhr_core.h"
#include "io.h"
#include "inverters.h"
#include "random.h"
#include "memory.h"

/* use inverse power method to find smallest eigenvalue */
int min_eigval(spinor_operator H, geometry_descriptor *type, double *min) {
    spinor_field *s0, *s1, *s2, *stmp;
    mshift_par par;
    int count = 0;
    double eig = 0;
    double old = 0;
    double tmp = 1;
    double norm;

    s0 = alloc_spinor_field(2, type);
    s1 = s0;
    s2 = s1 + 1;
    gaussian_spinor_field(s1);

    *min = 0;
    par.n = 1;
    par.shift = min;
    par.err2 = 1.0e-6;
    par.max_iter = 0;

    while (tmp > 1.0e-4) {
        norm = sqnorm_spinor_field(s1);
        norm = 1.0 / sqrt(norm);
        mul_spinor_field(s1, norm, s1);

        count += cg_mshift(&par, H, s1, s2);

        norm = prod_re_spinor_field(s1, s2);
        old = eig;
        eig = norm;
        tmp = fabs((eig - old) / eig);

        stmp = s2;
        s2 = s1;
        s1 = stmp;
    }

    *min = (1.0 / eig);
    lprintf("EVLIMIT", 10, "min_eigval = %1.8e [MVM = %d]\n", *min, count);

    free_spinor_field(s0);
    return count;
}

/* use power method to find largest eigenvalue */
int max_eigval(spinor_operator H, geometry_descriptor *type, double *max) {
    spinor_field *s0, *s1, *s2, *stmp;
    int count = 0;
    double eig = 0;
    double old = 0;
    double tmp = 1;
    double norm;

    s0 = alloc_spinor_field(2, type);
    s1 = s0;
    s2 = s1 + 1;
    gaussian_spinor_field(s1);

    while (tmp > 1.0e-4) {
        norm = sqnorm_spinor_field(s1);
        norm = 1.0 / sqrt(norm);
        mul_spinor_field(s1, norm, s1);

        H(s2, s1);
        count++;

        norm = prod_re_spinor_field(s1, s2);
        old = eig;
        eig = norm;
        tmp = fabs((eig - old) / eig);

        stmp = s2;
        s2 = s1;
        s1 = stmp;
    }

    *max = (1.0 * eig);
    lprintf("EVLIMIT", 10, "max_eigval = %1.8e [MVM = %d]\n", *max, count);

    free_spinor_field(s0);
    return count;
}
