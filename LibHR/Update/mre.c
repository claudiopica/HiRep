/***************************************************************************
 * Copyright (c) 2014, Martin Hansen                                        *
 * All rights reserved.                                                     *
 *                                                                          *
 * Chronological inverter using the MRE algorithm                           *
 * arXiv: hep-lat/9509012                                                   *
 ***************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "io.h"
#include "Inverters/linear_algebra.h"

//#define cabs(a) sqrt(_complex_prod_re(a,a))
#define MAX 15

// Global variables
static spinor_field *v[MAX];
static spinor_field *Dv;
static int num_init = 0;

// Variables used in in the LU solver
static hr_complex A[MAX][MAX];
static hr_complex b[MAX];
static hr_complex x[MAX];
static hr_complex y[MAX];
static int mutate[MAX];

static void gram_schmidt(mre_par *par, int p, int max) {
    hr_complex rij;
    double rii;

    for (int i = 0; i < max; i++) {
        copy_spinor_field(v[i], &par->s[p][i]);
    }

    for (int i = 0; i < max; i++) {
        rii = sqnorm_spinor_field(v[i]);
        rii = 1.0 / sqrt(rii);
        mul_spinor_field(v[i], rii, v[i]);

        for (int j = i + 1; j < max; j++) {
            rij = prod_spinor_field(v[i], v[j]);
            _complex_minus(rij, rij);
            mulc_add_assign_spinor_field(v[j], rij, v[i]);
        }
    }
}

static void lu_solve(int max) {
    double big;
    int row;
    hr_complex ctmp;
    int itmp;

    // Setup mutate
    for (int i = 0; i < max; i++) {
        mutate[i] = i;
    }

    // LU factorization
    for (int i = 0; i < max; i++) {
        for (int j = 0; j <= i; j++) {
            for (int k = 0; k < j; k++) {
                _complex_mul_sub_assign(A[j][i], A[j][k], A[k][i]);
            }
        }

        big = cabs(A[i][i]);
        row = i;

        for (int j = i + 1; j < max; j++) {
            for (int k = 0; k < i; k++) {
                _complex_mul_sub_assign(A[j][i], A[j][k], A[k][i]);
            }

            if (cabs(A[j][i]) > big) {
                big = cabs(A[j][i]);
                row = j;
            }
        }

        if (big < 1.0e-14) {
            lprintf("MRE", 10, "LU decomposition failed: matrix is singular\n");
            return;
        }

        if (row != i) {
            for (int k = 0; k < max; k++) {
                ctmp = A[row][k];
                A[row][k] = A[i][k];
                A[i][k] = ctmp;
            }

            itmp = mutate[row];
            mutate[row] = mutate[i];
            mutate[i] = itmp;
        }

        for (int k = i + 1; k < max; k++) {
            _complex_div(ctmp, A[k][i], A[i][i]);
            A[k][i] = ctmp;
        }
    }

    // Forward substitution
    for (int i = 0; i < max; i++) {
        y[i] = conj(b[mutate[i]]);
        for (int k = 0; k < i; k++) {
            _complex_mul_sub_assign(y[i], A[i][k], y[k]);
        }
    }

    // Backward substitution
    for (int i = max - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int k = i + 1; k < max; k++) {
            _complex_mul_sub_assign(x[i], A[i][k], x[k]);
        }
        _complex_div(ctmp, x[i], A[i][i]);
        x[i] = ctmp;
    }
}

#if 0
//these are not used anymore
static int factorial(int k)
{
	int f = 1;
	for (; k > 1; k--)
		f *= k;
	return f;
}

static int coefficient(int k, int n)
{
	int c = factorial(n) / (factorial(k) * factorial(n - k));
	if ((k - 1) % 2)
		c = -c;
	return c;
}
#endif

void mre_init(mre_par *par, int max, double prec) {
    if (max == 0) {
        par->max = 0;
        par->init = 0;
        return;
    } else {
        par->max = (max < MAX) ? max : MAX;
        par->init = 1;
    }

    par->s[0] = alloc_spinor_field(par->max, &glat_default);
    par->s[1] = alloc_spinor_field(par->max, &glat_default);
    par->num[0] = 0;
    par->num[1] = 0;

    if (num_init == 0) { Dv = alloc_spinor_field(1, &glat_default); }

    for (int i = num_init; i < par->max; i++) {
        v[i] = alloc_spinor_field(1, &glat_default);
        num_init++;
    }

    if (prec > 1e-14) { lprintf("MRE", 10, "WARNING: Inverter precision should be at least 1e-14 to ensure reversibility!\n"); }

    lprintf("MRE", 10, "Enabled chronological inverter with %d past solutions\n", par->max);
}

void mre_store(mre_par *par, int p, spinor_field *in) {
    spinor_field *tmp;

    if (num_init == 0 || par->init == 0 || par->max <= 0) { return; }

    if (p > 1) {
        lprintf("MRE", 10, "Cannot store solution: p = %d is not valid\n", p);
        return;
    }

    // Shift all vectors by one and store the new vector at position zero
    tmp = &par->s[p][par->max - 1];

    for (int i = (par->max - 1); i > 0; i--) {
        par->s[p][i] = par->s[p][i - 1];
    }

    par->s[p][0] = *tmp;
    copy_spinor_field(&par->s[p][0], in);
    par->num[p]++;
}

void mre_guess(mre_par *par, int p, spinor_field *out, spinor_operator DD, spinor_field *pf) {
    int max;

    if (num_init == 0 || par->init == 0 || par->max <= 0) {
        zero_spinor_field(out);
        return;
    }

    if (p > 1) {
        lprintf("MRE", 10, "Cannot guess solution: p = %d is not valid\n", p);
        return;
    }

    zero_spinor_field(out);
    max = (par->num[p] > par->max) ? par->max : par->num[p];
    gram_schmidt(par, p, max);

    for (int i = 0; i < max; i++) {
        DD(Dv, v[i]);

        for (int j = 0; j < max; j++) {
            A[j][i] = prod_spinor_field(Dv, v[j]);
        }

        b[i] = prod_spinor_field(v[i], pf);
    }

    lu_solve(max);

    for (int i = 0; i < max; i++) {
        //		mul_add_assign_spinor_field(out, coefficient(i+1, max), s[p][i]);
        mulc_add_assign_spinor_field(out, x[i], v[i]);
    }
}
