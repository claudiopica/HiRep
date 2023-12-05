/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "libhr_core.h"
#include "random.h"
#include "io.h"
#include "memory.h"
#include "Utils/single_double_utils.h"
#include <math.h>
#include <assert.h>

// Uncomment to use BiCGstab instead of MINRES (does not work with multishift!)
//#define USE_BICGSTAB

//TODO: why is this here?
#undef NDEBUG

#define _print_par(c) printf("[%d] " #c " = %e\n", cgiter, c)

/* Calculate the Givens rotation (c,s) that maps:
* |  c   s | | a | = | 1/r |
* | -s^+ c | | b |   | 0 |
*/
#define _Givens_rot(r, c, s, a, b)              \
    {                                           \
        if (fabs(a) < 1.e-15) {                 \
            (r) = 1. / (b);                     \
            (c) = 0.;                           \
            (s) = 1.;                           \
        } else {                                \
            (s) = (b) / (a);                    \
            (c) = 1. / sqrt(1. + (s) * (s));    \
            (s) *= (c);                         \
            (r) = 1. / ((c) * (a) + (s) * (b)); \
        }                                       \
    }

/*
* performs the multi-shifted QMR inversion for g5-hermitean matrices:
* out[i] = (M-(par->shift[i]))^-1 in
* returns the number of cg iterations done.
*/
static int g5QMR_mshift_core(short *valid, mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out) {
    spinor_field **q1, **q2;
    spinor_field *p1, *p2, *Mp;
    spinor_field *sptmp, *memall;

    double alpha, beta, delta, rho, innorm2;
    double *r, *s1, *s2, *c1, *c2;
    double maxm;

    int cgiter;
    unsigned int notconverged;

    unsigned short *flags;

    /* fare qualche check sugli input */
    /* par->n deve essere almeno 1! */
    assert(par->n > 0);
    _TWO_SPINORS_MATCHING(in, &out[0]);
    _ARRAY_SPINOR_MATCHING(out, par->n);

    /* allocate spinors fields and aux real variables */
    /* implementation note: to minimize the number of malloc calls
  * objects of the same type are allocated together
  */
    memall = alloc_spinor_field(2 * (par->n) + 3, in->type);
    q1 = (spinor_field **)malloc(sizeof(spinor_field *) * par->n);
    q2 = (spinor_field **)malloc(sizeof(spinor_field *) * par->n);
    for (int i = 0; i < par->n; i++) {
        q1[i] = memall + i;
        q2[i] = memall + par->n + i;
    }
    p1 = memall + 2 * par->n;
    p2 = p1 + 1;
    Mp = p2 + 1;

    r = (double *)malloc(sizeof(double) * 5 * (par->n));
    s1 = r + (par->n);
    s2 = s1 + (par->n);
    c1 = s2 + (par->n);
    c2 = c1 + (par->n);

    flags = (unsigned short *)malloc(sizeof(unsigned short) * (par->n));

    /* init recursion */
    cgiter = 0;
    notconverged = par->n;

    copy_spinor_field(p2, in); /* trial solution = 0 */
    innorm2 = sqnorm_spinor_field(in);
    if (par->n == 1) { /* in this case is not a multishift and we use as starting vector out[0] */
        M(Mp, &out[0]);
        mul_add_assign_spinor_field(Mp, -par->shift[0], &out[0]);
        sub_spinor_field(p2, p2, Mp);
    }
    rho = sqrt(sqnorm_spinor_field(p2));
    lprintf("INVERTER", 50, "g5QMR: rho init: %1.8e\n", rho * rho / innorm2);

    mul_spinor_field(p2, 1. / rho, p2);
    zero_spinor_field(p1);
    for (int i = 0; i < (par->n); ++i) {
        r[i] = rho;
        c2[i] = c1[i] = 1.;
        s1[i] = s2[i] = 0.;
        if (par->n != 1) { zero_spinor_field(&out[i]); /* if no multishift we start with the trial solution */ }
        zero_spinor_field(q1[i]);
        zero_spinor_field(q2[i]);
        flags[i] = 1;
    }
    delta = g5_prod_re_spinor_field(p2, p2);
    error(fabs(delta) < 1.e-28, 1, "g5QMR " __FILE__, "g5 norm of the initial guess is zero, change the initial guess");
    beta = delta;

    /* cg recursion */
    do {
        ++cgiter;

        M(Mp, p2);
        mul_add_assign_spinor_field(Mp, -par->shift[0], p2);

        /* compute alpha */
        alpha = g5_prod_re_spinor_field(p2, Mp) / delta;

        /* update p1, p2 */
        mul_add_assign_spinor_field(Mp, -beta, p1);
        mul_add_assign_spinor_field(Mp, -alpha, p2);
        sptmp = p1;
        p1 = p2;
        p2 = Mp;
        Mp = sptmp;

        /* update rho */
        rho = sqrt(sqnorm_spinor_field(p2));

        maxm = 1.e-10; /* to check if the error is going down */

        for (int i = 0; i < (par->n); ++i) { /* update solutions */
            if (flags[i]) {
                double a, t, e, d, m;
                a = (alpha - par->shift[i] + par->shift[0]);
                t = s1[i] * beta;
                e = c2[i] * c1[i] * beta + s2[i] * a;
                m = c2[i] * a - s2[i] * c1[i] * beta;
                s1[i] = s2[i];
                c1[i] = c2[i];
                _Givens_rot(d, c2[i], s2[i], m, rho);

                maxm = (maxm > fabs(m)) ? maxm : fabs(m);

                /* update q */
                lc_spinor_field(Mp, -t * d, q1[i], -e * d, q2[i]);
                mul_add_assign_spinor_field(Mp, d, p1);
                sptmp = q1[i];
                q1[i] = q2[i];
                q2[i] = Mp;
                Mp = sptmp; /* swap q1[i]<-q2[i]<-Mp and Mp point to q1 */

                /* update solution */
                mul_add_assign_spinor_field(&out[i], c2[i] * r[i], q2[i]);

                /* update residuum */
                r[i] *= -s2[i];

                if (fabs(r[i] * r[i]) * (double)(1 + cgiter) < par->err2 * innorm2) {
                    flags[i] = 0;
                    --notconverged;
                    lprintf("INVERTER", 50, "g5QMR: vect %d converged at iter: %d\n", i, cgiter);
                } else {
                    lprintf("INVERTER", 50, "g5QMR: vect %d iter: %d res: %1.8e s2: %1.8e m=%1.8e\n", i, cgiter,
                            fabs(r[i] * r[i]) * (double)(1 + cgiter) / innorm2, s2[i], m);
                }
            }
        }

        if (fabs(rho) > 1.e-13) {
            double olddelta;

            /*normalize p2 */
            mul_spinor_field(p2, 1. / rho, p2);

            /* update delta and beta */
            olddelta = delta;
            delta = g5_prod_re_spinor_field(p2, p2);

            lprintf("INVERTER", 40, "g5QMR: delta=%1.8e rho=%1.8e [iter=%d]\n", delta, rho, cgiter);
            if (fabs(delta) < 1.e-6 || maxm < 1.e-3) { /* the method has failed ! */
                lprintf("INVERTER", 40, "g5QMR: method failed! delta=%e [iter=%d]\n", delta, cgiter);
                notconverged = 0; /* exit loop */
            }

            /* update beta */
            beta = rho * delta / olddelta;

        } else { /* system has converged */
            notconverged = 0;
            lprintf("INVERTER", 40, "g5QMR: rho < 1.e-13 ! (system converged)\n");
        }

#ifndef NDEBUG
        if (cgiter % 100 == 0) {
            lprintf("INVERTER", 40, "g5QMR: [%d] res[0]=%e notconverged=%d\n", cgiter,
                    fabs(r[0] * r[0]) * (double)(1 + cgiter) / innorm2, notconverged);
        }
#endif

    } while ((par->max_iter == 0 || cgiter < par->max_iter) && notconverged);

    /* test results */
    for (int i = 0; i < par->n; ++i) {
        double norm;
        M(Mp, &out[i]);
        ++cgiter;
        if (par->shift[i] != 0.) { mul_add_assign_spinor_field(Mp, -par->shift[i], &out[i]); }
        sub_spinor_field(Mp, Mp, in);
        norm = sqnorm_spinor_field(Mp) / innorm2;
        valid[i] = 1;
        if (fabs(norm) > par->err2) {
            valid[i] = 0;
            lprintf("INVERTER", 30, "g5QMR failed on vect %d: err2 = %1.8e > %1.8e\n", i, norm, par->err2);
        } else {
            lprintf("INVERTER", 20, "g5QMR inversion: err2 = %1.8e < %1.8e\n", norm, par->err2);
        }
    }

    /* free memory */
    free_spinor_field(memall);
    free(q1);
    free(q2);
    free(r);

    free(flags);

    /* return number of cg iter */
    if (par->n == 1) { ++cgiter; }
    return cgiter;
}

static int g5QMR_core_flt(short *valid, double err2, int max_iter, spinor_operator_flt M, spinor_field_flt *in,
                          spinor_field_flt *out) {
    spinor_field_flt *q1, *q2;
    spinor_field_flt *p1, *p2, *Mp;
    spinor_field_flt *sptmp, *memall;

    double alpha, beta, delta, rho, innorm2;
    double r, s1, s2, c1, c2;
    double maxm;

    int cgiter;
    unsigned int notconverged;

    unsigned short flag;

    /* fare qualche check sugli input */
    _TWO_SPINORS_MATCHING(in, out);

    /* allocate spinors fields and aux real variables */
    /* implementation note: to minimize the number of malloc calls
  * objects of the same type are allocated together
  */
    memall = alloc_spinor_field_flt(2 + 3, in->type);
    q1 = memall;
    q2 = q1 + 1;
    p1 = q2 + 1;
    p2 = p1 + 1;
    Mp = p2 + 1;

#ifndef NDEBUG
    spinor_field_flt *sdbg = alloc_spinor_field_flt(1, in->type);
#endif

    /* init recursion */
    cgiter = 0;
    notconverged = 1;

    copy_spinor_field_flt(p2, in); /* trial solution = 0 */
    innorm2 = sqnorm_spinor_field_flt(in);
    M(Mp, out);
    sub_spinor_field_flt(p2, p2, Mp);
    rho = sqrt(sqnorm_spinor_field_flt(p2));
    lprintf("INVERTER", 60, "g5QMR_core_flt: rho init: %1.8e\n", rho * rho / innorm2);

    mul_spinor_field_flt(p2, 1.f / ((float)(rho)), p2);
    zero_spinor_field_flt(p1);
    r = rho;
    c2 = c1 = 1.;
    s1 = s2 = 0.;
    zero_spinor_field_flt(q1);
    zero_spinor_field_flt(q2);
    flag = 1;
    delta = g5_prod_re_spinor_field_flt(p2, p2);
    beta = delta;

    /* cg recursion */
    do {
        ++cgiter;

        M(Mp, p2);

        /* compute alpha */
        alpha = g5_prod_re_spinor_field_flt(p2, Mp) / delta;

        /* update p1, p2 */
        mul_add_assign_spinor_field_flt(Mp, ((float)(-beta)), p1);
        mul_add_assign_spinor_field_flt(Mp, -((float)(alpha)), p2);
        sptmp = p1;
        p1 = p2;
        p2 = Mp;
        Mp = sptmp;

        /* update rho */
        rho = sqrt(sqnorm_spinor_field_flt(p2));

        maxm = 1.e-10; /* to check if the error is going down */

        if (flag) {
            double a, t, e, d, m;
            a = alpha;
            t = s1 * beta;
            e = c2 * c1 * beta + s2 * a;
            m = c2 * a - s2 * c1 * beta;
            s1 = s2;
            c1 = c2;
            if (fabs(m) < 1.e-15) {
                d = 1. / rho;
                c2 = 0.;
                s2 = 1.;
            } else {
                s2 = rho / m;
                c2 = 1. / sqrt(1. + s2 * s2);
                s2 *= c2;
                d = 1. / (c2 * m + s2 * rho);
            }

            maxm = (maxm > fabs(m)) ? maxm : fabs(m);

            /* update q */
            lc_spinor_field_flt(Mp, ((float)(-t * d)), q1, ((float)(-e * d)), q2);
            mul_add_assign_spinor_field_flt(Mp, ((float)(d)), p1);
            sptmp = q1;
            q1 = q2;
            q2 = Mp;
            Mp = sptmp; /* swap q1<-q2<-Mp and Mp point to q1 */

            /* update solution */
            mul_add_assign_spinor_field_flt(out, ((float)(c2 * r)), q2);

            /* update residuum */
            r *= -s2;

            if (fabs(r * r) * (double)(1 + cgiter) < err2 * innorm2) {
                flag = 0;
                --notconverged;
                lprintf("INVERTER", 60, "g5QMR_core_flt: vect converged at iter: %d\n", cgiter);
            } else {
                lprintf("INVERTER", 60, "g5QMR_core_flt: iter: %d res: %1.8e s2: %1.8e m=%1.8e\n", cgiter,
                        fabs(r * r) * (double)(1 + cgiter) / innorm2, s2, m);
            }
        }

        if (fabs(rho) > 1.e-5) {
            double olddelta;

            /*normalize p2 */
            mul_spinor_field_flt(p2, ((float)(1. / rho)), p2);

            /* update delta and beta */
            olddelta = delta;
            delta = g5_prod_re_spinor_field_flt(p2, p2);

            lprintf("INVERTER", 60, "g5QMR_core_flt: delta=%1.8e rho=%1.8e [iter=%d]\n", delta, rho, cgiter);
            if (fabs(delta) < 1.e-6 || maxm < 1.e-3) { /* the method has failed ! */
                lprintf("INVERTER", 60, "g5QMR_core_flt: method failed! delta=%e [iter=%d]\n", delta, cgiter);
                notconverged = 0; /* exit loop */
            }

            /* update beta */
            beta = rho * delta / olddelta;

        } else { /* system has converged */
            notconverged = 0;
            lprintf("INVERTER", 60, "g5QMR_core_flt: rho < 1.e-5 ! (system converged)\n");
        }

#ifndef NDEBUG
        if (cgiter % 1 == 0) {
            M(sdbg, out);
            sub_assign_spinor_field_flt(sdbg, in);
            lprintf("INVERTER", 20, "g5QMR_core_flt: [%d] res=%e notconverged=%d\n", cgiter,
                    sqnorm_spinor_field_flt(sdbg) / innorm2, notconverged);
        }
#endif

    } while ((max_iter == 0 || cgiter < max_iter) && notconverged);

    /* test results */
    double norm;
    M(Mp, out);
    ++cgiter;
    sub_spinor_field_flt(Mp, Mp, in);
    norm = sqnorm_spinor_field_flt(Mp) / innorm2;
    *valid = 1;
    if (fabs(norm) > err2) {
        *valid = 0;
        lprintf("INVERTER", 50, "g5QMR_core_flt failed: err2 = %1.8e > %1.8e\n", norm, err2);
    } else {
        lprintf("INVERTER", 50, "g5QMR_core_flt inversion: err2 = %1.8e < %1.8e\n", norm, err2);
    }

    /* free memory */
    free_spinor_field_flt(memall);

#ifndef NDEBUG
    free_spinor_field_flt(sdbg);
#endif

    /* return number of cg iter */
    ++cgiter;
    return cgiter;
}

static double sh;
static spinor_operator g5Herm;
static void Herm(spinor_field *out, spinor_field *in) {
    g5Herm(out, in);
    if (sh != 0.) { mul_add_assign_spinor_field(out, -sh, in); }
    g5_spinor_field(out, out);
}

static double sh_flt;
static spinor_operator_flt g5Herm_flt;
static void Herm_flt(spinor_field_flt *out, spinor_field_flt *in) {
    g5Herm_flt(out, in);
    if (sh_flt != 0.) { mul_add_assign_spinor_field_flt(out, ((float)(-sh_flt)), in); }
    g5_spinor_field_flt(out, out);
}

int g5QMR_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out) {
    int cgiter;
    int n;
    mshift_par orig;
    short *valid;
    int loccg;
    int msiter;

    orig = *par; /* save par */
    valid = malloc(sizeof(short) * orig.n);
    cgiter = g5QMR_mshift_core(valid, par, M, in, out);
    msiter = cgiter;

#ifdef USE_BICGSTAB
    if (par->n == 1) {
        if (valid[0] == 0) { cgiter += BiCGstab(par, M, in, out); }
    } else {
#endif
        /* if some vector has not converged try non-multishift */
        par->n = 1;
        for (n = 0; n < orig.n; ++n) {
            if (valid[n] == 0) {
#ifndef NDEBUG
                lprintf("INVERTER", 20, "Doing non multishift on vector %d\n", n);
#endif
                par->shift = orig.shift + n;
                par->max_iter = 0;
                loccg = g5QMR_mshift_core(valid + n, par, M, in, out + n);
                cgiter += loccg;

                if (valid[n] == 0) {
                    /* revert to MINRES */
                    MINRES_par Mpar;

                    lprintf("INVERTER", 20, "g5QMR_mshift: using MINRES on vect %d\n", n);

                    Mpar.err2 = par->err2;
                    Mpar.max_iter = 0;

                    g5Herm = M;
                    sh = par->shift[0];

                    g5_spinor_field(in, in); /* multiply input by g5 for MINRES */
                    loccg = MINRES(&Mpar, &Herm, in, &out[n], (n == 0) ? &out[n] : &out[n - 1]);
                    g5_spinor_field(in, in); /* restore input vector */

                    cgiter += loccg;
                }
            }
        }
#ifdef USE_BICGSTAB
    }
#endif

    /* this is for debug purposes */
    *par = orig; /* restore par */
    free(valid);
    lprintf("INVERTER", 10, "g5QMR_mshift: cgiter (mshift,tot) = %d ; %d\n", msiter, cgiter);

    return cgiter;
}

int g5QMR_fltacc(g5QMR_fltacc_par *par, spinor_operator M, spinor_operator_flt M_flt, spinor_field *in, spinor_field *out) {
    int cgiter = 0, cgiter_flt = 0, cgiter_minres = 0, k;
    short valid;
    spinor_field_flt *in_flt, *out_flt, *res_flt;
    spinor_field *res;
    double err2, innorm2;

    res = alloc_spinor_field(1, in->type);
    in_flt = alloc_spinor_field_flt(3, in->type);
    out_flt = in_flt + 1;
    res_flt = out_flt + 1;
    innorm2 = sqnorm_spinor_field(in);

    /* 0. out = 0 ; res = in */
    copy_spinor_field(res, in);
    zero_spinor_field(out);

    if (par->err2_flt <= 0.) {
        mshift_par mpar;
        double shift = 0.;
        mpar.n = 1;
        mpar.shift = &shift;
        mpar.err2 = par->err2;
        cgiter = g5QMR_mshift(&mpar, M, in, out);
        return cgiter;
    }

    do {
        /* 1. M^{-1} res = out + M^{-1} in */
        /* 2. out + M^{-1} res -> out */
        /* 3. res = in - M.out */

        assign_sd2s(res_flt, res);
        zero_spinor_field_flt(out_flt);

        cgiter_flt += g5QMR_core_flt(&valid, par->err2_flt, par->max_iter_flt, M_flt, res_flt, out_flt);
        if (valid == 0) {
            lprintf("INVERTER", 20, "g5QMR_fltacc: using MINRES\n");
            MINRES_par Mpar;
            Mpar.err2 = par->err2_flt;
            Mpar.max_iter = 0;
            g5Herm_flt = M_flt;
            sh_flt = 0.;
            g5_spinor_field_flt(res_flt, res_flt); /* multiply input by g5 for MINRES */
            cgiter_minres += MINRES_flt(&Mpar, &Herm_flt, res_flt, out_flt, out_flt);
            g5_spinor_field_flt(res_flt, res_flt); /* restore input vector */
        }
#ifdef WITH_FUSE_MASTER_FOR
        _FUSE_MASTER_FOR(out->type, ix) {
            _FUSE_IDX(out->type, ix);
#else
        _MASTER_FOR(out->type, ix) {
#endif
            for (k = 0; k < 8 * NF; k++) {
                ((double *)_FIELD_AT(out, ix))[k] += (double)((float *)_FIELD_AT(out_flt, ix))[k];
            }
        }

        M(res, out);
        sub_spinor_field(res, in, res);
        err2 = sqnorm_spinor_field(res);
        cgiter++;

        lprintf("INVERTER", 20, "g5QMR_fltacc [%d]: res=%e\n", cgiter, err2 / innorm2);

    } while (err2 / innorm2 > par->err2);

    lprintf("INVERTER", 20, "g5QMR_fltacc: res=%e < %e \n", err2 / innorm2, par->err2);
    lprintf("INVERTER", 10, "g5QMR_fltacc: cgiter (g5QMR_flt,g5QMR,MINRES) = %d ; %d ; %d\n", cgiter_flt, cgiter,
            cgiter_minres);

    free_spinor_field(res);
    free_spinor_field_flt(in_flt);

    return cgiter_flt + cgiter_minres;
}
