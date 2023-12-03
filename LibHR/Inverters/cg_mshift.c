
/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "libhr_core.h"
#include "random.h"
#include "io.h"
#include "memory.h"
#include "utils.h"
#include <math.h>
#include <assert.h>

/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
static int cg_mshift_core(short int *sflags, mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out) {
    spinor_field *k, *r, *Mk;
    spinor_field *p;
    double omega, oldomega, gamma;
    double alpha, lambda, delta;
    double innorm2;
    double *z1, *z2, *z3;

    int i;
    int cgiter;
    unsigned short notconverged;

    /* fare qualche check sugli input */
    assert(par->n > 0);
#ifndef CHECK_SPINOR_MATCHING
    for (i = 0; i < par->n; ++i) {
        _TWO_SPINORS_MATCHING(in, &out[i]);
    }
#endif

    /* allocate spinors fields and aux real variables */
    p = alloc_spinor_field(3 + par->n, in->type);
    k = p + par->n;
    r = k + 1;
    Mk = r + 1;

    z1 = malloc(sizeof(*z1) * (par->n));
    z2 = malloc(sizeof(*z2) * (par->n));
    z3 = malloc(sizeof(*z3) * (par->n));

    /* init recursion */
    cgiter = 0;
    omega = 1.;
    gamma = 0.;
    innorm2 = sqnorm_spinor_field(in);
    if (par->n == 1) { /* non multishift case */
        /* use out[0] as initial guess */
        M(Mk, &out[0]);
        ++cgiter;
        mul_add_assign_spinor_field(Mk, -par->shift[0], &out[0]);
        sub_spinor_field(r, in, Mk);

    } else { /* initial guess = 0 for multishift */
        copy_spinor_field(r, in);
    }
    copy_spinor_field(k, r);
    delta = sqnorm_spinor_field(r);
    for (i = 0; i < (par->n); ++i) {
        z1[i] = z2[i] = 1.;
        copy_spinor_field(&p[i], r);
        if (par->n != 1) { zero_spinor_field(&out[i]); }
        sflags[i] = 1;
    }

    /* cg recursion */
    do {
        M(Mk, k);
        alpha = prod_re_spinor_field(k, Mk);
        oldomega = omega;
        omega = -delta / alpha;
        for (i = 0; i < (par->n); ++i) {
            if (sflags[i]) {
                z3[i] = oldomega * z1[i] * z2[i] /
                        (omega * gamma * (z1[i] - z2[i]) + z1[i] * oldomega * (1. + par->shift[i] * omega));
                mul_add_assign_spinor_field(&out[i], -omega * z3[i] / z2[i], &p[i]);
            }
        }
        mul_add_assign_spinor_field(r, omega, Mk);
        lambda = sqnorm_spinor_field(r);
        gamma = lambda / delta;
        delta = lambda;

        mul_spinor_field(k, gamma, k);
        add_assign_spinor_field(k, r);
        notconverged = 0; /* assume that all vectors have converged */
        for (i = 0; i < (par->n); ++i) {
            /* check convergence of vectors */
            if (delta * z3[i] * z3[i] < par->err2 * innorm2) { sflags[i] = 0; }
            if (sflags[i]) {
                notconverged++;
                mul_spinor_field(&p[i], gamma * z3[i] * z3[i] / (z2[i] * z2[i]), &p[i]);
                mul_add_assign_spinor_field(&p[i], z3[i], r);
                z1[i] = z2[i];
                z2[i] = z3[i];
            }
        }

        /* Uncomment this to print cg recursion parameters 
       lprintf("CGTEST",0,"[ %d ] alpha=%e\n",cgiter,alpha);
       lprintf("CGTEST",0,"[ %d ] omega=%e\n",cgiter,omega);
       lprintf("CGTEST",0,"[ %d ] still runnning=%d\n",cgiter,notconverged);
	   for (i=0;i<par->n;++i) lprintf("CGTEST",0,"z3[%d]=%e; test=[%e,%e]\n",i,z3[i],delta*z3[i]*z3[i],par->err2*innorm2);
       lprintf("CGTEST",0,"\n[ %d ] gamma=%e\n",cgiter,gamma);
       lprintf("CGTEST",0,"[ %d ] delta=%e\n",cgiter,delta);
    */

        ++cgiter;
    } while ((par->max_iter == 0 || cgiter < par->max_iter) && notconverged);

    /* test results */
    for (i = 0; i < par->n; ++i) {
        double norm;
        M(Mk, &out[i]);
        ++cgiter;
        mul_add_assign_spinor_field(Mk, -par->shift[i], &out[i]);
        sub_spinor_field(Mk, Mk, in);
        norm = sqnorm_spinor_field(Mk) / sqnorm_spinor_field(in);
        sflags[i] = 1;
        if (fabs(norm) > par->err2) {
            sflags[i] = 0;
            lprintf("INVERTER", 30, "CG failed on vect %d: err2 = %1.8e > %1.8e\n", i, norm, par->err2);
        } else {
            lprintf("INVERTER", 20, "CG inversion: err2 = %1.8e < %1.8e\n", norm, par->err2);
        }
    }

    /* free memory */
    free_spinor_field(p);
    free(z1);
    free(z2);
    free(z3);

    /* return number of cg iter */
    return cgiter;
}

int cg_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out) {
#ifdef TIMING
#ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
#endif
    Timer clock;
    timer_set(&clock);
#endif

    //  int ix = ipt(0,0,0,0);
    //  printf( "test in 0 0 0 0 re=%f im=%f\n",   creal((*_FIELD_AT(in,ix)).c[0].c[0]), cimag( (*_FIELD_AT(in,ix)).c[0].c[0]   ) );

    int cgiter, msiter;
    int i;
    mshift_par par_save = *par;
    short int sflags[par->n];

    cgiter = cg_mshift_core(sflags, par, M, in, out);
    msiter = cgiter; /* save multishift iterations for logging */

    par->n = 1;
    for (i = 0; i < par_save.n; ++i) {
        int rep = 0;
        while (sflags[i] == 0) {
            par->shift = par_save.shift + i;
            cgiter += cg_mshift_core(sflags + i, par, M, in, out + i);
            if ((++rep) % 5 == 0) { lprintf("INVERTER", -10, "CG_mshift recursion = %d (precision too high?)\n", rep); }
        }
    }

    *par = par_save;

    lprintf("INVERTER", 10, "CG_mshift: MVM = %d/%d\n", msiter, cgiter);

#ifdef TIMING
#ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
#endif
    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("TIMING", 0, "cg_mshift %.6f s\n", elapsed_sec);
#endif

    return cgiter;
}
