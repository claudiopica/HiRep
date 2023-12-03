/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "libhr_core.h"
#include <math.h>

#define _print_par(c) //  printf("[%d] " #c " = %e\n",cgiter,c)

/*
 * performs the multi-shifted BiCGstab inversion for the hermitean matrix M:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
int HBiCGstab_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out) {
    spinor_field *s;
    spinor_field *r, *r1, *o, *Ms, *Mo, *o0;
    spinor_field *sptmp;

    double delta, phi;
    double *z1, *z2, *z3, *alpha, *beta, *chi, *rho;
    double rtmp1, rtmp2, rtmp3, oldalpha, oldchi;

    int i;
    int cgiter;
    char *sflags;
    unsigned short notconverged;

    /* fare qualche check sugli input */
    /* par->n deve essere almeno 2! */
    /*
    printf("numero vettori n=%d\n",par->n);
    for (i=0; i<(par->n); ++i) {
    printf("shift[%d]=%f\n",i,par->shift[i]);
    printf("out[%d]=%p\n",i,out[i]);      
    }
  */

    /* allocate spinors fields and aux real variables */
    /* implementation note: to minimize the number of malloc calls
   * objects of the same type are allocated together
   */
#ifndef CHECK_SPINOR_MATCHING
    for (i = 0; i < par->n; ++i) {
        _TWO_SPINORS_MATCHING(in, &out[i]);
    }
#endif

    s = alloc_spinor_field((par->n) + 6, in->type);
    r = s + par->n;
    r1 = r + 1;
    o = r1 + 1;
    Ms = o + 1;
    Mo = Ms + 1;
    o0 = Mo + 1;

    z1 = (double *)malloc(sizeof(double) * 7 * (par->n));
    z2 = z1 + (par->n);
    z3 = z2 + (par->n);
    alpha = z3 + (par->n);
    beta = alpha + (par->n);
    chi = beta + (par->n);
    rho = chi + (par->n);

    sflags = (char *)malloc(sizeof(char) * (par->n));

    /* init recursion */
    cgiter = 0;
    notconverged = 1;
    beta[0] = 1.;
    copy_spinor_field(r, in);
    for (i = 0; i < (par->n); ++i) {
        rho[i] = z1[i] = z2[i] = beta[0];
        alpha[i] = 0.;
        copy_spinor_field(&s[i], in);
        zero_spinor_field(&out[i]);
        sflags[i] = 1;
    }
    chi[0] = 0;
    /* choose omega so that delta and phi are not zero */
    /* copy_spinor_field(o, in);  omega = in ; this may be changed */
    /* gaussian_spinor_field(o0); */
    copy_spinor_field(o0, in);
    copy_spinor_field(o, o0);
    delta = prod_re_spinor_field(o, r);

    M(Ms, &s[0]);
    phi = prod_re_spinor_field(o0, Ms) / delta; /* o = in (see above) */

    /*  _print_par(delta);
  _print_par(phi);
  _print_par(alpha[0]);
  _print_par(beta[0]);*/

    /* cg recursion */
    do {
        ++cgiter;

        /* compute beta[0] and store in ctmp1 the old value (needed in the shifted update */
        rtmp1 = beta[0];
        beta[0] = -1. / phi; /* b=1/phi */

        /* compute omega and chi[0] */
        mul_spinor_field(o, beta[0], Ms);
        add_assign_spinor_field(o, r);

        M(Mo, o);

        oldchi = chi[0];
        chi[0] = prod_re_spinor_field(Mo, o) / sqnorm_spinor_field(Mo);

        /* compute r1 */
        mul_spinor_field(r1, chi[0], Mo);
        sub_spinor_field(r1, o, r1);

        /* update delta and alpha[0] */
        oldalpha = alpha[0];
        alpha[0] = -beta[0] / (chi[0] * delta);
        delta = prod_re_spinor_field(o0, r1); /* in = omega al passo zero */
        alpha[0] *= delta;

        /* compute new out[0] */
        rtmp2 = -beta[0];
        lc_add_assign_spinor_field(&out[0], rtmp2, &s[0], chi[0], o);

        /* compute new s[0] */
        rtmp3 = -alpha[0] * chi[0];
        lc_spinor_field(Mo, alpha[0], &s[0], rtmp3, Ms); /* use Mo as temporary storage */
        add_spinor_field(&s[0], r1, Mo);

        /* assign r<-r1 */
        sptmp = r;
        r = r1;
        r1 = sptmp; /*exchange pointers */

        /* update phi */
        M(Ms, &s[0]);
        phi = prod_re_spinor_field(o0, Ms) / delta;

        /*    _print_par(delta);
    _print_par(phi);
    _print_par(alpha[0]);
    _print_par(beta[0]);
    _print_par(chi[0]);*/

        for (i = 1; i < (par->n); ++i) { /* update shifted quantities i=1..n-1 */
            if (sflags[i]) {
                /* compute z3[i]/z2[i] */
                z3[i] =
                    z1[i] * rtmp1 / (beta[0] * oldalpha * (z1[i] - z2[i]) + z1[i] * rtmp1 * (1. + par->shift[i - 1] * beta[0]));
                _print_par(z3[i]);
                /*compute beta[i] */
                beta[i] = beta[0] * z3[i];
                _print_par(beta[i]);
                /*compute chi[i] */
                chi[i] = chi[0] / (1. - chi[0] * par->shift[i - 1]);
                _print_par(chi[i]);
                /* compute alpha[i] */
                alpha[i] = alpha[0] * z3[i] * z3[i];
                _print_par(alpha[i]);
                /* compute z3[i] */
                z3[i] *= z2[i];
                _print_par(z3[i]);
                /* update solution */
                rtmp2 = -beta[i];
                rtmp3 = chi[i] * rho[i] * z3[i];
                lc_add_assign_spinor_field(&out[i], rtmp2, &s[i], rtmp3, o);
                /* update s[i] */
                rtmp2 = chi[i] / beta[i] * rho[i];
                rtmp3 = rtmp2 * z2[i];
                rtmp2 *= -z3[i];
                lc_add_assign_spinor_field(&s[i], rtmp2, o, rtmp3, r1); /* not done yet */

                rho[i] /= (1. - rho[0] * par->shift[i - 1]); /* update rho */
                _print_par(rho[i]);

                rtmp2 = z3[i] * rho[i];
                lc_spinor_field(Mo, rtmp2, r, alpha[i], &s[i]); /* use Mo as temporary storage */
                copy_spinor_field(&s[i], Mo);

                /* change pointers instead */
                /*
	  sptmp=s[i];
	  s[i]=Mo;
	  Mo=sptmp;
	*/

                /* shift z2<-z3; z1<-z2 */
                z2[i] = z3[i];
                z1[i] = z2[i];
            }
        }

        rtmp1 = sqnorm_spinor_field(r);

        if (rtmp1 < par->err2) { notconverged = 0; }

        //    printf("[ %d ] residuo=%e\n",cgiter,rtmp1);

        /* Uncomment this to print cg recursion parameters
       printf("[ %d ] alpha=%e\n",cgiter,alpha);
       printf("[ %d ] omega=%e\n",cgiter,omega);
       printf("[ %d ] still runnning=%d\n",cgiter,notconverged);
       for (i=0;i<par->n;++i) printf("z3[%d]=%e; ",i,z3[i]);
       printf("\n[ %d ] gamma=%e\n",cgiter,gamma);
       printf("[ %d ] delta=%e\n",cgiter,delta);
    */

    } while ((par->max_iter == 0 || cgiter < par->max_iter) && notconverged);

    /* test results */
#ifndef NDEBUG
    for (i = 0; i < par->n; ++i) {
        double norm;
        M(Ms, &out[i]);
        if (i != 0) { mul_add_assign_spinor_field(Ms, -par->shift[i - 1], &out[i]); }
        mul_add_assign_spinor_field(Ms, -1.0, in);
        norm = sqnorm_spinor_field(Ms);
        if (fabs(norm) > 5. * par->err2) { printf("BiCGstab Failed: err2[%d] = %e\n", i, norm); }
    }
#endif

    /* free memory */
    free_spinor_field(s);
    free(z1);
    free(sflags);

    /* return number of cg iter */
    return cgiter;
}
