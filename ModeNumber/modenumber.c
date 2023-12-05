#include "libhr.h"
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>

//local header
#include "modenumber.h"

#define NMAX 1000

static int init = 0;

static mshift_par par;

static double epsilon;
static double delta;
static int order;
static double c[NMAX + 1];
static double star;

static double mass = 0.;

static spinor_field *w0, *w1, *w2, *x, *eta;

static int nhits;

/*
static int nsources;
static double ratio;
*/

void init_modenumber(double m, double inv, int nh, char *approxfile) {
    error(init == 1, 1, "modenumber.c", "Already initialized!");

    mass = m;

    lprintf("MODENUMBER", 0, "Mass = %e\n", mass);

    nhits = nh;

    lprintf("MODENUMBER", 0, "Number of random spinors = %d\n", nhits);

    par.n = 1;
    par.shift = (double *)malloc(sizeof(double) * (par.n));
    par.err2 = inv;
    par.max_iter = 0;
    par.shift[0] = 0.;

    lprintf("MODENUMBER", 0, "Error2 cg_mshift = %e\n", par.err2);

    FILE *file = fopen(approxfile, "r");
    error(file == NULL, 1, "init_modenumber [modenumber.c]", "Failed to open approximation file\n");

    int ret;

    ret = fscanf(file, "%d", &order);
    error(ret == 0 || order <= 0 || feof(file), 1, "init_modenumber [modenumber.c]", "Wrong format for approximation file\n");
    lprintf("MODENUMBER", 0, "Chebychev approximation: order = %d\n", order);

    ret = fscanf(file, "%lf", &epsilon);
    error(ret == 0 || epsilon <= 0. || feof(file), 1, "init_modenumber [modenumber.c]",
          "Wrong format for approximation file\n");
    lprintf("MODENUMBER", 0, "Chebychev approximation: epsilon = %e\n", epsilon);

    ret = fscanf(file, "%lf", &delta);
    error(ret == 0 || delta <= 0. || feof(file), 1, "init_modenumber [modenumber.c]", "Wrong format for approximation file\n");
    lprintf("MODENUMBER", 0, "Chebychev approximation: delta = %e\n", delta);

    ret = fscanf(file, "%lf", &star);
    error(ret == 0 || star <= 0. || feof(file), 1, "init_modenumber [modenumber.c]", "Wrong format for approximation file\n");
    lprintf("MODENUMBER", 0, "M/Mstar = %e\n", sqrt(star));

    for (int i = 0; i <= order; i++) {
        ret = fscanf(file, "%lf", &c[i]);
        error(ret == 0 || feof(file), 1, "init_modenumber [modenumber.c]", "Wrong format for approximation file\n");
        lprintf("MODENUMBER", 0, "c[%d] = %e\n", i, c[i]);
    }

    fclose(file);

    /*  
  nsources = ns;

  lprintf("PROJCORRELATOR",0,"Number of sources = %d\n",nsources);
  
  ratio = r;

  lprintf("PROJCORRELATOR",0,"Ratio of eigenvalues for correlator = %e\n",ratio);
*/

    x = alloc_spinor_field(7, &glattice);
    w0 = x + 2;
    w1 = x + 3;
    w2 = x + 4;
    eta = x + 5;

    init = 1;
}

void free_modenumber() {
    error(init == 0, 1, "modenumber.c", "Not initialized!");
    free(par.shift);
    free_spinor_field(x);
    init = 0;
}

static void H2X(spinor_field *out, spinor_field *in) {
    g5Dphi_sq(mass, out, in);
}

/*
static double H2X_shift=0.;
static void g5H2X(spinor_field *out, spinor_field *in){
  g5Dphi_sq(mass, out, in);
  mul_add_assign_spinor_field(out,H2X_shift,in);
  g5_assign_spinor_field(out);
}
*/

static void operatorX(spinor_field *out, spinor_field *in, double M2) {
    par.shift[0] = -M2;
    cg_mshift(&par, &H2X, in, out);
    mul_spinor_field(out, -2. * M2, out);
    add_assign_spinor_field(out, in);
}

static void operatorX2(spinor_field *out, spinor_field *in, double M2) {
    /*H2X_shift = M2;*/
    par.shift[0] = -M2;

    zero_spinor_field(&x[0]);
    zero_spinor_field(&x[1]);

    cg_mshift(&par, &H2X, in, &x[0]);
    cg_mshift(&par, &H2X, &x[0], &x[1]);

    /* x[0] = H2X^{-1} in = g5H2X^{-1} g5 in */
    /* x[1] = H2X^{-2} in = g5H2X^{-1} g5 x[0] */
    /*
  g5_assign_spinor_field(in);
  g5QMR_mshift(&par, &g5H2X, in, &x[0]);
  g5_assign_spinor_field(in);
  g5_assign_spinor_field(&x[0]);  
  g5QMR_mshift(&par, &H2X, &x[0], &x[1]);
  g5_assign_spinor_field(&x[0]);
  */

    copy_spinor_field(out, in);
    lc_add_assign_spinor_field(out, -4. * M2, &x[0], 4. * M2 * M2, &x[1]);
}

static void operatorZ(spinor_field *out, spinor_field *in, double M2) {
    /*double z=(2.*x*x-1.-epsilon)/(1.-epsilon);*/
    operatorX2(out, in, M2);
    mul_spinor_field(out, 2. / (1. - epsilon), out);
    mul_add_assign_spinor_field(out, -(1. + epsilon) / (1. - epsilon), in);
}

static void operatorH(spinor_field *out, spinor_field *in, double M2) {
    spinor_field *tmp;

    /*
  double b0,b1,b2;
  double z=(2.*x*x-1.-epsilon)/(1.-epsilon);
  
  b0=b1=0.;
  for(int n=order; n>=0; n--) {
    b2=b1;
    b1=b0;
    b0 = c[n] + 2.*z*b1 - b2;
  }
  return .5 - .5*x*(b0 - b1*z);
  */

    mul_spinor_field(w1, c[order], in);

    operatorZ(w0, in, M2);
    mul_spinor_field(w0, 2. * c[order], w0);
    mul_add_assign_spinor_field(w0, c[order - 1], in);

    for (int n = order - 2; n >= 0; n--) {
        tmp = w2;
        w2 = w1;
        w1 = w0;
        w0 = tmp;
        operatorZ(w0, w1, M2);

        mul_spinor_field(w0, 2., w0);
        lc_add_assign_spinor_field(w0, c[n], in, -1., w2);
    }

    operatorZ(w2, w1, M2);
    sub_assign_spinor_field(w0, w2);
    operatorX(out, w0, M2);
    mul_spinor_field(out, -.5, out);

    mul_add_assign_spinor_field(out, .5, in);
}

double ModeNumber(double M2) {
    /*double norm;*/
    double ret = 0.;
    double M2star;

    error(nhits <= 0, 1, "modenumber.c", "[ModeNumber] nhits must be positive!");

    M2star = M2 / star;

    for (int i = 0; i < nhits; i++) {
        create_z2_volume_source(&eta[0]);

        operatorH(&eta[1], &eta[0], M2star);
        operatorH(&eta[0], &eta[1], M2star);

        ret += sqnorm_spinor_field(&eta[0]);
    }

    return ret / nhits;
}
