#include "global.h"
#include "linear_algebra.h"
#include "suN.h"
#include "observables.h"
#include "inverters.h"
#include <malloc.h>

/*
 * Computes the pi correlator
 * \sum_{\vec x} Tr(H^-1(\vec x, t, 0)^2 )
 */
void pi_correlator(float *out, suNf_spinor *qp) {
  int t,x,y,z;
  for (t=0; t<T; ++t) {
    double hc=0.;
    for (x=0; x<L; ++x) {
      for (y=0; y<L; ++y) {
	for (z=0; z<L; ++z) {
	  suNf_spinor *s = &(qp[ipt[t][x][y][z]]);
	  hc+=_spinor_prod_re_f(*s,*s);
	}
      }
    }
    out[t] = hc;
  }
}

void pi_correlator_QMR(float *out, suNf_spinor_dble *qp) {
  int t,x,y,z;
  for (t=0; t<T; ++t) {
    double hc=0.;
    for (x=0; x<L; ++x) {
      for (y=0; y<L; ++y) {
	for (z=0; z<L; ++z) {
	  suNf_spinor_dble *s = &(qp[ipt[t][x][y][z]]);
	  hc+=_spinor_prod_re_f(*s,*s);
	}
      }
    }
    out[t] = hc;
  }
}

/*
 * Computes the rho correlator
 * \sum_{\vec x} Tr([g_3 g_5 H^-1(\vec x, t, 0)]^2 )
 */
void rho_correlator(float *out, suNf_spinor **qp) {
  int t,x,y,z, i;
  for (t=0; t<T; ++t) {
    double hc=0.;
    for (x=0; x<L; ++x) {
      for (y=0; y<L; ++y) {
	for (z=0; z<L; ++z) {
	  for (i=0; i<4*NF; ++i) {
	    suNf_spinor *s = &(qp[i][ipt[t][x][y][z]]);/*MOLTIPLICARE PER G_3G_5 */
	    hc+=_spinor_prod_re_f(*s,*s);
	  }
	}
      }
    }
    out[t] = hc;
  }
}
