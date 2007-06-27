#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include <malloc.h>

spinor_operator loc_H;
static suNf_spinor tmpspinor[VOLUME];

static float hmass;
static void H(suNf_spinor *out, suNf_spinor *in){
  g5Dphi(hmass,out,in);
}

static void D(suNf_spinor *out, suNf_spinor *in){
  Dphi(hmass,out,in);
}


/*prende spinori lunghi VOLUME/2 !!! */
static void D_pre(suNf_spinor *out, suNf_spinor *in){
  Dphi_(OE,tmpspinor,in);
  Dphi_(EO,out,tmpspinor);
  spinor_field_mul_add_assign_f(out,-hmass,in);
}

/*
 * Computes the matrix elements (H^-1)_{x,0}
 * H is the hermitean dirac operator
 * out is a vector of 4*NF spinor fields
 */
void quark_propagator_QMR(unsigned int source, int nm, float *mass, suNf_spinor_dble **out) {
  mshift_par QMR_par;
  int i;
  double *shift;
  suNf_spinor *in;

  /* allocate input spinor field */
  in = (suNf_spinor *)malloc(sizeof(suNf_spinor)*VOLUME/2);
  shift=(double*)malloc(sizeof(double)*(nm-1));

  set_spinor_len(VOLUME/2);

  /* the source is on the first even site */
  spinor_field_zero_f(in);
  *(((float *) in)+2*source)=1.; /* put in source */

  hmass=4.+mass[0];
  hmass*=hmass;

  QMR_par.n = nm;
  for(i=0;i<nm-1;++i){
    shift[i]=hmass-(4.+mass[i+1])*(4.+mass[i+1]);
  }
  QMR_par.shift = shift;
  QMR_par.err2 = 1.e-9;
  QMR_par.max_iter = 0;

  g5QMR_mshift(&QMR_par, &D_pre, in, out);

  /* compute res */
  for(i=0;i<nm;++i){
    assign_sd2s(VOLUME/2,tmpspinor,out[i]);
    Dphi_(OE,(suNf_spinor*)(out[i]+(VOLUME/2)),tmpspinor);
    assign_s2sd(VOLUME/2,out[i]+(VOLUME/2),(suNf_spinor*)(out[i]+(VOLUME/2)));
    spinor_field_mul_dble_f(out[i],-(4.+mass[i]),out[i]);
    spinor_field_minus_dble_f(out[i]+(VOLUME/2),out[i]+(VOLUME/2));
  }
  
  set_spinor_len(VOLUME);

  /* free input spinor field */
  free(in);
  free(shift);

}

/*
 * Computes the matrix elements (D^-1)_{x,0}
 * H is the hermitean dirac operator
 * out is a vector of 4*NF spinor fields
 */
void quark_propagator(unsigned int source, int nm, float *mass, suNf_spinor **out) {
  static MINRES_par MINRESpar;
  int i;
  suNf_spinor *in;

  /* allocate input spinor field */
  in = (suNf_spinor *)malloc(sizeof(suNf_spinor)*VOLUME);

  set_spinor_len(VOLUME);

  /* the source is on the first even site */
  spinor_field_zero_f(in);
  *(((float *) in)+2*source)=1.; /* put in source */

  hmass=mass[0];

  MINRESpar.spinorlen=VOLUME;
  MINRESpar.err2 = 1.e-10;
  MINRESpar.max_iter = 0;

  MINRES(&MINRESpar, &H, in, out[0],0);
  for(i=1;i<nm;++i){
    hmass=mass[0]-mass[i];
    MINRES(&MINRESpar, &H, in, out[i],out[i-1]);
  }

  /* free input spinor field */
  free(in);

}


