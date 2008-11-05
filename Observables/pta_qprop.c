/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"



static double hmass, hmass_pre;

static void H(spinor_field *out, spinor_field *in){
  g5Dphi(hmass,out,in);
}

static void D(spinor_field *out, spinor_field *in){
  Dphi(hmass,out,in);
}

static void g5D_pre(spinor_field *out, spinor_field *in){
  g5Dphi_eopre(hmass_pre,out,in);
}

/*
 * Computes the matrix elements (H^-1)_{x,0}
 * using the QMR inverter with even/odd preconditioning
 * In order to use this inverter without look-ahead on point
 * like sources we add a background noise to the source which is then 
 * subtracted at the end at the cost of one more inversion.
 */
void pta_qprop_QMR_eo(spinor_field **pta_qprop, int nm, double *mass, double acc) {
	mshift_par QMR_par;
	int i, x0;
	int source;
	double *shift;
	spinor_field *in=0;
	spinor_field *resdn=0;
	spinor_field *resd=0;
	spinor_field *res=0;
	spinor_field qprop_mask;
	int cgiter=0;
	double norm;

#ifdef NDEBUG
  spinor_field *test=0;
  test=alloc_spinor_field_f(1,&glattice);
  spinor_field *test_e=0;
  test_e=alloc_spinor_field_f(1,&glat_even);
#endif

	res = alloc_spinor_field_f(1,&glattice);

	in=alloc_spinor_field_f(2*nm+1,&glat_even);
	resdn=in+1;
	resd=resdn+nm;

	/* set up inverters parameters */
	shift=(double*)malloc(sizeof(double)*(nm));
	hmass_pre=mass[0]; /* we can put any number here!!! */
	for(i=0;i<nm;++i){
		shift[i]=(4.+mass[i])*(4.+mass[i])-(4.+hmass_pre)*(4.+hmass_pre);
	}
	QMR_par.n = nm;
	QMR_par.shift = shift;
	QMR_par.err2 = .5*acc;
	QMR_par.max_iter = 0;

  x0=0;
  if(glat_even.master_end[0]<glat_even.master_start[0])
    x0 = glat_even.master_start[1];
  else
    x0 = glat_even.master_start[0];

	/* noisy background */
	gaussian_spinor_field(in);
	if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0) {
    for (source=0;source<NF*4*2;++source)
	  	*( (double *)_FIELD_AT(in,x0) + source ) =0.; /* zero in source */
	}
	norm=sqrt(spinor_field_sqnorm_f(in));
	spinor_field_mul_f(in,1./norm,in);

	/* invert noise */
	cgiter+=g5QMR_mshift(&QMR_par, &g5D_pre, in, resdn);

	/* now loop over sources */
	for (source=0;source<4*NF;++source){
		/* put in source on an EVEN site */
	  if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0)
	    *( (double *)_FIELD_AT(in,x0) + 2*source ) =1.;

		cgiter+=g5QMR_mshift(&QMR_par, &g5D_pre, in, resd);

		for(i=0;i<QMR_par.n;++i){

#ifdef NDEBUG
			/* this is a test of the solution */
			hmass_pre=mass[i];
			g5D_pre(test_e,&resd[i]);
			++cgiter;
			spinor_field_sub_f(test_e,test_e,in);
			norm=spinor_field_sqnorm_f(test_e);
			lprintf("PROPAGATOR",0,"g5QMR_eo residuum of source [%d,%d] = %e\n",i,source,norm);
			hmass_pre=mass[0];
#endif /* NDEBUG */

			spinor_field_sub_f(&resd[i],&resd[i],&resdn[i]); /* compute difference */

			/* compute solution */
			qprop_mask=pta_qprop[i][source];
			qprop_mask.type=&glat_even;
			spinor_field_mul_f(&qprop_mask,(4.+mass[i]),&resd[i]);
			qprop_mask.type=&glat_odd;
			Dphi_(&qprop_mask,&resd[i]);
      spinor_field_minus_f(&qprop_mask,&qprop_mask);
			if(source&1) ++cgiter; /* count only half of calls. works because the number of sources is even */

#ifdef NDEBUG
			/* this is a test of the solution */
			hmass=mass[i];
			H(test,&pta_qprop[i][source]);
			++cgiter;
			if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0)
			  *( (double *)_FIELD_AT(test,x0) + 2*source ) -=1.;
			norm=spinor_field_sqnorm_f(test);
			lprintf("PROPAGATOR",0,"g5QMR_eo residuum of source [%d,%d] = %e\n",i,source,norm);
			hmass=mass[0];
#endif /* NDEBUG */
		}

		/* remove source */
		if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0)
		  *( (double *)_FIELD_AT(in,x0) + 2*source ) =0.;
	}

	lprintf("PROPAGATOR",10,"QMR_eo MVM = %d\n",cgiter);

	/* free memory */

	free_spinor_field(in);
	free_spinor_field(res);
	free(shift);
#ifdef NDEBUG
	free_spinor_field(test);
	free_spinor_field(test_e);
#endif
}



/*
 * Computes the matrix elements (H^-1)_{x,0}
 * using the QMR inverter
 * In order to use this inverter without look-ahead on point
 * like sources we add a background noise to the source which is then 
 * subtracted at the end at the cost of one more inversion.
 */
void pta_qprop_QMR(spinor_field **pta_qprop, int nm, double *mass, double acc) {
  mshift_par QMR_par;
  int i, x0;
	int source;
  double *shift;
  spinor_field *in=0;
  spinor_field *resdn=0;
  spinor_field *resd=0;
	int cgiter=0;
	double norm;

#ifdef NDEBUG
  spinor_field *test=0;
  test=alloc_spinor_field_f(1,&glattice);
#endif

  /* allocate input spinor field */
	in=alloc_spinor_field_f(1+2*nm,&glattice);
	resdn=in+1;
	resd=resdn+nm;

	/* set up inverters parameters */
  shift=(double*)malloc(sizeof(double)*(nm));
  hmass=mass[0]; /* we can put any number here!!! */
  for(i=0;i<nm;++i){
    shift[i]=hmass-mass[i];
  }
  QMR_par.n = nm;
  QMR_par.shift = shift;
  QMR_par.err2 = .5*acc;
  QMR_par.max_iter = 0;

  /* noisy background */
	gaussian_spinor_field(in);
  x0 = ipt(0,0,0,0);
	if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0) {
    for (source=0;source<NF*4*2;++source)
	  	*( (double *)_FIELD_AT(in,x0) + source ) =0.; /* zero in source */
	}
	norm=sqrt(spinor_field_sqnorm_f(in));
	spinor_field_mul_f(in,1./norm,in);
	
	/* invert noise */
	spinor_field_g5_f(in,in);
  cgiter+=g5QMR_mshift(&QMR_par, &D, in, resdn);
  spinor_field_g5_f(in,in);

	/* now loop over sources */
	for (source=0;source<4*NF;++source){
		/* put in source */
	  if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0)
	    *( (double *)_FIELD_AT(in,x0) + 2*source ) =1.;
		spinor_field_g5_f(in,in);

		cgiter+=g5QMR_mshift(&QMR_par, &D, in, resd);

		for(i=0;i<QMR_par.n;++i){

#ifdef NDEBUG
			/* this is a test of the solution */
			hmass=mass[i];
			D(test,&resd[i]);
			++cgiter;
			spinor_field_sub_f(test,test,in);
			norm=spinor_field_sqnorm_f(test);
			lprintf("PROPAGATOR",0,"g5QMR residuum of source [%d,%d] = %e\n",i,source,norm);
			hmass=mass[0];
#endif /* NDEBUG */

			spinor_field_sub_f(&pta_qprop[i][source],&resd[i],&resdn[i]); /* compute difference */

#ifdef NDEBUG
			/* this is a test of the solution */
			hmass=mass[i];
			H(test,&pta_qprop[i][source]);
			++cgiter;
			if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0)
			  *( (double *)_FIELD_AT(test,x0) + 2*source ) -=1.;
			norm=spinor_field_sqnorm_f(test);
			lprintf("PROPAGATOR",0,"g5QMR residuum of source [%d,%d] = %e\n",i,source,norm);
			hmass=mass[0];
#endif /* NDEBUG */
		}

		/* remove source */
		spinor_field_g5_f(in,in);
		if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0)
		  *( (double *)_FIELD_AT(in,x0) + 2*source ) =0.;
	}

	lprintf("PROPAGATOR",10,"QMR MVM = %d\n",cgiter);
  
  /* free memory */
	free_spinor_field(in);
  free(shift);
#ifdef NDEBUG
	free_spinor_field(test);
#endif
}


/*
 * Computes the matrix elements (H^-1)_{x,0}
 * H is the hermitean dirac operator
 * out is a vector of nm spinor fields
 */
void pta_qprop_MINRES(spinor_field **pta_qprop, int nm, double *mass, double acc) {
  static MINRES_par MINRESpar;
  int i,x0,cgiter,source;
  spinor_field *in;
#ifdef NDEBUG
  double norm;
#endif

  /* allocate input spinor field */
  in = alloc_spinor_field_f(1,&glattice);


  /* the source is on the first even site */
  spinor_field_zero_f(in);
  x0 = ipt(0,0,0,0);

  cgiter=0;

	for (source=0;source<4*NF;++source){
		/* put in source */
	  if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0)
	    *( (double *)_FIELD_AT(in,x0) + 2*source ) =1.;

#ifdef NDEBUG
    norm=spinor_field_sqnorm_f(in);
    lprintf("PROPAGATOR",0,"norm of source [%d] = %e\n",source,norm);
#endif

    hmass=mass[0];

    MINRESpar.err2 = acc;
    MINRESpar.max_iter = 0;

    cgiter+=MINRES(&MINRESpar, &H, in, &pta_qprop[0][source],0);
    for(i=1;i<nm;++i){
      hmass=mass[i];
      cgiter+=MINRES(&MINRESpar, &H, in, &pta_qprop[i][source],&pta_qprop[i-1][source]);
    }

		/* remove source */
		if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0)
		  *( (double *)_FIELD_AT(in,x0) + 2*source ) =0.;
  }
  
  lprintf("PROPAGATOR",10,"MINRES MVM = %d",cgiter);

  /* free input spinor field */
  free_spinor_field(in);

}
