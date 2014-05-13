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
#include "io.h"



static double hmass, hmass_pre;

static void H_pta(spinor_field *out, spinor_field *in){
  g5Dphi(hmass,out,in);
}

static void D_pta(spinor_field *out, spinor_field *in){
  Dphi(hmass,out,in);
}

static void D_pta_pre(spinor_field *out, spinor_field *in){
  Dphi_eopre(hmass_pre,out,in);
}

/*
 * Computes the matrix elements (H^-1)_{x,x0}
 * using the QMR inverter with even/odd preconditioning
 * In order to use this inverter without look-ahead on point
 * like sources we add a background noise to the source which is then 
 * subtracted at the end at the cost of one more inversion.
 */
void pta_qprop_QMR_eo(int g0[4], spinor_field **pta_qprop, int nm, double *mass, double acc) {
	mshift_par QMR_par;
	int i, x0, C0[4], c0[4];
	int source;
	double *shift;
	spinor_field *in=0;
	spinor_field *resdn=0;
	spinor_field *resd=0;
	spinor_field *res=0;
	spinor_field qprop_mask;
	int cgiter=0;
	double norm;

#ifndef NDEBUG
	spinor_field *test_e=0;
	spinor_field *test=0;
	test=alloc_spinor_field_f(1,&glattice);
	test_e=alloc_spinor_field_f(1,&glat_even);
#endif


	res = alloc_spinor_field_f(1,&glattice);

	in=alloc_spinor_field_f(2*nm+1,&glat_even);
	resdn=in+1;
	resd=resdn+nm;

	C0[0]=g0[0]/T; C0[1]=g0[1]/X; C0[2]=g0[2]/Y; C0[3]=g0[3]/Z;
	c0[0]=g0[0]%T; c0[1]=g0[1]%X; c0[2]=g0[2]%Y; c0[3]=g0[3]%Z;
	x0=ipt(c0[0],c0[1],c0[2],c0[3]);
	
	/* set up inverters parameters */
	shift=(double*)malloc(sizeof(double)*(nm));
	hmass_pre=mass[0]; /* we can put any number here!!! */
	for(i=0;i<nm;++i){
		shift[i]=(4.+hmass_pre)*(4.+hmass_pre)-(4.+mass[i])*(4.+mass[i]);
	}
	QMR_par.n = nm;
	QMR_par.shift = shift;
	QMR_par.err2 = .5*acc;
	QMR_par.max_iter = 0;

	/* noisy background */
	gaussian_spinor_field(in);
	if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3]) {
		for (source=0;source<NF*4*2;++source) {
			*( (double *)_FIELD_AT(in,x0) + source ) =0.; /* zero in source */
		}
	}
	norm=sqrt(spinor_field_sqnorm_f(in));
	spinor_field_mul_f(in,1./norm,in);

#ifndef NDEBUG
	spinor_field_zero_f(test_e);
	if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3]) {
		for (source=0;source<NF*4*2;++source) {
			*( (double *)_FIELD_AT(test_e,x0) + source ) =1.;
		}
	}
	norm=sqrt(spinor_field_sqnorm_f(test_e));
	error(norm==0.,1,"pta_qprop_QMR_eo [pta_qprop.c]",
        "The origin is an odd point.");
#endif /* NDEBUG */
  
	/* invert noise */
	for(i=0;i<QMR_par.n;++i){
#ifndef NDEBUG
		norm=spinor_field_sqnorm_f(&resdn[i]);
		lprintf("PROPAGATOR",0,"Sqnorm resdn[%d] = %e\n",i,norm);
#endif /* NDEBUG */
		spinor_field_zero_f(&resdn[i]);
	}
	spinor_field_g5_assign_f(in);
	cgiter+=g5QMR_mshift(&QMR_par, &D_pta_pre, in, resdn);
	spinor_field_g5_assign_f(in);

	/* now loop over sources */
	for (source=0;source<4*NF;++source){
		/* put in source on an EVEN site */
		if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3])
			*( (double *)_FIELD_AT(in,x0) + 2*source ) =1.;
		spinor_field_g5_assign_f(in);

		for(i=0;i<QMR_par.n;++i){
#ifndef NDEBUG
	  		norm=spinor_field_sqnorm_f(&resd[i]);
			lprintf("PROPAGATOR",0,"Sqnorm resd[%d] = %e\n",i,norm);
#endif /* NDEBUG */
  	  		spinor_field_zero_f(&resd[i]);
		}
		cgiter+=g5QMR_mshift(&QMR_par, &D_pta_pre, in, resd);

		for(i=0;i<QMR_par.n;++i){

#ifndef NDEBUG
			/* this is a test of the inverter on the difference vector */
			hmass_pre=mass[i];
			D_pta_pre(test_e,&resd[i]);
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
			/*this below is a formal operation that needs to be done every time we change type of spinor
			 but since the shift for the even spinor is always zero we can skip it*/
			/* qprop_mask.ptr=pta_qprop[i][source].ptr+glat_even.master_shift; */


			spinor_field_mul_f(&qprop_mask,(4.+mass[i]),&resd[i]);
			qprop_mask.type=&glat_odd;
			qprop_mask.ptr=pta_qprop[i][source].ptr+glat_odd.master_shift;
			
			Dphi_(&qprop_mask,&resd[i]);
			spinor_field_minus_f(&qprop_mask,&qprop_mask);
			if(source&1) ++cgiter; /* count only half of calls. works because the number of sources is even */

#ifndef NDEBUG
			/* this is a test of the solution */
			hmass=mass[i];
			H_pta(test,&pta_qprop[i][source]);
			++cgiter;
			if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3])
			  *( (double *)_FIELD_AT(test,x0) + 2*source ) -=1.;
			norm=spinor_field_sqnorm_f(test);
			lprintf("PROPAGATOR",0,"g5QMR_eo residuum of source [%d,%d] = %e\n",i,source,norm);
			hmass=mass[0];
#endif /* NDEBUG */
		}

		/* remove source */
		spinor_field_g5_assign_f(in);
		if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3])
			*( (double *)_FIELD_AT(in,x0) + 2*source ) =0.;
	}

	lprintf("PROPAGATOR",10,"QMR_eo MVM = %d\n",cgiter);

	/* free memory */

	free_spinor_field_f(in);
	free_spinor_field_f(res);
	free(shift);
#ifndef NDEBUG
	free_spinor_field_f(test_e);
	free_spinor_field_f(test);
#endif
}


/*
 * Computes the matrix elements (H^-1)_{x,0}
 * using the QMR inverter
 * In order to use this inverter without look-ahead on point
 * like sources we add a background noise to the source which is then 
 * subtracted at the end at the cost of one more inversion.
 */
void pta_qprop_QMR(int g0[4], spinor_field **pta_qprop, int nm, double *mass, double acc) {
  mshift_par QMR_par;
	int i, x0, C0[4], c0[4];
	int source;
  double *shift;
  spinor_field *in=0;
  spinor_field *resdn=0;
  spinor_field *resd=0;
	int cgiter=0;
	double norm;

#ifndef NDEBUG
  spinor_field *test=0;
  test=alloc_spinor_field_f(1,&glattice);
#endif

  /* allocate input spinor field */
	in=alloc_spinor_field_f(1+2*nm,&glattice);
	resdn=in+1;
	resd=resdn+nm;

	C0[0]=g0[0]/T; C0[1]=g0[1]/X; C0[2]=g0[2]/Y; C0[3]=g0[3]/Z;
	c0[0]=g0[0]%T; c0[1]=g0[1]%X; c0[2]=g0[2]%Y; c0[3]=g0[3]%Z;
	x0=ipt(c0[0],c0[1],c0[2],c0[3]);

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
	if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3]) {
    for (source=0;source<NF*4*2;++source)
	  	*( (double *)_FIELD_AT(in,x0) + source ) =0.; /* zero in source */
	}
	norm=sqrt(spinor_field_sqnorm_f(in));
	spinor_field_mul_f(in,1./norm,in);
	
	/* invert noise */
	spinor_field_g5_assign_f(in);
  cgiter+=g5QMR_mshift(&QMR_par, &D_pta, in, resdn);
  spinor_field_g5_assign_f(in);

	/* now loop over sources */
	for (source=0;source<4*NF;++source){
		/* put in source */
	  if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3])
	    *( (double *)_FIELD_AT(in,x0) + 2*source ) =1.;
		spinor_field_g5_assign_f(in);

		cgiter+=g5QMR_mshift(&QMR_par, &D_pta, in, resd);

		for(i=0;i<QMR_par.n;++i){

#ifndef NDEBUG
			/* this is a test of the inverter on the difference vector */
			hmass=mass[i];
			D_pta(test,&resd[i]);
			++cgiter;
			spinor_field_sub_f(test,test,in);
			norm=spinor_field_sqnorm_f(test);
			lprintf("PROPAGATOR",0,"g5QMR residuum of source [%d,%d] = %e\n",i,source,norm);
			hmass=mass[0];
#endif /* NDEBUG */

			spinor_field_sub_f(&pta_qprop[i][source],&resd[i],&resdn[i]); /* compute difference */

#ifndef NDEBUG
			/* this is a test of the solution */
			hmass=mass[i];
			H_pta(test,&pta_qprop[i][source]);
			++cgiter;
    	if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3])
			  *( (double *)_FIELD_AT(test,x0) + 2*source ) -=1.;
			norm=spinor_field_sqnorm_f(test);
			lprintf("PROPAGATOR",0,"g5QMR residuum of source [%d,%d] = %e\n",i,source,norm);
			hmass=mass[0];
#endif /* NDEBUG */
		}

		/* remove source */
		spinor_field_g5_assign_f(in);
  	if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3])
		  *( (double *)_FIELD_AT(in,x0) + 2*source ) =0.;
	}

	lprintf("PROPAGATOR",10,"QMR MVM = %d\n",cgiter);
  
  /* free memory */
	free_spinor_field_f(in);
  free(shift);
#ifndef NDEBUG
	free_spinor_field_f(test);
#endif
}


/*
 * Computes the matrix elements (H^-1)_{x,0}
 * H is the hermitean dirac operator
 * out is a vector of nm spinor fields
 */
void pta_qprop_MINRES(int g0[4], spinor_field **pta_qprop, int nm, double *mass, double acc) {
  static MINRES_par MINRESpar;
	int i, x0, C0[4], c0[4];
  int cgiter,source;
  spinor_field *in;
#ifndef NDEBUG
  double norm;
#endif

  /* allocate input spinor field */
  in = alloc_spinor_field_f(1,&glattice);

  /* the source is on the first even site */
  spinor_field_zero_f(in);
  
	C0[0]=g0[0]/T; C0[1]=g0[1]/X; C0[2]=g0[2]/Y; C0[3]=g0[3]/Z;
	c0[0]=g0[0]%T; c0[1]=g0[1]%X; c0[2]=g0[2]%Y; c0[3]=g0[3]%Z;
	x0=ipt(c0[0],c0[1],c0[2],c0[3]);

  cgiter=0;

	for (source=0;source<4*NF;++source){
		/* put in source */
  	if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3])
	    *( (double *)_FIELD_AT(in,x0) + 2*source ) =1.;

#ifndef NDEBUG
    norm=spinor_field_sqnorm_f(in);
    lprintf("PROPAGATOR",0,"norm of source [%d] = %e\n",source,norm);
#endif

    hmass=mass[0];

    MINRESpar.err2 = acc;
    MINRESpar.max_iter = 0;

    cgiter+=MINRES(&MINRESpar, &H_pta, in, &pta_qprop[0][source],0);
    for(i=1;i<nm;++i){
      hmass=mass[i];
      cgiter+=MINRES(&MINRESpar, &H_pta, in, &pta_qprop[i][source],&pta_qprop[i-1][source]);
    }

		/* remove source */
	  if(COORD[0]==C0[0] && COORD[1]==C0[1] && COORD[2]==C0[2] && COORD[3]==C0[3])
		  *( (double *)_FIELD_AT(in,x0) + 2*source ) =0.;
  }
  
  lprintf("PROPAGATOR",10,"MINRES MVM = %d",cgiter);

  /* free input spinor field */
  free_spinor_field_f(in);

}
