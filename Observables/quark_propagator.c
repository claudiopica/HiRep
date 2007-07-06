#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "update.h"
#include "error.h"
#include <malloc.h>
#include <math.h>

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
void old_quark_propagator_QMR(unsigned int source, int nm, float *mass, suNf_spinor_dble **out) {
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
 * Computes the matrix elements (H^-1)_{x,0}
 * H is the hermitean dirac operator
 * out is a vector of nm spinor fields
 */
void quark_propagator_QMR(FILE *propfile, unsigned int ssite, int nm, float *mass) {
  mshift_par QMR_par;
  int i;
	int source;
  double *shift;
  suNf_spinor *in=0;
  suNf_spinor_dble **resdn=0;
  suNf_spinor_dble **resd=0;
	int cgiter=0;
	double norm;
  suNf_spinor *test=0;

  /* allocate input spinor field */
  in = (suNf_spinor *)malloc(sizeof(suNf_spinor)*VOLUME);
	resdn=(suNf_spinor_dble**)malloc(sizeof(suNf_spinor_dble*)*nm);
	resdn[0]=(suNf_spinor_dble*)malloc(sizeof(suNf_spinor_dble)*nm*VOLUME);
	for(i=1;i<nm;++i)
		resdn[i]=resdn[i-1]+VOLUME;
	resd=(suNf_spinor_dble**)malloc(sizeof(suNf_spinor_dble*)*nm);
	resd[0]=(suNf_spinor_dble*)malloc(sizeof(suNf_spinor_dble)*nm*VOLUME);
	for(i=1;i<nm;++i)
		resd[i]=resd[i-1]+VOLUME;
  test = (suNf_spinor *)malloc(sizeof(suNf_spinor)*VOLUME);

  set_spinor_len(VOLUME);

	/* set up inverters parameters */
  shift=(double*)malloc(sizeof(double)*(nm));
  hmass=mass[0]; /* we can put any number here!!! */
  for(i=0;i<nm;++i){
    shift[i]=hmass-mass[i];
  }
  QMR_par.n = nm;
  QMR_par.shift = shift;
  QMR_par.err2 = .5e-10;
  QMR_par.max_iter = 0;

  /* noisy background */
	gaussian_spinor_field(in);
	for (source=0;source<NF*4*2;++source) {
		*(((float *) in)+(NF*8*ssite+source))=0.; /* zero in source */
	}
	norm=sqrt(spinor_field_sqnorm_f(in));
	spinor_field_mul_f(in,1./norm,in);
	
	/* invert noise */
  cgiter+=g5QMR_mshift(&QMR_par, &D, in, resdn);

	/* now loop over sources */
	for (source=0;source<4*NF;++source){
		/* put in source */
		*(((float *) in)+(NF*8*ssite+2*source))=1.;
		cgiter+=g5QMR_mshift(&QMR_par, &D, in, resd);

		for(i=0;i<QMR_par.n;++i){
			spinor_field_sub_dble_f(resd[i],resd[i],resdn[i]); /* compute difference */
			assign_sd2s(VOLUME,(suNf_spinor*)resd[i],resd[i]);

			/* this is a test of the solution */
			D(test,(suNf_spinor *)resd[i]);
			spinor_field_mul_add_assign_f(test,-QMR_par.shift[i],(suNf_spinor*)resd[i]);
			*(((float *) test)+(NF*8*ssite+2*source))-=1.;
			norm=spinor_field_sqnorm_f(test);
			if(norm>QMR_par.err2*2.)
				printf("test of difference g5QMR[%d] = %e\n",i,norm);

			/* write propagator on file */
			/* multiply by g_5 to macth the MINRES version */
			spinor_field_g5_f((suNf_spinor*)resd[i],(suNf_spinor*)resd[i]);
			error(fwrite(resd[i],(size_t) sizeof(suNf_spinor),(size_t)(VOLUME),propfile)!=(VOLUME),1,"Main",
					"Failed to write quark propagator to file");
		}

		/* remove source */
		*(((float *) in)+(NF*8*ssite+2*source))=0.;
	}

	printf("Quark Propagator QMR: %d\n",cgiter);
  
  /* free memory */
  free(in);
	free(resd[0]);
	free(resd);
	free(resdn[0]);
	free(resdn);
  free(shift);
	free(test);

}


/*
 * Computes the matrix elements (H^-1)_{x,0}
 * H is the hermitean dirac operator
 * out is a vector of nm spinor fields
 */
void quark_propagator(unsigned int source, int nm, float *mass, suNf_spinor **out) {
  static MINRES_par MINRESpar;
  int i,cgiter;
  suNf_spinor *in;

  /* allocate input spinor field */
  in = (suNf_spinor *)malloc(sizeof(suNf_spinor)*VOLUME);

  set_spinor_len(VOLUME);

  /* the source is on the first even site */
  spinor_field_zero_f(in);
  *(((float *) in)+2*source)=1.; /* put in source */

  hmass=mass[0];

  MINRESpar.err2 = 1.e-10;
  MINRESpar.max_iter = 0;

	cgiter=0;
  cgiter+=MINRES(&MINRESpar, &H, in, out[0],0);
  for(i=1;i<nm;++i){
    hmass=mass[i];
    cgiter+=MINRES(&MINRESpar, &H, in, out[i],out[i-1]);
  }
	printf("Quark Propagator MINRES: %d\n",cgiter);

  /* free input spinor field */
  free(in);

}
