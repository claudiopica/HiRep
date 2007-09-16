/******************************************************************************
*
* File dublin_mesons.c
*
* This program computes the generic two-point mesonic correlation function.
* It uses the algorithm of the group of Doublin (see Foley, Juge, O'Cais,
* Peardon, Ryan, Skullerud, hep-lat/0505023), for the all-to-all propagators.
* The inverse of g5D is splitted in the contribution of the low-lying
* eigenvalues and the rest. The former part is exactly computed. The latter
* one is stochastically estimated.
* The technique of time-dilution is implemeted, in order to reduce the
* stochastic fluctuations.
*
*
* Parameters.
* n_eigenvalues : the number of eigenvalues needed to compute the firts part
*                 of the quark propagator
* n_hp_eigenvalues : the number of eigenvalues (it must be less or equal than
*                 n_eigenvalues) to be computed with the given precision
*                 (see omega1, omega2)
* n_global_noisy_sources_per_point : the number of stochastic sources needed
*                 to estimate the second part of the quark propagator
* omega1 : the absolute precision required to compute 'n_hp_eigenvalues'
*                 eigenvalues (see eva.c)
* omega2 : the relative precision required to compute 'n_hp_eigenvalues'
*                 eigenvalues (see eva.c)
* acc : the precision to be used in the inversion routines
*
*
* Only the function
*       void dublin_meson_correlators(double** correlator[],
*            char corr_name[][256], int n_corr, int nm, double *mass)
* is accessible from the extern.
* It computes the correlators defined in the program, for all the masses
* mass[0]...mass[nm]. The correlators are stored in the array
* correlator[0...n_corr-1][0...nm-1][0...T-1]. The first index identify
* which correlator (pi, rho, eta ...). The string array corr_name[0...n_corr]
* contains the name of the computed correlators.
*
*
* Author: Agostino Patella
*
******************************************************************************/



#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "update.h"
#include "error.h"
#include "logger.h"
#include "memory.h"
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <string.h>



#define TIME_DILUTION
/* #define TESTINGMODE */
/* #define NO_NOISERECYCLING */


#ifndef TIME_DILUTION
#define NO_DILUTION
#endif

#define SOUCE_SINK_INDEX(p,q) ( (p)*n_sources + (q) )
#define SPIN_2D_INDEX(i,j) ( (i)*4 + (j) )
#define SOURCE_INDEX(m,r,t) ( (m)*n_diluted_noisy_sources + (r)*n_dilution_slices + (t) )
#define SINK_INDEX(m,r,t) ( ((r)*n_dilution_slices + (t))*n_masses + (m) )
#define DILUTION_INDEX(r,t) ( (r)*n_dilution_slices + (t) )
#define EV_INDEX(m,a) ( (m)*n_eigenvalues + (a) )


static int loglevel = 0;

static int n_masses;
static double hmass;
static double *shift;
static double *mass;

static const float omega1 = 1e-6;
static const float omega2 = 1e-3;
static const float acc = 1e-9;

enum{ n_eigenvalues = 12 }; /* N_{ev} */
enum{ n_hp_eigenvalues = 12 };

enum{ n_global_noisy_sources_per_point = 2 }; /* N_r */
enum{ n_points = 2 };
enum{ n_global_noisy_sources = n_global_noisy_sources_per_point * n_points }; /* N_{gns} */
#ifdef TIME_DILUTION
enum{ n_dilution_slices = T }; /* N_d */
#endif
#ifdef NO_DILUTION
enum{ n_dilution_slices = 1 }; /* N_d */
#endif
enum{ n_diluted_noisy_sources = n_global_noisy_sources*n_dilution_slices };

enum{ n_sources = n_eigenvalues + n_diluted_noisy_sources };

enum{ n_gamma_matrices = 16 };
enum{ n_correlators = 17 };

static double *d;          /* [i*n_eigenvalues+j], i<n_masses, j<n_eigenvalues */
static suNf_spinor **ev;   /* [i*n_eigenvalues+j], i<n_masses, j<n_eigenvalues */

static suNf_spinor **noisy_sources; /* [i*n_diluted_noisy_sources+j], i<n_masses, j<n_diluted_noisy_sources */
static suNf_spinor **noisy_sinks;   /* [i*n_masses+j], i<n_masses, j<n_diluted_noisy_sources */

static complex meson[n_gamma_matrices][n_sources*n_sources][T];



static void H(suNf_spinor *out, suNf_spinor *in);
static void D(suNf_spinor *out, suNf_spinor *in);

static void all_to_all_quark_propagator_init();

static void z2_spinor(suNf_spinor *source);

static void get_sources(suNf_spinor **source);
static void get_time_diluted_sources(suNf_spinor **source);
static void get_sinks_QMR(suNf_spinor *source, suNf_spinor **sink);
static void get_sinks_MINRES(suNf_spinor *source, suNf_spinor **sink);
static void project_sources(suNf_spinor **source, suNf_spinor **vectors);
static void project_time_diluted_sources(suNf_spinor **source, suNf_spinor **vectors);
static void project_sinks(suNf_spinor *sink, suNf_spinor **vectors);

static void source_sink_contraction(complex out[][16], suNf_spinor *source, suNf_spinor *sink, double z);

static void triplet_correlator_re(double* out, complex A[][T], complex B[][T]);
static void triplet_correlator_im(double* out, complex A[][T], complex B[][T]);
static void hairpin_re(double* out, complex A[][T], complex B[][T]);
static void hairpin_im(double* out, complex A[][T], complex B[][T]);


static void id_trace_H(complex* out, complex* smat); static int ID = -1;
static void g0_trace_H(complex* out, complex* smat); static int G0 = -1;
static void g5_trace_H(complex* out, complex* smat); static int G5 = -1;
static void g0g5_trace_H(complex* out, complex* smat); static int G0G5 = -1;
static void g1_trace_H(complex* out, complex* smat); static int G1 = -1;
static void g2_trace_H(complex* out, complex* smat); static int G2 = -1;
static void g3_trace_H(complex* out, complex* smat); static int G3 = -1;
static void g0g1_trace_H(complex* out, complex* smat); static int G0G1 = -1;
static void g0g2_trace_H(complex* out, complex* smat); static int G0G2 = -1;
static void g0g3_trace_H(complex* out, complex* smat); static int G0G3 = -1;
static void g5g1_trace_H(complex* out, complex* smat); static int G5G1 = -1;
static void g5g2_trace_H(complex* out, complex* smat); static int G5G2 = -1;
static void g5g3_trace_H(complex* out, complex* smat); static int G5G3 = -1;
static void g0g5g1_trace_H(complex* out, complex* smat); static int G0G5G1 = -1;
static void g0g5g2_trace_H(complex* out, complex* smat); static int G0G5G2 = -1;
static void g0g5g3_trace_H(complex* out, complex* smat); static int G0G5G3 = -1;


#ifdef TESTINGMODE
static complex trace_ave, trace_err, mesonID_ave[T],  mesonID_err[T], mesonG5_ave[T], mesonG5_err[T];
static int static_counter = 0;
#endif


void dublin_meson_correlators(double** correlator[], char corr_name[][256], int n_corr, int nm, double *mptr) {
	int m, r, p, q, t, counter;
#ifdef TESTINGMODE
	int x, i;
	double oldhmass, norm, ave;
	complex ctmp, trace, mesonID[T], mesonG5[T];
#endif
	int status;
	
	suNf_spinor *source;
	suNf_spinor *sink;
	
	complex ss[T][16];
	double tmpcorr[T];

	suNf_spinor *test=0;


	/* static parameters */

	n_masses = nm;
	shift=(double*)malloc(sizeof(double)*n_masses);
	mass = mptr;
	hmass = mass[0]; /* we can put any number here!!! */
	for(m = 0; m < n_masses; ++m){
		shift[m] = hmass - mass[m];
	}

	lprintf("DUBLIN_MESON_CORRELATORS",loglevel+1,"Number of masses = %d\n",n_masses);
	lprintf("DUBLIN_MESON_CORRELATORS",loglevel+1,"Highest mass = %f\n",hmass);
	for(m = 0;m < n_masses; m++){
		lprintf("DUBLIN_MESON_CORRELATORS",loglevel+1,"Shifts: %f\n",shift[m]);
	}

	error(n_corr < n_correlators,1,"meson_correlators [dublin_mesons.c]","Bad dimension for correlator[][]");


	/* allocate memory */

	set_spinor_len(VOLUME);
	test = alloc_spinor_field_f();
	all_to_all_quark_propagator_init();


	/* compute the lowest n_eigenvalues eigenvalues/vectors */

	if(n_eigenvalues > 0) {
		dirac_eva(n_hp_eigenvalues,n_eigenvalues,100,10*n_eigenvalues,omega1,omega2,n_masses,mass,ev,d,&status);
	}

	
	/* generate random sources & sinks */

	for(r = 0; r < n_global_noisy_sources; r++) {

#ifdef NO_DILUTION
		get_sources(noisy_sources + SOURCE_INDEX(0,r,0));
		for(m = 1; m < n_masses; m++)
			memcpy(noisy_sources[SOURCE_INDEX(m,r,0)],noisy_sources[SOURCE_INDEX(0,r,0)],sizeof(suNf_spinor)*VOLUME);
#endif
#ifdef TIME_DILUTION
		get_time_diluted_sources(noisy_sources + SOURCE_INDEX(0,r,0));
		for(m = 1; m < n_masses; m++)
		for(t = 0; t < n_dilution_slices; t++)
			memcpy(noisy_sources[SOURCE_INDEX(m,r,t)],noisy_sources[SOURCE_INDEX(0,r,t)],sizeof(suNf_spinor)*VOLUME);
#endif

#ifdef TESTINGMODE
		for(m = 0;m < n_masses; m++)
		for(t = 0; t < n_dilution_slices; t++) {
			norm = spinor_field_sqnorm_f(noisy_sources[SOURCE_INDEX(m,r,t)]);
			ave = 0.;
			for(x = 0; x < sizeof(suNf_spinor)*VOLUME/sizeof(double); x++)
				ave += ((double*)(noisy_sources[SOURCE_INDEX(m,r,t)]))[x];
			lprintf("DUBLIN_MESON_CORRELATORS",loglevel,"Testing noisy sources [m=%d,r=%d,t=%d]. norm/(4*Nf*VOL3)=%f ave=%e\n",m,r,t,norm/(4*VOL3*NF),.5*ave/norm);
		}
#endif
		
#ifdef NO_DILUTION
		get_sinks_MINRES(noisy_sources[SOURCE_INDEX(0,r,0)], noisy_sinks+SINK_INDEX(0,r,0));
		
		if(n_eigenvalues > 0) {
			for(m = 0; m < n_masses; m++) {
				project_sources(noisy_sources + SOURCE_INDEX(m,r,0), ev + EV_INDEX(m,0));
				project_sinks(noisy_sinks[SINK_INDEX(m,r,0)], ev + EV_INDEX(m,0));
			}
		}
#endif
#ifdef TIME_DILUTION
		for(t = 0; t < n_dilution_slices; t++)
			get_sinks_MINRES(noisy_sources[SOURCE_INDEX(0,r,t)], noisy_sinks+SINK_INDEX(0,r,t));
		
		if(n_eigenvalues > 0) {
			for(m = 0; m < n_masses; m++) {
				project_time_diluted_sources(noisy_sources + SOURCE_INDEX(m,r,0), ev + EV_INDEX(m,0));
				for(t = 0; t < n_dilution_slices; t++)
					project_sinks(noisy_sinks[SINK_INDEX(m,r,t)], ev + EV_INDEX(m,0));
			}
		}
#endif


#ifdef TESTINGMODE
		for(m = 0; m < n_masses; m++)
		for(t = 0; t < n_dilution_slices; t++) {
			norm = 0.;
			for(p = 0; p < n_eigenvalues; p++) {
				ctmp.re = spinor_field_prod_re_f(noisy_sources[SOURCE_INDEX(m,r,t)],ev[EV_INDEX(m,p)]);
				ctmp.im = spinor_field_prod_im_f(noisy_sources[SOURCE_INDEX(m,r,t)],ev[EV_INDEX(m,p)]);
				norm += ctmp.re*ctmp.re+ctmp.im*ctmp.im;
			}
			norm = sqrt(norm);
			lprintf("DUBLIN_MESON_CORRELATORS",loglevel,"Testing sources projection [m=%d,r=%d,t=%d] = %e\n",m,r,t,norm);
		}

		for(m = 0; m < n_masses; m++)
		for(t = 0; t < n_dilution_slices; t++) {
			norm = 0.;
			for(p = 0; p < n_eigenvalues; p++) {
				ctmp.re = spinor_field_prod_re_f(noisy_sinks[SINK_INDEX(m,r,t)],ev[EV_INDEX(m,p)]);
				ctmp.im = spinor_field_prod_im_f(noisy_sinks[SINK_INDEX(m,r,t)],ev[EV_INDEX(m,p)]);
				norm += ctmp.re*ctmp.re+ctmp.im*ctmp.im;
			}
			norm = sqrt(norm);
			lprintf("DUBLIN_MESON_CORRELATORS",loglevel,"Testing sinks projection [m=%d,r=%d,t=%d] = %e\n",m,r,t,norm);
		}

		oldhmass = hmass;
		for(m = 0; m < n_masses; m++)
		for(t = 0; t < n_dilution_slices; t++) {
			hmass = mass[m];
			H(test,noisy_sinks[SINK_INDEX(m,r,t)]);
			spinor_field_mul_add_assign_f(test,-1.,noisy_sources[SOURCE_INDEX(m,r,t)]);
			norm=spinor_field_sqnorm_f(test);
			/*if(norm>acc)*/
			lprintf("DUBLIN_MESON_CORRELATORS",loglevel,"Testing inversion [m=%d,r=%d,t=%d] = %e\n",m,r,t,norm);
		}
		hmass = oldhmass;
#endif


	}
	lprintf("DUBLIN_MESON_CORRELATORS",loglevel,"Generated %d noisy sources and relative sinks.\n",n_diluted_noisy_sources);


	/* two-point functions */

	for(m = 0; m < n_masses; m++) {

#ifdef TESTINGMODE
		for(t = 0; t < T; t++) {
			mesonID[t].re = 0.;
			mesonID[t].im = 0.;
			mesonG5[t].re = 0.;
			mesonG5[t].im = 0.;
		}
#endif
		
		for(p = 0; p < n_sources; p++) {
			
			if(p < n_eigenvalues) source = ev[EV_INDEX(m,p)];
			else source = noisy_sources[SOURCE_INDEX(m,0,p-n_eigenvalues)];
			
			for(q = 0; q < n_sources; q++) {
				
				if(q < n_eigenvalues) {
					sink = ev[EV_INDEX(m,q)];
					source_sink_contraction(ss, source, sink, 1.0f/d[EV_INDEX(m,q)]);
				} else {
					sink = noisy_sinks[SINK_INDEX(m,0,q-n_eigenvalues)];
					source_sink_contraction(ss, source, sink, 2.0f/n_global_noisy_sources);
				}
#ifdef TESTINGMODE
				ctmp.re = ctmp.im = 0.;
				for(t = 0; t < T; t++)
				for(i = 0; i < 4; i++) {
					ctmp.re += ss[t][SPIN_2D_INDEX(i,i)].re*VOL3;
					ctmp.im += ss[t][SPIN_2D_INDEX(i,i)].im*VOL3;
				}
				if(q < n_eigenvalues) {
					ctmp.re -= spinor_field_prod_re_f(source,sink)* 1.0f/d[EV_INDEX(m,q)];
					ctmp.im -= spinor_field_prod_im_f(source,sink)* 1.0f/d[EV_INDEX(m,q)];
				} else {
					ctmp.re -= spinor_field_prod_re_f(source,sink)* 2.0f/n_global_noisy_sources;
					ctmp.im -= spinor_field_prod_im_f(source,sink)* 2.0f/n_global_noisy_sources;
				}
				norm = sqrt(ctmp.re*ctmp.re+ctmp.im*ctmp.im);
				lprintf("DUBLIN_MESON_CORRELATORS",loglevel,"Testing contraction [m=%d,p=%d,q=%d] = %e\n",m,p,q,norm);
#endif
				
				for(t = 0; t < T; t++) {
					counter = 0;
					
					id_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					ID = counter;
					counter++;

					g0_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G0 = counter;
					counter++;

					g5_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G5 = counter;
					counter++;

					g0g5_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G0G5 = counter;
					counter++;

					g1_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G1 = counter;
					counter++;

					g2_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G2 = counter;
					counter++;

					g3_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G3 = counter;
					counter++;

					g0g1_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G0G1 = counter;
					counter++;

					g0g2_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G0G2 = counter;
					counter++;

					g0g3_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G0G3 = counter;
					counter++;

					g5g1_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G5G1 = counter;
					counter++;

					g5g2_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G5G2 = counter;
					counter++;

					g5g3_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G5G3 = counter;
					counter++;

					g0g5g1_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G0G5G1 = counter;
					counter++;

					g0g5g2_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G0G5G2 = counter;
					counter++;

					g0g5g3_trace_H(meson[counter][SOUCE_SINK_INDEX(p,q)]+t, ss[t]);
					G0G5G3 = counter;
					counter++;
					
				}
#ifdef TESTINGMODE
				ctmp.re = ctmp.im = 0.;
				for(t = 0; t < T; t++) {
					ctmp.re += meson[G5][SOUCE_SINK_INDEX(p,q)][t].re*VOL3;
					ctmp.im += meson[G5][SOUCE_SINK_INDEX(p,q)][t].im*VOL3;
				}
				if(q < n_eigenvalues) {
					ctmp.re -= spinor_field_prod_re_f(source,sink)* 1.0f/d[EV_INDEX(m,q)];
					ctmp.im -= spinor_field_prod_im_f(source,sink)* 1.0f/d[EV_INDEX(m,q)];
				} else {
					ctmp.re -= spinor_field_prod_re_f(source,sink)* 2.0f/n_global_noisy_sources;
					ctmp.im -= spinor_field_prod_im_f(source,sink)* 2.0f/n_global_noisy_sources;
				}
				norm = sqrt(ctmp.re*ctmp.re+ctmp.im*ctmp.im);
				lprintf("DUBLIN_MESON_CORRELATORS",loglevel,"Testing meson[G5] [m=%d,p=%d,q=%d] = %e\n",m,p,q,norm);

				ctmp.re = ctmp.im = 0.;
				for(t = 0; t < T; t++) {
					ctmp.re += meson[ID][SOUCE_SINK_INDEX(p,q)][t].re*VOL3;
/* 					ctmp.im += meson[ID][SOUCE_SINK_INDEX(p,q)][t].im*VOL3; */
				}
				if(q < n_eigenvalues) {
					ctmp.re -= spinor_field_g5_prod_re_f(source,sink)* 1.0f/d[EV_INDEX(m,q)];
/* 					ctmp.im -= spinor_field_g5_prod_im_f(source,sink)* 1.0f/d[EV_INDEX(m,q)]; */
				} else {
					ctmp.re -= spinor_field_g5_prod_re_f(source,sink)* 2.0f/n_global_noisy_sources;
/* 					ctmp.im -= spinor_field_g5_prod_im_f(source,sink)* 2.0f/n_global_noisy_sources; */
				}
/* 				norm = sqrt(ctmp.re*ctmp.re+ctmp.im*ctmp.im); */
				norm = sqrt(ctmp.re*ctmp.re);
				lprintf("DUBLIN_MESON_CORRELATORS",loglevel,"Testing meson[ID] [m=%d,p=%d,q=%d] = %e\n",m,p,q,norm);


				if(q < n_eigenvalues) {
					for(t = 0; t < T; t++) {
						mesonID[t].re += VOL3* meson[ID][SOUCE_SINK_INDEX(p,q)][t].re;
						mesonID[t].im += VOL3* meson[ID][SOUCE_SINK_INDEX(p,q)][t].im;
						mesonG5[t].re += VOL3* meson[G5][SOUCE_SINK_INDEX(p,q)][t].re;
						mesonG5[t].im += VOL3* meson[G5][SOUCE_SINK_INDEX(p,q)][t].im;
					}
				} else {
					for(t = 0; t < T; t++) {
						mesonID[t].re += VOL3* .5* meson[ID][SOUCE_SINK_INDEX(p,q)][t].re;
						mesonID[t].im += VOL3* .5* meson[ID][SOUCE_SINK_INDEX(p,q)][t].im;
						mesonG5[t].re += VOL3* .5* meson[G5][SOUCE_SINK_INDEX(p,q)][t].re;
						mesonG5[t].im += VOL3* .5* meson[G5][SOUCE_SINK_INDEX(p,q)][t].im;
					}
				}
#endif
				lprintf("DUBLIN_MESON_CORRELATORS",loglevel+2,"Written %d gamma matrices [m=%d,p=%d,q=%d].\n",counter,m,p,q);
				
			}
		}

#ifdef TESTINGMODE
		static_counter++;
		trace.re = trace.im = 0.;
		for(t = 0; t < T; t++) {
			trace.re += mesonID[t].re;
			trace.im += mesonID[t].im;
			mesonID_ave[t].re += mesonID[t].re;
			mesonID_ave[t].im += mesonID[t].im;
			mesonG5_ave[t].re += mesonG5[t].re;
			mesonG5_ave[t].im += mesonG5[t].im;
			mesonID_err[t].re += mesonID[t].re*mesonID[t].re;
			mesonID_err[t].im += mesonID[t].im*mesonID[t].im;
			mesonG5_err[t].re += mesonG5[t].re*mesonG5[t].re;
			mesonG5_err[t].im += mesonG5[t].im*mesonG5[t].im;
			lprintf("DUBLIN_MESON_CORRELATORS",loglevel," T2 mesonID [m=%d,t=%d] = %f(%f) %f(%f)\n",m,t,
				mesonID_ave[t].re/static_counter,
				sqrt(mesonID_err[t].re-mesonID_ave[t].re*mesonID_ave[t].re/(static_counter*static_counter))/static_counter,
				mesonID_ave[t].im/static_counter,
				sqrt(mesonID_err[t].im-mesonID_ave[t].im*mesonID_ave[t].im/(static_counter*static_counter))/static_counter);
			lprintf("DUBLIN_MESON_CORRELATORS",loglevel," T3 mesonG5 [m=%d,t=%d] = %f(%f) %f(%f)\n",m,t,
				mesonG5_ave[t].re/static_counter,
				sqrt(mesonG5_err[t].re-mesonG5_ave[t].re*mesonG5_ave[t].re/(static_counter*static_counter))/static_counter,
				mesonG5_ave[t].im/static_counter,
				sqrt(mesonG5_err[t].im-mesonG5_ave[t].im*mesonG5_ave[t].im/(static_counter*static_counter))/static_counter);
		}
		trace_ave.re += trace.re;
		trace_ave.im += trace.im;
		trace_err.re += trace.re*trace.re;
		trace_err.im += trace.im*trace.im;
		lprintf("DUBLIN_MESON_CORRELATORS",loglevel," T1 Trace(invD) [m=%d] = %f(%f) %f(%f)\n",m,
				trace_ave.re/static_counter,
				sqrt(trace_err.re-trace_ave.re*trace_ave.re/(static_counter*static_counter))/static_counter,
				trace_ave.im/static_counter,
				sqrt(trace_err.im-trace_ave.im*trace_ave.im/(static_counter*static_counter))/static_counter);
#endif

	
		counter = 0;
		
		
		/* scalar C=+1; isotriplet=a isosinglet=f */
		
		if(m == 0) strcpy(corr_name[counter],"a");
		triplet_correlator_re(correlator[counter][m], meson[ID], meson[ID]);
		counter++;

		if(m == 0) strcpy(corr_name[counter],"f");
		hairpin_re(correlator[counter][m], meson[ID], meson[ID]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] = correlator[counter-1][m][t] + correlator[counter][m][t];
		counter++;

		
		/* pseudoscalar C=+1; isotriplet=pi isosinglet=eta */
		
		if(m == 0) strcpy(corr_name[counter],"pi");
		triplet_correlator_re(correlator[counter][m], meson[G5], meson[G5]);
		counter++;

		if(m == 0) strcpy(corr_name[counter],"eta");
		hairpin_re(correlator[counter][m], meson[G5], meson[G5]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] = correlator[counter-1][m][t] + correlator[counter][m][t];
		counter++;
		
		
		/* vector C=-1; isotriplet=rho isosinglet=phi */
		
		if(m == 0) strcpy(corr_name[counter],"rho");
		triplet_correlator_re(correlator[counter][m], meson[G1], meson[G1]);
		triplet_correlator_re(tmpcorr, meson[G2], meson[G2]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] += tmpcorr[t];
		triplet_correlator_re(tmpcorr, meson[G3], meson[G3]);
		for(t = 0; t < T; t++) {
			correlator[counter][m][t] += tmpcorr[t];
			correlator[counter][m][t] /= 3.0;
		}
		counter++;
		
		if(m == 0) strcpy(corr_name[counter],"phi");
		hairpin_re(correlator[counter][m], meson[G1], meson[G1]);
		hairpin_re(tmpcorr, meson[G2], meson[G2]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] += tmpcorr[t];
		hairpin_re(tmpcorr, meson[G3], meson[G3]);
		for(t = 0; t < T; t++) {
			correlator[counter][m][t] += tmpcorr[t];
			correlator[counter][m][t] /= 3.0;
			correlator[counter][m][t] = correlator[counter-1][m][t] + correlator[counter][m][t];
		}		
		counter++;
		
		
		/* pseudovector C=-1; isotriplet=b isosinglet=h */
		
		if(m == 0) strcpy(corr_name[counter],"b");
		triplet_correlator_re(correlator[counter][m], meson[G0G5G1], meson[G0G5G1]);
		triplet_correlator_re(tmpcorr, meson[G0G5G2], meson[G0G5G2]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] += tmpcorr[t];
		triplet_correlator_re(tmpcorr, meson[G0G5G3], meson[G0G5G3]);
		for(t = 0; t < T; t++) {
			correlator[counter][m][t] += tmpcorr[t];
			correlator[counter][m][t] /= 3.0;
		}
		counter++;
		
		if(m == 0) strcpy(corr_name[counter],"h");
		hairpin_re(correlator[counter][m], meson[G0G5G1], meson[G0G5G1]);
		hairpin_re(tmpcorr, meson[G0G5G2], meson[G0G5G2]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] += tmpcorr[t];
		hairpin_re(tmpcorr, meson[G0G5G3], meson[G0G5G3]);
		for(t = 0; t < T; t++) {
			correlator[counter][m][t] += tmpcorr[t];
			correlator[counter][m][t] /= 3.0;
			correlator[counter][m][t] = correlator[counter-1][m][t] + correlator[counter][m][t];
		}		
		counter++;

		
		/* another pseudoscalar C=+1; isotriplet=pi isosinglet=eta */
		
		if(m == 0) strcpy(corr_name[counter],"pi2");
		triplet_correlator_re(correlator[counter][m], meson[G0G5], meson[G0G5]);
		counter++;

		if(m == 0) strcpy(corr_name[counter],"eta2");
		hairpin_re(correlator[counter][m], meson[G0G5], meson[G0G5]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] = correlator[counter-1][m][t] + correlator[counter][m][t];
		counter++;
		
		
		/* another vector C=-1; isotriplet=rho isosinglet=phi */
		
		if(m == 0) strcpy(corr_name[counter],"rho2");
		triplet_correlator_re(correlator[counter][m], meson[G0G1], meson[G0G1]);
		triplet_correlator_re(tmpcorr, meson[G0G2], meson[G0G2]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] += tmpcorr[t];
		triplet_correlator_re(tmpcorr, meson[G0G3], meson[G0G3]);
		for(t = 0; t < T; t++) {
			correlator[counter][m][t] += tmpcorr[t];
			correlator[counter][m][t] /= 3.0;
		}
		counter++;
		
		if(m == 0) strcpy(corr_name[counter],"phi2");
		hairpin_re(correlator[counter][m], meson[G0G1], meson[G0G1]);
		hairpin_re(tmpcorr, meson[G0G2], meson[G0G2]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] += tmpcorr[t];
		hairpin_re(tmpcorr, meson[G0G3], meson[G0G3]);
		for(t = 0; t < T; t++) {
			correlator[counter][m][t] += tmpcorr[t];
			correlator[counter][m][t] /= 3.0;
			correlator[counter][m][t] = correlator[counter-1][m][t] + correlator[counter][m][t];
		}		
		counter++;
		
		
		/* scalar C=-1; forbidden Minkowski mesons */
		
		if(m == 0) strcpy(corr_name[counter],"forbidden triplet 0+-");
		triplet_correlator_re(correlator[counter][m], meson[G0], meson[G0]);
		counter++;

		if(m == 0) strcpy(corr_name[counter],"forbidden singlet 0+-");
		hairpin_re(correlator[counter][m], meson[G0], meson[G0]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] = correlator[counter-1][m][t] + correlator[counter][m][t];
		counter++;
		
		
		/* pseudovector C=+1; forbidden Minkowski mesons */
		
		if(m == 0) strcpy(corr_name[counter],"forbidden triplet 1++");
		triplet_correlator_re(correlator[counter][m], meson[G5G1], meson[G5G1]);
		triplet_correlator_re(tmpcorr, meson[G5G2], meson[G5G2]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] += tmpcorr[t];
		triplet_correlator_re(tmpcorr, meson[G5G3], meson[G5G3]);
		for(t = 0; t < T; t++) {
			correlator[counter][m][t] += tmpcorr[t];
			correlator[counter][m][t] /= 3.0;
		}
		counter++;
		
		if(m == 0) strcpy(corr_name[counter],"forbidden singlet 1++");
		hairpin_re(correlator[counter][m], meson[G5G1], meson[G5G1]);
		hairpin_re(tmpcorr, meson[G5G2], meson[G5G2]);
		for(t = 0; t < T; t++)
			correlator[counter][m][t] += tmpcorr[t];
		hairpin_re(tmpcorr, meson[G5G3], meson[G5G3]);
		for(t = 0; t < T; t++) {
			correlator[counter][m][t] += tmpcorr[t];
			correlator[counter][m][t] /= 3.0;
			correlator[counter][m][t] = correlator[counter-1][m][t] + correlator[counter][m][t];
		}		
		counter++;
		
		
		/* pion decay amplitude */
		
		if(m == 0) strcpy(corr_name[counter],"pion decay amplitude");
		triplet_correlator_im(correlator[counter][m], meson[G5], meson[G0G5]);
		counter++;
		
		
	}
	
	free_field(test);
	free(shift);
}






static void H(suNf_spinor *out, suNf_spinor *in){
	g5Dphi(hmass,out,in);
}



static void D(suNf_spinor *out, suNf_spinor *in){
	Dphi(hmass,out,in);
}



static void all_to_all_quark_propagator_init() {
	int m;
	static int init_flag = 0;
	double requiredmemory = 0.0;
#ifdef TESTINGMODE
	int t;
#endif
	
	if(init_flag != 0) return;

	set_spinor_len(VOLUME);

	/* static complex meson[n_gamma_matrices][n_sources*n_sources][T]; */
	requiredmemory += sizeof(complex)*n_gamma_matrices*n_sources*n_sources*T;
	
	d = (double*)malloc(sizeof(double)*n_masses*n_eigenvalues);
	ev = (suNf_spinor**)malloc(sizeof(suNf_spinor*)*n_masses*n_eigenvalues);
	for(m = 0; m < n_masses*n_eigenvalues; m++)
		ev[m] = alloc_spinor_field_f();
	requiredmemory += sizeof(suNf_spinor)*VOLUME*n_eigenvalues*n_masses;
	
	noisy_sources = (suNf_spinor**)malloc(sizeof(suNf_spinor*)*n_masses*n_diluted_noisy_sources);
	noisy_sinks = (suNf_spinor**)malloc(sizeof(suNf_spinor*)*n_masses*n_diluted_noisy_sources);
	for(m = 0; m < n_masses*n_diluted_noisy_sources; m++) {
		noisy_sources[m] = alloc_spinor_field_f();
		noisy_sinks[m] = alloc_spinor_field_f();
	}
	requiredmemory += 2*sizeof(suNf_spinor)*VOLUME*n_masses*n_diluted_noisy_sources;
	
	requiredmemory /= 1024.0*1024.0;
	lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel,"Required memory = %f Mb\n",requiredmemory);

#ifdef TESTINGMODE
	trace_ave.re = trace_ave.im = 0.;
	trace_err.re = trace_err.im = 0.;
	for(t = 0; t < T; t++) {
		mesonID_ave[t].re = mesonID_ave[t].im = 0.;
		mesonID_err[T].re = mesonID_err[T].im = 0.;
		mesonG5_ave[T].re = mesonG5_ave[T].im = 0.;
		mesonG5_err[T].re = mesonG5_err[T].im = 0.;
	}
#endif
	
	init_flag = 1;
}



/*
suNf_spinor* get_stored_source(int random_index, int dilution_index) {
	return noisy_sources[random_index*T+dilution_index];
}



suNf_spinor** get_stored_sinks(int random_index, int dilution_index) {
	return noisy_sinks[random_index*T+dilution_index];
}
*/



static void z2_spinor(suNf_spinor *source) {	
	double plus = sqrt(.5);
	double minus = -sqrt(.5);
	unsigned int i;
	double *sf = (double*)source;
	for(i = 0; i < sizeof(suNf_spinor)/sizeof(double); i++) {
		sf[i] = ((rand() & 1) == 0) ? plus : minus;
	}
}



static void get_sources(suNf_spinor **source) {
	int x;
	
	for(x = 0; x < VOLUME; x++)
		z2_spinor(source[0]+x);
}



static void get_time_diluted_sources(suNf_spinor **source) {
	int t, x, index;
	
	for(t = 0; t < T; t++) {
		spinor_field_zero_f(source[t]);
		for(x = 0; x < VOL3; x++) {
			index = ipt_4d[t][x];
			z2_spinor(source[t] + index);
		}
	}
}



static void get_sinks_QMR(suNf_spinor *source, suNf_spinor **sink) {
	mshift_par QMR_par;
	int m;
	int cgiter=0;
	double norm;
	suNf_spinor *test=0;

	/* allocate test spinor field */
	test = alloc_spinor_field_f();

	/* set_spinor_len(VOLUME); */

	/* set up inverters parameters */
	QMR_par.n = n_masses;
	QMR_par.shift = shift;
	QMR_par.err2 = .5*acc;
	QMR_par.max_iter = 0;

	lprintf("GET_SINKS_QMR",loglevel+3,"Number of masses = %d\n",n_masses);
	lprintf("GET_SINKS_QMR",loglevel+3,"Highest mass = %f\n",hmass);
	for(m = 0; m < n_masses; m++){
		lprintf("GET_SINKS_QMR",loglevel+3,"Shifts: %f\n",shift[m]);
	}

	norm = sqrt(spinor_field_sqnorm_f(source));
	lprintf("GET_SINKS_QMR",loglevel+3,"Source norm = %f\n",norm);

	for(m=0;m<QMR_par.n;++m){
		spinor_field_zero_f(sink[m]);
	}

	D(test,source);
	++cgiter;
	norm=sqrt(spinor_field_sqnorm_f(test));
	lprintf("GET_SINKS_QMR",loglevel+3,"Test norm = %f\n",norm);
	for(m=0;m<QMR_par.n;++m){
		D(test,sink[m]);
		++cgiter;
		norm=sqrt(spinor_field_sqnorm_f(test));
		lprintf("GET_SINKS_QMR",loglevel+3,"Test norm = %f\n",norm);
	}

	cgiter+=g5QMR_mshift(&QMR_par, &D, source, sink);

	for(m=0;m<QMR_par.n;++m){
		/* this is a test of the solution */
		D(test,sink[m]);
		++cgiter;
		spinor_field_mul_add_assign_f(test,-QMR_par.shift[m],sink[m]);
		spinor_field_sub_f(test,test,sink[m]);
		norm=spinor_field_sqnorm_f(test);
		/*if(norm>acc)*/
		lprintf("GET_SINKS_QMR",loglevel+1,"g5QMR residuum of source [%d] = %e\n",m,norm);

		/* multiply by g_5 to match the MINRES version */
		spinor_field_g5_f(sink[m],sink[m]);
		
		/* convert to single precision */
		/* assign_sd2s(VOLUME,(suNf_spinor_flt*)resd[m],resd[m]);
		error(
				fwrite(resd[m],(size_t) sizeof(suNf_spinor_flt),(size_t)(VOLUME),propfile)!=(VOLUME),
				1,"quark_propagator_QMR","Failed to write quark propagator to file"); */
	}
	
	free_field(test);
}



static void get_sinks_MINRES(suNf_spinor *source, suNf_spinor **sink) {
	static MINRES_par MINRESpar;
	int m;
	int cgiter=0;
	double oldhmass;

	/* set_spinor_len(VOLUME); */

	/* set up inverters parameters */
	MINRESpar.err2 = acc;
	MINRESpar.max_iter = 0;

/* 	for(m=0;m<n_masses;++m){ */
/* 		spinor_field_zero_f(sink[m]); */
/* 	} */
	
	oldhmass = hmass;
	hmass = mass[0];
	cgiter = MINRES(&MINRESpar, &H, source, sink[0], 0);
	for(m = 1; m < n_masses; m++) {
		hmass = mass[m];
		cgiter += MINRES(&MINRESpar, &H, source, sink[m], sink[m-1]);
	}
	
	lprintf("GET_SINKS_MINRES",loglevel+1,"MINRES MVM = %d\n",cgiter);

	hmass = oldhmass;
}



static void project_sources(suNf_spinor **source, suNf_spinor **vectors) {
	int a, x;
	static complex alpha[n_eigenvalues];
	
	for(a = 0; a < n_eigenvalues; a++)
		alpha[a].re = alpha[a].im = 0.0;

	for(x = 0; x < VOLUME; x++) {
		for(a = 0; a < n_eigenvalues; a++) {
			_spinor_prod_assign_f(alpha[a],vectors[a][x],source[0][x]);
		}
	}

	for(a = 0; a < n_eigenvalues; a++) {
		for(x = 0; x < VOLUME; x++) {
			_spinor_project_f(source[0][x],alpha[a],vectors[a][x]);
		}
	}
}



static void project_time_diluted_sources(suNf_spinor **source, suNf_spinor **vectors) {
	int a, t, x, index;
	static complex alpha[n_eigenvalues];
	
	for(t = 0; t < T; t++) {
		for(a = 0; a < n_eigenvalues; a++)
			alpha[a].re = alpha[a].im = 0.0;

		for(x = 0; x < VOL3; x++) {
			index = ipt_4d[t][x];
			for(a = 0; a < n_eigenvalues; a++) {
				_spinor_prod_assign_f(alpha[a],vectors[a][index],source[t][index]);
			}
		}

		for(a = 0; a < n_eigenvalues; a++) {
			for(index = 0; index < VOLUME; index++) {
				_spinor_project_f(source[t][index],alpha[a],vectors[a][index]);
			}
		}
	}
}



static void project_sinks(suNf_spinor *sink, suNf_spinor **vectors) {
	int a, x;
	static complex alpha[n_eigenvalues];
	
	for(a = 0; a < n_eigenvalues; a++)
		alpha[a].re = alpha[a].im = 0.0;

	for(x = 0; x < VOLUME; x++)
	for(a = 0; a < n_eigenvalues; a++) {
		_spinor_prod_assign_f(alpha[a],vectors[a][x],sink[x]);
	}

	for(a = 0; a < n_eigenvalues; a++)
	for(x = 0; x < VOLUME; x++) {
		_spinor_project_f(sink[x],alpha[a],vectors[a][x]);
	}
}



static void source_sink_contraction(complex out[][16], suNf_spinor *source, suNf_spinor *sink, double z) {
	int i, j, t, x, index;
	suNf_vector *eta, *csi;
	
	for(t = 0; t < T; t++) {
		for(i = 0; i < 16; i++) {
/* 		memset(out[t], '\0', sizeof(complex)*16); */
			out[t][i].re = out[t][i].im = 0.;
		}
		
		for(x = 0; x < VOL3; x++) {
			index = ipt_4d[t][x];
			
			for(i = 0; i < 4; i++) {
				csi = (suNf_vector*)(sink+index) + i;
				for(j = 0; j < 4; j++) {
					eta = (suNf_vector*)(source+index) + j;
					_vector_prod_assign_f(out[t][SPIN_2D_INDEX(i,j)], *eta, *csi);
				}
			}
			
		}
		
		for(i = 0; i < 16; i++) {
			out[t][i].re *= z/VOL3;
			out[t][i].im *= z/VOL3;
		}
	}

	lprintf("SOURCE_SINK_CONTRACTION",loglevel+2,"Written in %p\n",out);

}



static void triplet_correlator_re(double* out, complex A[][T], complex B[][T]) {
	int j, k, rj, dj, rk, dk, t, dt, t1;
#ifndef NO_NOISERECYCLING
	const float z = (0.5f*n_global_noisy_sources)/(n_global_noisy_sources-1);
#endif
	complex *a;
	complex *b;
	
	for(dt = 0; dt < T; dt++)
		out[dt] = 0.0;
	
	for(j = 0; j < n_eigenvalues; j++) {
	
		for(k = 0; k < n_eigenvalues; k++) {
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += a->re*b->re - a->im*b->im;
				}
			}
		}
		
		for(rk = 0; rk < n_global_noisy_sources/2; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += a->re*b->re - a->im*b->im;
				}
			}
		}
		
	}
	
	for(rj = n_global_noisy_sources/2; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(k = 0; k < n_eigenvalues; k++) {
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += a->re*b->re - a->im*b->im;
				}
			}
		}
		
	}
	
#ifdef NO_NOISERECYCLING
	for(rj = n_global_noisy_sources/2; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(rk = 0; rk < n_global_noisy_sources/2; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += a->re*b->re - a->im*b->im;
				}
			}
		}
	
	}
	
#else
	for(rj = 1; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(rk = 0; rk < rj; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += z*(a->re*b->re - a->im*b->im);
				}
			}
		}
		
	}
#endif
	
	for(dt = 0; dt < T; dt++)
		out[dt] *= -1./T;
	
}



static void triplet_correlator_im(double* out, complex A[][T], complex B[][T]) {
	int j, k, rj, dj, rk, dk, t, dt, t1;
#ifndef NO_NOISERECYCLING
	const float z = (0.5f*n_global_noisy_sources)/(n_global_noisy_sources-1);
#endif
	complex *a;
	complex *b;
	
	for(dt = 0; dt < T; dt++)
		out[dt] = 0.0;
	
	for(j = 0; j < n_eigenvalues; j++) {
	
		for(k = 0; k < n_eigenvalues; k++) {
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += a->re*b->im + a->im*b->re;
				}
			}
		}
		
		for(rk = 0; rk < n_global_noisy_sources/2; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += a->re*b->im + a->im*b->re;
				}
			}
		}
		
	}
	
	for(rj = n_global_noisy_sources/2; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(k = 0; k < n_eigenvalues; k++) {
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += a->re*b->im + a->im*b->re;
				}
			}
		}
		
	}
	
	
#ifdef NO_NOISERECYCLING
	for(rj = n_global_noisy_sources/2; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(rk = 0; rk < n_global_noisy_sources/2; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += a->re*b->im + a->im*b->re;
				}
			}
		}
	
	}
	
#else
	for(rj = 1; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(rk = 0; rk < rj; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,k)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,j)]+t1;
					out[dt] += z*(a->re*b->im + a->im*b->re);
				}
			}
		}

	}
#endif
	
	for(dt = 0; dt < T; dt++)
		out[dt] *= -1./T;
	
}



static void hairpin_re(double* out, complex A[][T], complex B[][T]) {
	int j, k, rj, dj, rk, dk, t, dt, t1;
#ifndef NO_NOISERECYCLING
	const float z = (0.5f*n_global_noisy_sources)/(n_global_noisy_sources-1);
#endif
	complex *a;
	complex *b;
	
	for(dt = 0; dt < T; dt++)
		out[dt] = 0.0;
	
	for(j = 0; j < n_eigenvalues; j++) {
	
		for(k = 0; k < n_eigenvalues; k++) {
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += a->re*b->re - a->im*b->im;
				}
			}
		}
		
		for(rk = 0; rk < n_global_noisy_sources/2; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += a->re*b->re - a->im*b->im;
				}
			}
		}
		
	}
	
	for(rj = n_global_noisy_sources/2; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(k = 0; k < n_eigenvalues; k++) {
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += a->re*b->re - a->im*b->im;
				}
			}
		}
		
	}
	
	
#ifdef NO_NOISERECYCLING
	for(rj = n_global_noisy_sources/2; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(rk = 0; rk < n_global_noisy_sources/2; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += a->re*b->re - a->im*b->im;
				}
			}
		}
		
	}
	
#else
	for(rj = 1; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(rk = 0; rk < rj; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += z*(a->re*b->re - a->im*b->im);
				}
			}
		}
		
	}
#endif
	
	for(dt = 0; dt < T; dt++)
		out[dt] *= 1./T;
	
}



static void hairpin_im(double* out, complex A[][T], complex B[][T]) {
	int j, k, rj, dj, rk, dk, t, dt, t1;
#ifndef NO_NOISERECYCLING
	const float z = (0.5f*n_global_noisy_sources)/(n_global_noisy_sources-1);
#endif
	complex *a;
	complex *b;
	
	for(dt = 0; dt < T; dt++)
		out[dt] = 0.0;
	
	for(j = 0; j < n_eigenvalues; j++) {
	
		for(k = 0; k < n_eigenvalues; k++) {
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += a->re*b->im + a->im*b->re;
				}
			}
		}
		
		for(rk = 0; rk < n_global_noisy_sources/2; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += a->re*b->im + a->im*b->re;
				}
			}
		}
		
	}
	
	for(rj = n_global_noisy_sources/2; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(k = 0; k < n_eigenvalues; k++) {
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += a->re*b->im + a->im*b->re;
				}
			}
		}
		
	}
	
	
#ifdef NO_NOISERECYCLING
	for(rj = n_global_noisy_sources/2; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(rk = 0; rk < n_global_noisy_sources/2; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += a->re*b->im + a->im*b->re;
				}
			}
		}
		
	}
	
#else
	for(rj = 1; rj < n_global_noisy_sources; rj++)
	for(dj = 0; dj < n_dilution_slices; dj++) {
		j = n_eigenvalues + DILUTION_INDEX(rj,dj);
		
		for(rk = 0; rk < rj; rk++)
		for(dk = 0; dk < n_dilution_slices; dk++) {
			k = n_eigenvalues + DILUTION_INDEX(rk,dk);
			
			for(t = 0; t < T; t++) {
				a = A[SOUCE_SINK_INDEX(j,j)]+t;
				for(dt = 0; dt < T; dt++) {
					t1 = (t+dt)%T;
					b = B[SOUCE_SINK_INDEX(k,k)]+t1;
					out[dt] += z*(a->re*b->im + a->im*b->re);
				}
			}
		}
		
	}
#endif
	
	for(dt = 0; dt < T; dt++)
		out[dt] *= 1./T;
	
}



#define GAMMA_TRACE_RE_DEFINITION \
static void NAME(complex* out, complex* smat) { \
	out->re = _S1_*smat[SPIN_2D_INDEX(0,_C1_-1)].re \
		+ _S2_*smat[SPIN_2D_INDEX(1,_C2_-1)].re \
		+ _S3_*smat[SPIN_2D_INDEX(2,_C3_-1)].re \
		+ _S4_*smat[SPIN_2D_INDEX(3,_C4_-1)].re; \
	out->im = _S1_*smat[SPIN_2D_INDEX(0,_C1_-1)].im \
		+ _S2_*smat[SPIN_2D_INDEX(1,_C2_-1)].im \
		+ _S3_*smat[SPIN_2D_INDEX(2,_C3_-1)].im \
		+ _S4_*smat[SPIN_2D_INDEX(3,_C4_-1)].im; \
}



#define GAMMA_TRACE_IM_DEFINITION \
static void NAME(complex* out, complex* smat) { \
	out->im = _S1_*smat[SPIN_2D_INDEX(0,_C1_-1)].re \
		+ _S2_*smat[SPIN_2D_INDEX(1,_C2_-1)].re \
		+ _S3_*smat[SPIN_2D_INDEX(2,_C3_-1)].re \
		+ _S4_*smat[SPIN_2D_INDEX(3,_C4_-1)].re; \
	out->re = - _S1_*smat[SPIN_2D_INDEX(0,_C1_-1)].im \
		- _S2_*smat[SPIN_2D_INDEX(1,_C2_-1)].im \
		- _S3_*smat[SPIN_2D_INDEX(2,_C3_-1)].im \
		- _S4_*smat[SPIN_2D_INDEX(3,_C4_-1)].im; \
}



#define NAME id_trace_H

#define _C1_ 1
#define _C2_ 2
#define _C3_ 3
#define _C4_ 4

#define _S1_ 1
#define _S2_ 1
#define _S3_ -1
#define _S4_ -1

GAMMA_TRACE_RE_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g0_trace_H

#define _C1_ 3
#define _C2_ 4
#define _C3_ 1
#define _C4_ 2

#define _S1_ -1
#define _S2_ -1
#define _S3_ 1
#define _S4_ 1

GAMMA_TRACE_RE_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g5_trace_H

#define _C1_ 1
#define _C2_ 2
#define _C3_ 3
#define _C4_ 4

#define _S1_ 1
#define _S2_ 1
#define _S3_ 1
#define _S4_ 1

GAMMA_TRACE_RE_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g0g5_trace_H

#define _C1_ 3
#define _C2_ 4
#define _C3_ 1
#define _C4_ 2

#define _S1_ -1
#define _S2_ -1
#define _S3_ -1
#define _S4_ -1

GAMMA_TRACE_RE_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g1_trace_H

#define _C1_ 4
#define _C2_ 3
#define _C3_ 2
#define _C4_ 1

#define _S1_ 1
#define _S2_ 1
#define _S3_ 1
#define _S4_ 1

GAMMA_TRACE_IM_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g2_trace_H

#define _C1_ 4
#define _C2_ 3
#define _C3_ 2
#define _C4_ 1

#define _S1_ -1
#define _S2_ 1
#define _S3_ -1
#define _S4_ 1

GAMMA_TRACE_RE_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g3_trace_H

#define _C1_ 3
#define _C2_ 4
#define _C3_ 1
#define _C4_ 2

#define _S1_ -1
#define _S2_ 1
#define _S3_ -1
#define _S4_ 1

GAMMA_TRACE_IM_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g0g1_trace_H

#define _C1_ 2
#define _C2_ 1
#define _C3_ 4
#define _C4_ 3

#define _S1_ -1
#define _S2_ -1
#define _S3_ -1
#define _S4_ -1

GAMMA_TRACE_IM_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g0g2_trace_H

#define _C1_ 2
#define _C2_ 1
#define _C3_ 4
#define _C4_ 3

#define _S1_ 1
#define _S2_ -1
#define _S3_ 1
#define _S4_ -1

GAMMA_TRACE_RE_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g0g3_trace_H

#define _C1_ 1
#define _C2_ 2
#define _C3_ 3
#define _C4_ 4

#define _S1_ 1
#define _S2_ -1
#define _S3_ 1
#define _S4_ -1

GAMMA_TRACE_IM_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g5g1_trace_H

#define _C1_ 4
#define _C2_ 3
#define _C3_ 2
#define _C4_ 1

#define _S1_ -1
#define _S2_ -1
#define _S3_ 1
#define _S4_ 1

GAMMA_TRACE_IM_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g5g2_trace_H

#define _C1_ 4
#define _C2_ 3
#define _C3_ 2
#define _C4_ 1

#define _S1_ 1
#define _S2_ -1
#define _S3_ -1
#define _S4_ 1

GAMMA_TRACE_RE_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g5g3_trace_H

#define _C1_ 3
#define _C2_ 4
#define _C3_ 1
#define _C4_ 2

#define _S1_ 1
#define _S2_ -1
#define _S3_ -1
#define _S4_ 1

GAMMA_TRACE_IM_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g0g5g1_trace_H

#define _C1_ 2
#define _C2_ 1
#define _C3_ 4
#define _C4_ 3

#define _S1_ 1
#define _S2_ 1
#define _S3_ -1
#define _S4_ -1

GAMMA_TRACE_IM_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g0g5g2_trace_H

#define _C1_ 2
#define _C2_ 1
#define _C3_ 4
#define _C4_ 3

#define _S1_ -1
#define _S2_ 1
#define _S3_ 1
#define _S4_ -1

GAMMA_TRACE_RE_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_



#define NAME g0g5g3_trace_H

#define _C1_ 1
#define _C2_ 2
#define _C3_ 3
#define _C4_ 4

#define _S1_ -1
#define _S2_ 1
#define _S3_ 1
#define _S4_ -1

GAMMA_TRACE_IM_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_


