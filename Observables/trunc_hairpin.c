/******************************************************************************
*
* File hairpin_mesons.c
*
* This program computes the hairpin term of the isosinglet mesonic correlators.
* It uses the algorithm of the group of Dublin (see Foley, Juge, O'Cais,
* Peardon, Ryan, Skullerud, hep-lat/0505023), for the all-to-all propagators.
* The inverse of g5D is splitted in the contribution of the low-lying
* eigenvalues and the rest. The former part is exactly computed. The latter
* one is stochastically estimated.
* The technique of dilution is implemeted, in order to reduce the stochastic
* fluctuations.
*
*
* Parameters.
* n_eigenvalues : the number of eigenvalues needed to compute the firts part
*                 of the quark propagator
* nevt :          nevt parameter of eva
* n_global_noisy_sources_per_point : the number of stochastic sources needed
*                 to estimate the second part of the quark propagator
* omega1 :        the absolute precision required to compute the eigenvalues
* omega2 :        the relative precision required to compute the eigenvalues
* acc :           the precision to be used in the inversion routines
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
#include "random.h"
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>



/*#define TESTINGMODE*/

#define QMR_INVERTER
#define EO_PRE

#define EV_INDEX(m,a) ( (m)*nevt + (a) )

static int loglevel = 0;
static int init_flag = 0;

static int n_masses = 0;
static double hmass;
static double *shift;
static double *mass;

static const float omega1 = 1e-16;
static const float omega2 = 1e-10;
static const float acc = 1e-16;

static int n_eigenvalues = 0; /* N_{ev} */
static int nevt = 0;

static int order;
static int truncate;

static int n_global_noisy_sources_inversion;
static int n_global_noisy_sources_correction;

static int dilution;
static int n_dilution_slices;

static double *d;              /* [i*nevt+j], i<n_masses, j<n_eigenvalues */
static spinor_field **ev;       /* [i*nevt+j], i<n_masses, j<n_eigenvalues */
static spinor_field *trashevs_array;
static spinor_field *goodevs_array;
static complex *alpha;

static spinor_field *noisy_sources; /* [i], i<n_masses */
static spinor_field *noisy_sinks;   /* [i], i<n_masses */

enum {NO_DILUTION, TIME_DILUTION, TIME_SPIN_DILUTION, EXACT} dilution_mode;

typedef void (*inverting_operator)(spinor_field *source, spinor_field *sink);


#ifdef QMR_INVERTER
static spinor_field *QMR_source;
static spinor_field *QMR_sinks, *QMR_sinks_trunc;
#endif



#ifdef TESTINGMODE
static void H(spinor_field *out, spinor_field *in);
#endif

#ifdef QMR_INVERTER
static void QMR_init(int flag);
static void D_qmr(spinor_field *out, spinor_field *in);
static void D_qmr_eo(spinor_field *out, spinor_field *in);
static void D_qmr_oe(spinor_field *out, spinor_field *in);


static void create_sinks_QMR(spinor_field *source, spinor_field *sink, int max_iter);
static void create_sinks_QMR_inversion(spinor_field *source, spinor_field *sink);
static void create_sinks_QMR_correction(spinor_field *source, spinor_field *sink);
#define GET_SINKS create_sinks_QMR_inversion
#define CORRECT_SINKS create_sinks_QMR_correction

#endif


static void create_diluted_source(spinor_field *source, int di, int dilution);
/*static void project_spinor(spinor_field *sp, spinor_field **vectors);*/
static void generate_propagator(inverting_operator D, int n_src, complex** prop);
static void hopping_propagator(complex** prop);
static void add_source_sink_contraction(complex *out, spinor_field *source, spinor_field *sink, double z);



void traced_ata_qprop_dublin_trunc(complex** ev_prop, complex*** prop, int n_points) {
/*
	ev_prop[m][t*16+i]
	prop[x][m][t*16+i]
	x in [0, n_points-1]
	m in [0, n_masses-1]
	t in [0 , T-1]
	i in [0, 15]
*/

	int m, xp;
#ifdef TESTINGMODE
	int x;
	double oldhmass, norm, ave;
	complex ctmp;
#endif
	int status;

	spinor_field *test=0;

	/* allocate memory */

	test = alloc_spinor_field_f(1, &glattice);
#ifdef QMR_INVERTER
	QMR_init(0);
#endif

	/* EIGENVALUE CODE NOT WORKING */
	/*

	for(m = 0; m < n_masses; m++)
		memset(ev_prop[m], '\0', sizeof(complex)*T*16);
	*/
	/* compute the lowest n_eigenvalues eigenvalues/vectors */
	/*
	if(n_eigenvalues > 0) {
		dirac_eva(n_eigenvalues,nevt,50,5*n_eigenvalues,omega1,omega2,n_masses,mass,ev,d,&status);
	}

	for(m = 0; m < n_masses; m++)
	for(p = 0; p < n_eigenvalues; p++) {
		add_source_sink_contraction(ev_prop[m], ev[EV_INDEX(m,p)], ev[EV_INDEX(m,p)], 1.0f/d[EV_INDEX(m,p)]);
	}
	for(xp = 0; xp < n_points; xp++)
	for(m = 0; m < n_masses; m++) {
		memcpy(prop[xp][m], ev_prop[m], sizeof(complex)*T*16);
	}
	*/


	/* Generate first set of sources & sinks */

	lprintf("TRACED_ATA_QPROP_DUBLIN",loglevel,"Counter = 0\n");

	if (order>=0) hopping_propagator(prop[0]); /* Calculate hopping expansion exactly. Use order==-1 to disable. */
        generate_propagator(&GET_SINKS, n_global_noisy_sources_inversion, prop[0]); /* Invert using truncation. Use truncate==0 to disable truncation. */
        if (truncate>0) generate_propagator(&CORRECT_SINKS, n_global_noisy_sources_correction, prop[0]); /* Correct for truncation, if necessary */

	/* Generate remaining sets of independent sources & sinks.
           In the case of exact calculation (dilution==3), only 1 set of sources is needed, so simply copy exact result. */

        if(dilution==3) {
	  for(xp = 1; xp < n_points; xp++)
	    for(m = 0; m < n_masses; m++)
	      memcpy(prop[xp][m],prop[0][m],sizeof(complex)*16*T);

	} else {

	for(xp = 1; xp < n_points; xp++) {

	  lprintf("TRACED_ATA_QPROP_DUBLIN",loglevel,"Counter = %d\n",xp);

	  if (order>=0) hopping_propagator(prop[xp]);
	  generate_propagator(&GET_SINKS, n_global_noisy_sources_inversion, prop[xp]);
	  if (truncate>0) generate_propagator(&CORRECT_SINKS, n_global_noisy_sources_correction, prop[xp]);

	}

	lprintf("TRACED_ATA_QPROP_DUBLIN",loglevel,"Generated noisy sources and relative sinks.\n");

	free_spinor_field(test);

	}
}




#ifdef TESTINGMODE
static void H(spinor_field *out, spinor_field *in){
	g5Dphi(hmass,out,in);
}
#endif



#ifdef QMR_INVERTER
static void D_qmr(spinor_field *out, spinor_field *in){
	Dphi(hmass,out,in);
}
static void D_qmr_eo(spinor_field *out, spinor_field *in){
  Dphi_eopre(hmass,out,in);
}
static void D_qmr_oe(spinor_field *out, spinor_field *in){
  Dphi_oepre(hmass,out,in);
}

#endif


static void generate_propagator(inverting_operator D, int n_src, complex** prop) {

  int r, di, m;

  for(r = 0; r < n_src; r++) {

    /* generate random sources */
    
    for(di = 0; di < n_dilution_slices; di++) {
      	
      create_diluted_source(noisy_sources,di,dilution);

      for(m = 1; m < n_masses; m++)
	spinor_field_copy_f(noisy_sources+m, noisy_sources); /* Why do the noisy_sources need to be an array? */
      
#ifdef TESTINGMODE
      for(m = 0; m < n_masses; m++) {
	norm = spinor_field_sqnorm_f(noisy_sources+m);
	ave = 0.;
	_DECLARE_INT_ITERATOR(i);
	_MASTER_FOR(&glattice,i) {
	  for(x = 0; x < sizeof(suNf_spinor)/sizeof(double); x++)
	    ave += ((double*)_FIELD_AT(noisy_sources+m,i))[x];
	}
	lprintf("TRACED_ATA_QPROP_DUBLIN",loglevel,"Testing noisy sources [m=%d,r=%d,di=%d]. norm/(4*Nf*VOL3)=%f ave=%e\n",m,r,di,norm/(4*VOL3*NF),.5*ave/norm);
      }
#endif

      /* generate sinks; project sources and sinks */
      D(noisy_sources, noisy_sinks);

      /* EIGENVALUE CODE NEEDS TO BE FIXED...
      if(n_eigenvalues > 0) {
	for(m = 0; m < n_masses; m++) {
	  project_spinor(noisy_sources[m], ev + EV_INDEX(m,0));
	  project_spinor(noisy_sinks[m], ev + EV_INDEX(m,0));
	}
      }
      */  
      
#ifdef TESTINGMODE
      for(m = 0; m < n_masses; m++) {
	norm = 0.;
	for(p = 0; p < n_eigenvalues; p++) {
	  ctmp.re = spinor_field_prod_re_f(noisy_sources+m, ev[EV_INDEX(m,p)]);
	  ctmp.im = spinor_field_prod_im_f(noisy_sources+m, ev[EV_INDEX(m,p)]);
	  norm += ctmp.re*ctmp.re+ctmp.im*ctmp.im;
	}
	norm = sqrt(norm);
	lprintf("TRACED_ATA_QPROP_DUBLIN",loglevel,"Testing sources projection [m=%d,r=%d,di=%d] = %e\n",m,r,di,norm);
      }
      
      for(m = 0; m < n_masses; m++) {
	norm = 0.;
	for(p = 0; p < n_eigenvalues; p++) {
	  ctmp.re = spinor_field_prod_re_f(noisy_sinks+m, ev[EV_INDEX(m,p)]);
	  ctmp.im = spinor_field_prod_im_f(noisy_sinks+m, ev[EV_INDEX(m,p)]);
	  norm += ctmp.re*ctmp.re+ctmp.im*ctmp.im;
	}
	norm = sqrt(norm);
	lprintf("TRACED_ATA_QPROP_DUBLIN",loglevel,"Testing sinks projection [m=%d,r=%d,di=%d] = %e\n",m,r,di,norm);
      }
      
      oldhmass = hmass;
      for(m = 0; m < n_masses; m++) {
	hmass = mass[m];
	H(test,noisy_sinks+m);
	spinor_field_mul_add_assign_f(test,-1.,noisy_sources+m);
	norm=spinor_field_sqnorm_f(test);
	/*if(norm>acc)*/
	lprintf("TRACED_ATA_QPROP_DUBLIN",loglevel,"Testing inversion [m=%d,r=%d,di=%d] = %e\n",m,r,di,norm);
      }
      hmass = oldhmass;
#endif
      for(m = 0; m < n_masses; m++) {
	add_source_sink_contraction(prop[m], noisy_sources+m, noisy_sinks+m, 1.0f/n_src);
      }
      
    }
  }

}


static void hopping_propagator(complex** prop) {

  /* If order<0, do nothing */
  if (order>=0) {
    int i, m;
    /* 0th order */
    for(m = 0; m < n_masses; m++)
      for(i = 0; i < T; i++) {
	prop[m][16*i].re    =  NF*(1.0/(4.0+mass[m]));
	prop[m][16*i+5].re  =  NF*(1.0/(4.0+mass[m]));
	prop[m][16*i+10].re = -NF*(1.0/(4.0+mass[m]));
	prop[m][16*i+15].re = -NF*(1.0/(4.0+mass[m]));
      }
    
    /* 1st, 2nd, 3rd order are zero, start calculating if order>=3 */
    if (order>3) {
      
      int di=0;
      spinor_field *sinktmp;

      sinktmp = alloc_spinor_field_f(1, &glattice);

      /* Generate exact sources (not noisy!) */

      for(di = 0; di < GLB_X*GLB_Y*GLB_Z*GLB_T*sizeof(suNf_spinor)/sizeof(complex); di++) {
	create_diluted_source(noisy_sources,di,3);
	spinor_field_g5_f(noisy_sources, noisy_sources);
	for(m = 0; m < n_masses; m++) {
	  /*if (m!=0) spinor_field_copy_f(noisy_sources+m,noisy_sources);*/
	  /*	spinor_field_copy_f(sinktmp+m, noisy_sources+m);*/	
	  
	  spinor_field_mul_f(sinktmp,1.0/(4.0+mass[m]),noisy_sources);
	  spinor_field_g5_f(noisy_sources, noisy_sources);
	  for(i=0; i<order; i++) {
	    
	    Dphi(-4.0, noisy_sinks+m, sinktmp);
	    
	    spinor_field_mul_f(noisy_sinks+m,1.0/(4.0+mass[m]),noisy_sinks+m);
	    spinor_field_copy_f(sinktmp, noisy_sinks+m);
	    
	    /* FIX EIGENVALUES
	       if(n_eigenvalues > 0) {
	        for(m = 0; m < n_masses; m++) {
	         project_spinor(noisy_sources[m], ev + EV_INDEX(m,0));
	         project_spinor(noisy_sinks[m], ev + EV_INDEX(m,0));
	        }
	       }
	    */
	    spinor_field_g5_f(noisy_sources, noisy_sources);
	      
	    if ((i%2==1) && (i!=1)) add_source_sink_contraction(prop[m], noisy_sources, noisy_sinks+m, 1.0);
	  }
	}
      }
           
      free_spinor_field(sinktmp);
    }
  }
}


void ata_qprop_dublin_trunc_init(int nm, double *mptr, int pnev, int pnevt, int ord, int trnc, int pgns, int pgnsc, int pdil) {
  int m;
  double requiredmemory = 0.0;
  double reqmem_for_evs = 0.0;
  double reqmem_for_sources = 0.0;
  double reqmem_for_sinks = 0.0;
  
  if(init_flag != 0) return;

  /* static parameters */
  
  n_masses = nm;
  error(nm<1,1,"ata_qprop_dublin_init [trunc_hairpin.c]", "Bad choice for n_masses");
  mass=(double*)malloc(sizeof(double)*n_masses);
  shift=(double*)malloc(sizeof(double)*n_masses);
  hmass = mptr[0]; /* we can put any number for the index! */
  for(m = 0; m < n_masses; ++m) {
    mass[m] = mptr[m];
    shift[m] = hmass - mass[m];
  }
  
  n_eigenvalues = pnev;
  nevt = pnevt;
  error(pnev<0,1,"ata_qprop_dublin_init [trunc_hairpin.c]", "Bad choice for n_eigenvalues");
  error(pnevt<pnev,1,"ata_qprop_dublin_init [trunc_hairpin.c]", "Bad choice for nevt");
  
  dilution = pdil;
  if(dilution == NO_DILUTION)
    n_dilution_slices = 1;
  else if(dilution == TIME_DILUTION)
    n_dilution_slices = GLB_T;
  else if(dilution == TIME_SPIN_DILUTION)
    n_dilution_slices = 4*GLB_T;
  else if(dilution == EXACT)
    n_dilution_slices = GLB_X*GLB_Y*GLB_Z*GLB_T*sizeof(suNf_spinor)/sizeof(complex); /* variable name for total global volume? */
  else
    error(0==0,1,"ata_qprop_dublin_init [trunc_hairpin.c]", "Bad choice for dilution");

  truncate = trnc;
  error(truncate<0,1,"ata_qprop_dublin_init [trunc_hairpin.c]", "Bad choice for truncate");

  /* Use order=-1 to disable hopping expansion (order=0 is 0th order) */
  order = ord;
  error(order<-1,1,"ata_qprop_dublin_init [trunc_hairpin.c]", "Bad choice for order");

  n_global_noisy_sources_inversion = pgns;
  error(n_global_noisy_sources_inversion<0,1,"ata_qprop_dublin_init [trunc_hairpin.c]", "Bad choice for n_global_noisy_sources_inversion");

  n_global_noisy_sources_correction = pgnsc;
  error(n_global_noisy_sources_correction<0,1,"ata_qprop_dublin_init [trunc_hairpin.c]", "Bad choice for n_global_noisy_sources_correction");

  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Number of masses = %d\n",n_masses);
  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Highest mass = %f\n",hmass);
  for(m = 0;m < n_masses; m++){
    lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Shifts: %f\n",shift[m]);
  }
  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Number of eigenvalues = %d (%d)\n",n_eigenvalues,nevt);
  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Order of hopping expansion = %d\n",order);
  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Iterations after which to truncate = %d\n",truncate);
  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Number of global noisy sources for (truncated) inversion = %d\n",n_global_noisy_sources_inversion);
  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Number of global noisy sources for correction to truncation = %d\n",n_global_noisy_sources_correction);
  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Number of dilution slices = %d\n",n_dilution_slices);
  lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel+1,"Dilution level = %d\n",dilution);

  
   /* required memory --NOT UPDATED FOR TRUNCATION CODE OR MPI-- */
   
   reqmem_for_evs = sizeof(suNf_spinor)*VOLUME*n_eigenvalues*n_masses
     +sizeof(suNf_spinor)*VOLUME*(nevt-n_eigenvalues)
     +sizeof(complex)*n_eigenvalues;
   reqmem_for_sources = sizeof(suNf_spinor)*VOLUME*n_masses;
   reqmem_for_sinks = sizeof(suNf_spinor)*VOLUME*n_masses;
   requiredmemory = reqmem_for_evs + reqmem_for_sources + reqmem_for_sinks;

	/* EIGENVALUE CODE BROKEN - TO BE FIXED */

	/* eigenvalues */
	/*	if(nevt != 0)
   	d = (double*)malloc(sizeof(double)*n_masses*nevt);
   else
      d = NULL;
	*/
	/* eigenvectors */
	/*	if(n_eigenvalues != nevt)
   	trashevs_array = alloc_spinor_field_f(nevt-n_eigenvalues);
   else
      trashevs_array = NULL;

	if(nevt != 0) {
	   ev = (suNf_spinor**)malloc(sizeof(suNf_spinor*)*n_masses*nevt);
	   for(m = 0; m < n_masses; m++) {
		   for(a = 0; a < n_eigenvalues; a++)
			   ev[EV_INDEX(m,a)] = goodevs_array+(m*n_eigenvalues+a);
		   for(a = n_eigenvalues; a < nevt; a++)
			   ev[EV_INDEX(m,a)] = trashevs_array+(a-n_eigenvalues);
	   }
	} else
      ev = NULL;

	if(n_eigenvalues != 0)
   	goodevs_array = alloc_spinor_field_f(n_masses*n_eigenvalues);
   else
      goodevs_array = NULL;
	*/
	/* eigenvectors related */
	/*	if(n_eigenvalues != 0)
   	alpha = (complex*)malloc(sizeof(complex)*n_eigenvalues);
   else
      alpha = NULL;
	*/


   /* allocate sources & sinks */
   noisy_sources = alloc_spinor_field_f(n_masses,&glattice);
   noisy_sinks   = alloc_spinor_field_f(n_masses,&glattice);

#ifdef QMR_INVERTER
   QMR_source = alloc_spinor_field_f(1, &glattice);
   QMR_sinks  = alloc_spinor_field_f(n_masses, &glattice);
   QMR_sinks_trunc =  alloc_spinor_field_f(n_masses, &glattice);

   /* --MEMORY USE NOT UPDATED-- */		
   requiredmemory += 2*(1+n_masses)*sizeof(suNf_spinor)*VOLUME;
#endif
	
   requiredmemory /= 1024.0*1024.0;
   lprintf("ALL_TO_ALL_QUARK_PROPAGATOR_INIT",loglevel,"Required memory = %f Mb\n",requiredmemory);
	
   init_flag = 1;
}



void ata_qprop_dublin_trunc_free() {
	if(init_flag != 1) return;
	
	/* eigenvalues 
	if(nevt != 0) free(d);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed eigenvalues\n");
	*/
	/* eigenvectors 
	if(n_eigenvalues != nevt) free_field(trashevs_array);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed trash eigenvectors\n");
	if(n_eigenvalues != 0) free_field(goodevs_array);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed good eigenvectors\n");
	if(nevt != 0) free(ev);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed eigenvectors\n");
	*/
	/* eigevectors related 
	if(n_eigenvalues != 0) free(alpha);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed alphas\n");
	*/
	/* sources */
	free_spinor_field(noisy_sources);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed noisy sources\n");

	/* sinks */
	free_spinor_field(noisy_sinks);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed noisy sinks\n");

#ifdef QMR_INVERTER
	free_spinor_field(QMR_source);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed QMR sources\n");
	free_spinor_field(QMR_sinks);
	free_spinor_field(QMR_sinks_trunc);
	lprintf("DUBLIN_MESONS_FREE_MEMORY",loglevel+10,"Freed QMR sinks\n");
#endif
	
	init_flag = 0;
}

void ata_qprop_dublin_trunc_modify_par(int ord, int trnc, int pgns, int pgnsc, int pdil) {
  truncate = trnc;
  error(truncate<0,1,"ata_qprop_dublin_trunc_modify_par [trunc_hairpin.c]", "Bad choice for truncate");
  order = ord;
  error(order<-1,1,"ata_qprop_dublin_trunc_modify_par [trunc_hairpin.c]", "Bad choice for order");
  n_global_noisy_sources_inversion = pgns;
  error(n_global_noisy_sources_inversion<0,1,"ata_qprop_dublin_trunc_modify_par [trunc_hairpin.c]", "Bad choice for n_global_noisy_sources_inversion");
  n_global_noisy_sources_correction = pgnsc;
  error(n_global_noisy_sources_correction<0,1,"ata_qprop_dublin_trunc_modify_par [trunc_hairpin.c]", "Bad choice for n_global_noisy_sources_correction");
  dilution = pdil;
  if(dilution == NO_DILUTION)
    n_dilution_slices = 1;
  else if(dilution == TIME_DILUTION)
    n_dilution_slices = GLB_T;
  else if(dilution == TIME_SPIN_DILUTION)
    n_dilution_slices = 4*GLB_T;
  else if(dilution == EXACT)
    n_dilution_slices = GLB_X*GLB_Y*GLB_Z*GLB_T*sizeof(suNf_spinor)/sizeof(complex); /* variable name for total global volume? */
  else
    error(0==0,1,"generate_propagator [trunc_hairpin.c]", "Bad choice for dilution");
}

#ifdef QMR_INVERTER
static void QMR_init(int flag) {
	static int init_flag = 0;
	mshift_par QMR_par;
	int m;
	int cgiter=0;
	
	if(init_flag == 0 || flag != 0) {
		gaussian_spinor_field(QMR_source);
		init_flag = 1;
	}

	/* set up inverters parameters */
	QMR_par.n = n_masses;
	QMR_par.shift = shift;
	QMR_par.err2 = .5*acc;
	QMR_par.max_iter = 0;

	for(m = 0; m < n_masses; m++) {
		spinor_field_zero_f(QMR_sinks+m);
		spinor_field_zero_f(QMR_sinks_trunc+m);
	}

	spinor_field_g5_f(QMR_source, QMR_source);
	cgiter+=g5QMR_mshift_trunc(&QMR_par, truncate, &D_qmr, QMR_source, QMR_sinks_trunc, QMR_sinks);
	if (truncate==0) {
	  for(m = 0; m < n_masses; m++) {
	    spinor_field_copy_f(QMR_sinks_trunc+m,QMR_sinks+m);
	  }
	}

	spinor_field_g5_f(QMR_source, QMR_source);

	lprintf("QMR_INIT",loglevel+1,"QMR MVM = %d\n",cgiter);	
}
#endif


/* Not fully tested yet with dilution!=3 */
static void create_diluted_source(spinor_field *source, int di, int dilution) {
	
	if(dilution == NO_DILUTION) {
	  _DECLARE_INT_ITERATOR(i);
	  _MASTER_FOR(&glattice,i) {
	    ranz2((double*)_FIELD_AT(source,i),sizeof(suNf_spinor)/sizeof(double));
	  }
	} else if(dilution == TIME_DILUTION) {
	   int i, x[4];
	   spinor_field_zero_f(source);
	   x[0]=di;
	   for(i = 0; i < GLB_X*GLB_Y*GLB_Z; i++) {
	     x[1]=i % GLB_X;
	     x[2]=(i/GLB_X) % GLB_Y;
	     x[3]=(i/GLB_X/GLB_Y) % GLB_Z;
	     if(COORD[0]==x[0]/T && COORD[1]==x[1]/X && COORD[2]==x[2]/Y && COORD[3]==x[3]/Z) {
	       ranz2((double*)_FIELD_AT(source,ipt(x[0]%T, x[1]%X, x[2]%Y, x[3]%Z)),sizeof(suNf_spinor)/sizeof(double)); }
	   }
	} else if(dilution == TIME_SPIN_DILUTION) {
	  int i, x[4];
	  spinor_field_zero_f(source);
	  x[0] = di/4;
	  for(i = 0; i < GLB_X*GLB_Y*GLB_Z; i++) {
	    x[1]=i % GLB_X;
	    x[2]=(i/GLB_X) % GLB_Y;
	    x[3]=(i/GLB_X/GLB_Y) % GLB_Z;
	    if(COORD[0]==x[0]/T && COORD[1]==x[1]/X && COORD[2]==x[2]/Y && COORD[3]==x[3]/Z)
	      ranz2((double*)(&(_FIELD_AT(source,ipt(x[0]%T, x[1]%X, x[2]%Y, x[3]%Z)))->c[di%4]),sizeof(suNf_vector)/sizeof(double));
	  }
	} else if(dilution == EXACT) {
	  int site=di/(4*NF);
	  int component=di%(4*NF);
	  int x[4];
	  x[0]=site % GLB_T;
	  x[1]=(site/GLB_T) % GLB_X;
	  x[2]=(site/GLB_T/GLB_X) % GLB_Y;
	  x[3]=(site/GLB_T/GLB_X/GLB_Y) % GLB_Z;
	  spinor_field_zero_f(source);
	  if(COORD[0]==x[0]/T && COORD[1]==x[1]/X && COORD[2]==x[2]/Y && COORD[3]==x[3]/Z)
	    ((complex*)_FIELD_AT(source, ipt(x[0]%T, x[1]%X, x[2]%Y, x[3]%Z)))[component].re = 1.;  
	}
}


#ifdef QMR_INVERTER
static void create_sinks_QMR(spinor_field *source, spinor_field *sink, int max_iter) {
	mshift_par QMR_par;
	int i,m;
	int cgiter=0;

	spinor_field *sinktmp, *sink_trunc;

	/* set up inverters parameters */
	QMR_par.n = n_masses;
	QMR_par.shift = shift;
	QMR_par.err2 = .5*acc;
	QMR_par.max_iter = max_iter;

	sink_trunc=alloc_spinor_field_f(n_masses,&glattice);
	sinktmp=alloc_spinor_field_f(1,&glattice);

	for(m = 0; m < n_masses; m++) {
	  spinor_field_zero_f(sink + m);
	  spinor_field_zero_f(sink_trunc + m);
	}

	spinor_field_add_assign_f(source, QMR_source);
	
	spinor_field_g5_f(source,source);

#ifdef EO_PRE

	spinor_field *b, b_even, b_odd, sinktmp_even, sinktmp_odd, source_even, source_odd, *sink_even, *sink_odd, *sink_trunc_even, *sink_trunc_odd;
        double* shift_eo;

	/* Prepare even/odd masks */	
	b=alloc_spinor_field_f(1, &glattice);
	b_even=*b;
	b_even.type=&glat_even;
	b_odd=*b;
	b_odd.type=&glat_odd;

	sinktmp_even=*sinktmp;
	sinktmp_even.type=&glat_even;
	sinktmp_odd=*sinktmp;
	sinktmp_odd.type=&glat_odd;

	source_even=*source;
	source_even.type=&glat_even;
	source_odd=*source;
	source_odd.type=&glat_odd;
	
	sink_even=(spinor_field*)malloc(sizeof(spinor_field)*n_masses);
	sink_odd=(spinor_field*)malloc(sizeof(spinor_field)*n_masses);
	sink_trunc_even=(spinor_field*)malloc(sizeof(spinor_field)*n_masses);
	sink_trunc_odd=(spinor_field*)malloc(sizeof(spinor_field)*n_masses);

	for(m = 0; m < n_masses; m++) {
	  *(sink_even+m)=*(sink+m);
	  (sink_even+m)->type=&glat_even;
	  *(sink_odd+m)=*(sink+m);
	  (sink_odd+m)->type=&glat_odd;

	  *(sink_trunc_even+m)=*(sink_trunc+m);
	  (sink_trunc_even+m)->type=&glat_even;
	  *(sink_trunc_odd+m)=*(sink_trunc+m);
	  (sink_trunc_odd+m)->type=&glat_odd;
	
	}

        /* Start preconditioning & inversion */

	if (n_masses==1) {

	  /* Non-multishift EO preconditioning, following DeGrand & DeTar p. 174-175.
             Requires only one non-trivial inversion.                                  */
	  
	  Dphi_(&b_even, &source_odd);
          spinor_field_mul_f(&b_even, -1.0/(4.0+mass[0]), &b_even);
          spinor_field_add_assign_f(&b_even, &source_even);
	  /* NOT OK -> spinor_field_copy_f(&b_odd, &source_odd); */
	  /* SLOW?  -> spinor_field_mul_f(&b_odd, 1.0, &source_odd); */
	  spinor_field_zero_f(&b_odd); /* Is this better? */
	  spinor_field_add_assign_f(&b_odd, &source_odd);

          cgiter+=g5QMR_mshift_trunc(&QMR_par, truncate, &D_qmr_eo, &b_even, sink_trunc_even, sink_even);

          spinor_field_mul_f(sink_even,(4.0+mass[0]),sink_even);
          spinor_field_mul_f(sink_odd,1.0/(4.+mass[0]),sink_odd);
          Dphi_(sink_odd, sink_even);
          spinor_field_sub_assign_f(sink_odd, &b_odd);
          spinor_field_mul_f(sink_odd, -1.0/(4.0+mass[0]), sink_odd);

	  if(max_iter!=truncate) {
	    spinor_field_mul_f(sink_trunc_even,(4.0+mass[0]),sink_trunc_even);
	    spinor_field_mul_f(sink_trunc_odd,1.0/(4.+mass[0]),sink_trunc_odd);
	    Dphi_(sink_trunc_odd, sink_trunc_even);
	    spinor_field_sub_assign_f(sink_trunc_odd, &b_odd);
	    spinor_field_mul_f(sink_trunc_odd, -1.0/(4.0+mass[0]), sink_trunc_odd);
	  }
	  
	}
	
	else {

	  /* Multishift preconditioning scheme, requires two non-trivial inversions. */
	  
	  shift_eo=(double*)malloc(sizeof(double)*n_masses);
	  for(i=0; i<n_masses; i++)
	    shift_eo[i]=(4.+mass[0])*(4.+mass[0])-(4.+mass[i])*(4.+mass[i]);
	  QMR_par.shift=shift_eo;
	  
	  cgiter+=g5QMR_mshift_trunc(&QMR_par, truncate, &D_qmr_eo, &source_even, sink_trunc_even, sink_even);
	  cgiter+=g5QMR_mshift_trunc(&QMR_par, truncate, &D_qmr_oe, &source_odd, sink_trunc_odd, sink_odd);
	  
          free(shift_eo);

          for(m = 0; m < n_masses; m++) {
            spinor_field_copy_f(b,sink+m);
            spinor_field_mul_f(sink_even+m, 4.0+mass[m], &b_even);
            Dphi_(&sinktmp_even, &b_odd);
            spinor_field_sub_assign_f(sink_even+m, &sinktmp_even);
            spinor_field_mul_f(sink_odd+m, 4.0+mass[m], &b_odd);
            Dphi_(&sinktmp_odd, &b_even);
            spinor_field_sub_assign_f(sink_odd+m, &sinktmp_odd);

	    if(max_iter!=truncate) {
	      spinor_field_copy_f(b,sink_trunc+m);
	      spinor_field_mul_f(sink_trunc_even+m, 4.0+mass[m], &b_even);
	      Dphi_(&sinktmp_even, &b_odd);
	      spinor_field_sub_assign_f(sink_trunc_even+m, &sinktmp_even);
	      spinor_field_mul_f(sink_trunc_odd+m, 4.0+mass[m], &b_odd);
	      Dphi_(&sinktmp_odd, &b_even);
	      spinor_field_sub_assign_f(sink_trunc_odd+m, &sinktmp_odd);
	    }
          }
	}

	free_spinor_field(b);
	free(sink_even);
	free(sink_odd);
	free(sink_trunc_even);
	free(sink_trunc_odd);
	    
#else
	cgiter+=g5QMR_mshift_trunc(&QMR_par, truncate, &D_qmr, source, sink_trunc, sink);
#endif

	spinor_field_g5_f(source,source);

	for(m = 0; m < n_masses; m++) {

	  /* If truncation is present (truncate>0) and the inverter also runs to convergence
	     (max_iter==0), then this is the correction term, and the difference "exact minus
             trunc" should be computed. */
	  if ((truncate>0) && (max_iter==0)) {
	    spinor_field_sub_assign_f(sink + m, QMR_sinks + m);
	    spinor_field_sub_assign_f(sink_trunc + m, QMR_sinks_trunc + m);
	    spinor_field_sub_assign_f(sink + m, sink_trunc + m); 
	  }
	  else {
	    /* Otherwise only the truncated sink is needed.
               NB. If truncate==max_iter then sink==sink_trunc, so below statement is OK.
                   If no truncation, QMR_sinks==QMR_sinks_trunc, so below statement is again OK. */
	    spinor_field_sub_f(sink + m, sink + m, QMR_sinks_trunc + m);
	  }

	  for(i = 0; i<=order; i++) {
	    Dphi(-4.0, sinktmp, sink + m);
	    spinor_field_copy_f(sink + m, sinktmp);
	  }

	  spinor_field_mul_f(sink + m, -1.0/(4.0+mass[m])*pow(-1.0/(4.0+mass[m]), order), sink + m);

	}

	spinor_field_sub_assign_f(source, QMR_source);

	free_spinor_field(sinktmp);
	
	lprintf("GET_SINKS_QMR",loglevel+1,"QMR MVM = %d\n",cgiter);

}

static void create_sinks_QMR_inversion(spinor_field *source, spinor_field *sink) {
  create_sinks_QMR(source, sink, truncate);
}

static void create_sinks_QMR_correction(spinor_field *source, spinor_field *sink) {
    create_sinks_QMR(source, sink, 0);
}

#endif

/* NOT UPDATED YET
static void project_spinor(suNf_spinor *sp, suNf_spinor **vectors) {
	int a, x;
	
	for(a = 0; a < n_eigenvalues; a++)
		alpha[a].re = alpha[a].im = 0.0;

	for(x = 0; x < VOLUME; x++)
	for(a = 0; a < n_eigenvalues; a++) {
		_spinor_prod_assign_f(alpha[a],vectors[a][x],sp[x]);
	}

	for(a = 0; a < n_eigenvalues; a++)
	for(x = 0; x < VOLUME; x++) {
		_spinor_project_f(sp[x],alpha[a],vectors[a][x]);
	}
}
*/


static void add_source_sink_contraction(complex *out, spinor_field *source, spinor_field *sink, double z) {
  	int i, j, t, x, index;
	suNf_vector *eta, *csi;
	complex tmp;
#ifdef TESTINGMODE
	complex trace, g5trace;
	double norm;
#endif
	
#ifdef TESTINGMODE
	trace.re = trace.im = 0.;
	g5trace.re = g5trace.im = 0.;
#endif

	int point[4];

	for(t = 0; t < GLB_T; t++) {
	  for(x = 0; x < GLB_X*GLB_Y*GLB_Z; x++) {

	    point[0]=t;
	    point[1]=(x/GLB_X/GLB_Y) % GLB_Z;
	    point[2]=(x/GLB_X) % GLB_Y;
	    point[3]=x % GLB_X;
	    if(COORD[0]==point[0]/T && COORD[1]==point[1]/X && COORD[2]==point[2]/Y && COORD[3]==point[3]/Z) {
	      index = ipt(point[0]%T, point[1]%X, point[2]%Y, point[3]%Z);  

	      for(i = 0; i < 4; i++) {
		csi = (suNf_vector*)(_FIELD_AT(sink,index)) + i;
		for(j = 0; j < 4; j++) {
		  eta = (suNf_vector*)(_FIELD_AT(source,index)) + j;
		  tmp.re = tmp.im = 0.;
		  _vector_prod_assign_f(tmp, *eta, *csi);
		  out[t*16+SPIN_2D_INDEX(i,j)].re += tmp.re*z/VOL3;
		  out[t*16+SPIN_2D_INDEX(i,j)].im += tmp.im*z/VOL3;
#ifdef TESTINGMODE
		  if(i==j) {
		    trace.re += tmp.re;
		    trace.im += tmp.im;
		    g5trace.re += (i==0 || i==1) ? tmp.re : -tmp.re;
		    g5trace.im += (i==0 || i==1) ? tmp.im : -tmp.im;
		  }
#endif
		}
	      }
	    }
	  }
	}
  
	lprintf("ADD_SOURCE_SINK_CONTRACTION",loglevel+2,"Written in %p\n",out);

#ifdef TESTINGMODE
	trace.re -= spinor_field_prod_re_f(source,sink);
	trace.im -= spinor_field_prod_im_f(source,sink);
	g5trace.re -= spinor_field_g5_prod_re_f(source,sink);
	/*  Not defined! 
            g5trace.im -= spinor_field_g5_prod_im_f(source,sink);
                                                                       */
	norm = sqrt( trace.re*trace.re + trace.im*trace.im );
	lprintf("ADD_SOURCE_SINK_CONTRACTION",loglevel,"Testing trace = %e\n",norm);
	norm = sqrt( g5trace.re*g5trace.re + g5trace.im*g5trace.im );
	lprintf("ADD_SOURCE_SINK_CONTRACTION",loglevel,"Testing g5trace = %e\n",norm);
#endif

}
