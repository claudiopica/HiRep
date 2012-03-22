/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "dirac.h"
#include "representation.h"
#include "linear_algebra.h"
#include "logger.h"
#include "communications.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* State quantities for HMC */
extern suNg_av_field *momenta;
extern spinor_field *pf;
extern rhmc_par _update_par;
extern rational_app r_MD; /* used in the action MD evolution */
extern spinor_operator H2;
extern spinor_operator H;
extern integrator_par *integrator;
/* extern double minev, maxev; */ /* min and max eigenvalue of H^2 */


static short int init=0;

static suNg_field *u_gauge_old=NULL;
static scalar_field *la=NULL; /* local action field for Metropolis test */
static MINRES_par pfa;

void init_hmc(rhmc_par *par){
    
	if (init) {
	  /* already initialized */
	  lprintf("HMC",0,"WARNING: HMC already initialized!\nWARNING: Ignoring call to init_hmc.\n");
		return;
	}
    
	lprintf("HMC",0,"Initializing...\n");
    
	/* copy input parameters into the internal variable and make some tests */
	_update_par=*par;
	if (_update_par.nf%2 != 0) {
		lprintf("HMC",0,"The number of fermions even.\nTry with the RHMC algorithm.\n");
		error(1,1,"init_hmc","The HMC algorithm is not suitable for the parameters specified as input.\n");
	}
    
	lprintf("HMC",10,
	  "Number of Flavors = %d\n"
	  "beta = %.8f\n"
	  "Mass = %.8f\n"
	  "Metropolis test precision = %.8e\n"
	  "RHMC force precision = %.8e\n"
	  "MD trajectory length = %.8f\n"
	  "MD steps = %d\n"
	  "MD gauge substeps = %d\n"
	  ,_update_par.nf
	  ,_update_par.beta
	  ,_update_par.mass
	  ,_update_par.MT_prec
	  ,_update_par.force_prec
	  ,_update_par.tlen
	  ,_update_par.nsteps
	  ,_update_par.gsteps
	  );
	
	_update_par.n_pf = _update_par.nf/2;
    
	/* allocate space for the backup copy of gfield */
	if(u_gauge_old==NULL) u_gauge_old=alloc_gfield(&glattice);
	suNg_field_copy(u_gauge_old,u_gauge);
    
	/* allocate momenta */
	if(momenta==NULL) momenta = alloc_avfield(&glattice);
    
	/* allocate pseudofermions */
	/* we allocate one more pseudofermion for the computation 
	* of the final action density 
	*/
	if(pf==NULL) {
	  /* we need 1 more spinor field for Metropolis test action with MINRES */
	  pf=alloc_spinor_field_f(_update_par.n_pf+1,
#ifdef UPDATE_EO
      &glat_even /* even lattice for preconditioned dynamics */
#else
      &glattice /* global lattice */
#endif 
    );
  }
    
  /* allocate memory for the local action */
  if(la==NULL) la=alloc_sfield(&glattice);

  /* represent gauge field */
  represent_gauge_field();

  /* integrator */
  integrator = (integrator_par*)malloc(sizeof(integrator_par)*3);

  integrator[0].level = 0;
  integrator[0].tlen = _update_par.tlen;
  integrator[0].nsteps = _update_par.nsteps;
  integrator[0].force = &force_hmc;
  integrator[0].force_par = malloc(sizeof(force_hmc_par));
  ((force_hmc_par*)(integrator[0].force_par))->n_pf = _update_par.nf/2;
  ((force_hmc_par*)(integrator[0].force_par))->pf = pf;
  ((force_hmc_par*)(integrator[0].force_par))->hasenbusch = 0;
  ((force_hmc_par*)(integrator[0].force_par))->inv_err2 = _update_par.force_prec;
  ((force_hmc_par*)(integrator[0].force_par))->inv_err2_flt = _update_par.force_prec_flt;
  integrator[0].integrator = &O2MN_multistep;
  integrator[0].next = &integrator[1];

  integrator[1].level = 1;
  integrator[1].tlen = integrator[0].tlen/((double)(2*integrator[0].nsteps));
  integrator[1].nsteps = _update_par.gsteps;
  integrator[1].force = &force0;
  integrator[1].force_par = NULL;
  integrator[1].integrator = &O2MN_multistep;
  integrator[1].next = &integrator[2];

  integrator[2].level = 2;
  integrator[2].tlen = integrator[1].tlen/((double)(2*integrator[1].nsteps));
  integrator[2].nsteps = 1;
  integrator[2].force = NULL;
  integrator[2].force_par = NULL;
  integrator[2].integrator = &gauge_integrator;
  integrator[2].next = NULL;
  
  init_force_hmc();

    
  /* set up rational approx needed for HMC */
  pfa.err2=_update_par.MT_prec;
  pfa.err2*=pfa.err2;
  pfa.max_iter=0;
  
  init = 1;
	
  lprintf("HMC",0,"Initialization done.\n");
    
}

void free_hmc(){
   
  if (!init) {
    /* not initialized */
    lprintf("HMC",0,"WARNING: HMC not initialized!\nWARNNG: Ignoring call to free_hmc.\n");
    return;
  }
  
  /* free momenta */
  if(u_gauge_old!=NULL) free_gfield(u_gauge_old); u_gauge_old=NULL;
  if(momenta!=NULL) free_avfield(momenta); momenta=NULL;
  if(pf!=NULL) free_spinor_field_f(pf); pf=NULL;
  
  if(la!=NULL) free_sfield(la); la=NULL;
  
  if (integrator!=NULL) {
    int i=0;
    do {
      if(integrator[i].force_par != NULL) free(integrator[i].force_par);
    } while (integrator[i++].next != NULL);
    free(integrator);    
  }
  
  free_force_hmc();
  
  init = 0;
  
  lprintf("HMC",0,"Memory deallocated.\n");
    
}

int update_hmc(){
  double deltaH;
    /* double maxev,minev; */
    _DECLARE_INT_ITERATOR(i);
    
    if(!init) {
        /* not initialized */
        lprintf("HMC",0,"WARNING: HMC not initialized!\nWARNNG: Ignoring call to update_hmc.\n");
        return -1;
    }
        
    /* generate new momenta and pseudofermions */
    lprintf("HMC",30,"Generating gaussian momenta and pseudofermions...\n");
    gaussian_momenta(momenta);
    for (i=0;i<_update_par.n_pf;++i)
      gaussian_spinor_field(&pf[i]);
    
    /* compute starting action */
    lprintf("HMC",30,"Computing action density...\n");
    //    local_hmc_action_cpu(NEW, la, momenta, pf, pf);
    local_hmc_action(NEW, la, momenta, pf, pf);

    /* compute H2^{1/2}*pf = H*pf */
    lprintf("HMC",30,"Correcting pseudofermions distribution...\n");
#ifdef WITH_GPU //Make sure gauge field is on GPU
    gfield_copy_to_gpu_f(u_gauge_f); 
    gfield_copy_to_gpu_f_flt(u_gauge_f_flt);
#endif
    for (i=0;i<_update_par.n_pf;++i) {
        spinor_field_copy_f(&pf[_update_par.n_pf],&pf[i]);
        H.dbl(&pf[i], &pf[_update_par.n_pf]);
    }
    
    /* integrate molecular dynamics */
    lprintf("HMC",30,"MD integration...\n");
    (*(integrator[0].integrator))(momenta,&integrator[0]);
    
    /* project gauge field */
    project_gauge_field();
    represent_gauge_field();
    
    lprintf("HMC",30,"Computing new action density...\n");
    /* compute H2^{-1/2}*pf or H2^{-1}*pf */
    /* here we choose the first strategy which is more symmetric */
    /* for the HMC H2^-1/2 = H^-1 and we use MINRES for this inversion */
#ifdef WITH_GPU //Make sure gauge field is on GPU
    gfield_copy_to_gpu_f(u_gauge_f); 
    gfield_copy_to_gpu_f_flt(u_gauge_f_flt);
#endif
    for (i=0;i<_update_par.n_pf;++i) {
      spinor_field_copy_f(&pf[_update_par.n_pf],&pf[i]);
      MINRES(&pfa,H,&pf[_update_par.n_pf],&pf[i],0);
    }
    
    /* compute new action */
    local_hmc_action(DELTA, la, momenta, pf, pf);
    //    local_hmc_action_cpu(DELTA, la, momenta, pf, pf);
    
    /* Metropolis test */

#ifdef WITH_GPU
    deltaH = scalar_field_sum(la);
#else
    deltaH=0.;
    _MASTER_FOR(la->type,i) {
      deltaH+=*_FIELD_AT(la,i);
    }
#endif
    global_sum(&deltaH, 1);
    lprintf("HMC",10,"[DeltaS = %1.8e][exp(-DS) = %1.8e]\n",deltaH,exp(-deltaH));
    
    if(deltaH<0.) {
        suNg_field_copy(u_gauge_old,u_gauge);
    } else {
        double r;
        if (PID==0) { ranlxd(&r,1); if(r<exp(-deltaH)) r=1.0; else r=-1.0;}  /* make test on on PID 0 */
        bcast(&r,1);
        if(r>0.) {
            suNg_field_copy(u_gauge_old,u_gauge);
        } else {
            lprintf("HMC",10,"Configuration rejected.\n");
            suNg_field_copy(u_gauge,u_gauge_old);
            start_gf_sendrecv(u_gauge); /* this may not be needed if we always guarantee that we copy also the buffers */
            represent_gauge_field();         
            return 0;
        }
    }
    
    lprintf("HMC",10,"Configuration accepted.\n");
    
    return 1;
}


