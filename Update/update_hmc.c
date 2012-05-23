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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"

/* State quantities for HMC */
suNg_av_field *momenta=NULL;
spinor_field *pf=NULL;
hmc_par _update_par={0};
integrator_par *integrator = NULL;
double minev, maxev; /* min and max eigenvalue of H^2 */
int n_pf;
action_par hmc_action_par;
/* END of State */

static double static_mass=0.;

static void D(spinor_field *out, spinor_field *in){
  #ifdef UPDATE_EO
  Dphi_eopre(static_mass, out, in);
  #else
  Dphi(static_mass, out, in);
  #endif
}


static void D_flt(spinor_field_flt *out, spinor_field_flt *in){
  #ifdef UPDATE_EO
  Dphi_eopre_flt((float)(static_mass), out, in);
  #else
  Dphi_flt((float)(static_mass), out, in);
  #endif
}


static short int init=0;

static suNg_field *u_gauge_old=NULL;
static scalar_field *la=NULL; /* local action field for Metropolis test */

void init_hmc(hmc_par *par){
  
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
		"Metropolis test precision (F1==n) = %.8e\n"
		"HMC force precision (F1==n) = %.8e\n"
		"Metropolis test precision (F1==n) float = %.8e\n"
		"HMC force precision (F1==n) float = %.8e\n"
		,_update_par.nf
		,_update_par.beta
		,_update_par.mass
		,_update_par.n_MT_prec
		,_update_par.n_force_prec
		,_update_par.n_MT_prec_flt
		,_update_par.n_force_prec_flt
		);
	
	if(_update_par.hasenbusch)
	lprintf("HMC",10,
	  "Using Hasenbush accelerator\n"
	  "Metropolis test precision (F2==h) = %.8e\n"
	  "HMC force precision (F2==h) = %.8e\n"
	  "Hasenbush mass shift = %.8e\n"
	  ,_update_par.h_MT_prec
	  ,_update_par.h_force_prec
	  ,_update_par.hasen_dm
	  );
	
	lprintf("HMC",10,
		"MD trajectory length = %.8f\n"
		"MD steps = %d\n"
		"MD gauge substeps = %d\n"
		,_update_par.tlen
		,_update_par.nsteps
		,_update_par.gsteps
		);
	
	if(_update_par.hasenbusch)
	lprintf("HMC",10,
	  "MD Hasenbusch sub steps = %d\n"
	  ,_update_par.hsteps
	  );
	
	
	
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
	  n_pf = _update_par.hasenbusch?_update_par.nf:(_update_par.nf/2);
	  /* we need 1 more spinor field for Metropolis test action with MINRES */
	  pf=alloc_spinor_field_f(n_pf+1,
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
  
  int currentlevel;
  
  /* integrator */
  integrator = (integrator_par*)malloc(sizeof(integrator_par)*(_update_par.hasenbusch?4:3));
  
  currentlevel = 0;
  integrator[currentlevel].level = currentlevel;
  integrator[currentlevel].nsteps = _update_par.nsteps;
  integrator[currentlevel].force = &force_hmc;
  integrator[currentlevel].force_par = malloc(sizeof(force_hmc_par));
  ((force_hmc_par*)(integrator[currentlevel].force_par))->n_pf = _update_par.nf/2;
  ((force_hmc_par*)(integrator[currentlevel].force_par))->pf = pf;
  ((force_hmc_par*)(integrator[currentlevel].force_par))->mass = _update_par.mass;
  ((force_hmc_par*)(integrator[currentlevel].force_par))->hasenbusch = (_update_par.hasenbusch?2:0);
  #ifdef UPDATE_EO
  ((force_hmc_par*)(integrator[currentlevel].force_par))->b = (4.+_update_par.mass+_update_par.hasen_dm)*(4.+_update_par.mass+_update_par.hasen_dm)-(4.+_update_par.mass)*(4.+_update_par.mass);
  #else
  ((force_hmc_par*)(integrator[currentlevel].force_par))->b = _update_par.hasen_dm;
  #endif
  ((force_hmc_par*)(integrator[currentlevel].force_par))->inv_err2 = _update_par.n_force_prec;
  ((force_hmc_par*)(integrator[currentlevel].force_par))->inv_err2_flt = _update_par.n_force_prec_flt;
  integrator[currentlevel].integrator = &O2MN_multistep;
  integrator[currentlevel].next = &integrator[currentlevel+1];
  
  if(_update_par.hasenbusch) {
    currentlevel++;
    integrator[currentlevel].level = currentlevel;
    integrator[currentlevel].nsteps = _update_par.hsteps;
    integrator[currentlevel].force = &force_hmc;
    integrator[currentlevel].force_par = malloc(sizeof(force_hmc_par));
    ((force_hmc_par*)(integrator[currentlevel].force_par))->n_pf = _update_par.nf/2;
    ((force_hmc_par*)(integrator[currentlevel].force_par))->pf = pf+_update_par.nf/2;
    ((force_hmc_par*)(integrator[currentlevel].force_par))->mass = _update_par.mass+_update_par.hasen_dm;
    ((force_hmc_par*)(integrator[currentlevel].force_par))->hasenbusch = 0;
    ((force_hmc_par*)(integrator[currentlevel].force_par))->b = 0.;
    ((force_hmc_par*)(integrator[currentlevel].force_par))->inv_err2 = _update_par.h_force_prec;
    ((force_hmc_par*)(integrator[currentlevel].force_par))->inv_err2_flt = _update_par.h_force_prec_flt;
    integrator[currentlevel].integrator = &O2MN_multistep;
    integrator[currentlevel].next = &integrator[currentlevel+1];
  }
  
  currentlevel++;
  integrator[currentlevel].level = currentlevel;
  integrator[currentlevel].nsteps = _update_par.gsteps;
  integrator[currentlevel].force = &force0;
  integrator[currentlevel].force_par = (void*)malloc(sizeof(double));
  *((double*)integrator[currentlevel].force_par) = _update_par.beta;
  integrator[currentlevel].integrator = &O2MN_multistep;
  integrator[currentlevel].next = &integrator[currentlevel+1];
  
  currentlevel++;
  integrator[currentlevel].level = currentlevel;
  integrator[currentlevel].nsteps = 1;
  integrator[currentlevel].force = NULL;
  integrator[currentlevel].force_par = NULL;
  integrator[currentlevel].integrator = &gauge_integrator;
  integrator[currentlevel].next = NULL;
  
  
  init_force_hmc();
  
  
  hmc_action_par.beta = _update_par.beta;
  hmc_action_par.n_pf = n_pf;
  #ifndef ROTATED_SF
  hmc_action_par.SF_ct = _update_par.SF_ct;
  #endif
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
  if(pf!=NULL) free_spinor_field(pf); pf=NULL;
  
  if(la!=NULL) free_sfield(la); la=NULL;
  
  int i=0;
  while(1) {
    if(integrator[i].force_par != NULL) free(integrator[i].force_par);
    if(integrator[i].next == NULL) break;
    i++;
  }
  free(integrator);
  
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
  for (i=0;i<n_pf;++i)
    gaussian_spinor_field(&pf[i]);
  
  
  /* compute starting action */
  lprintf("HMC",30,"Computing action density...\n");
  local_hmc_action(NEW, &hmc_action_par, la, momenta, pf, pf);
  
  /* compute H2^{1/2}*pf = H*pf */
  lprintf("HMC",30,"Correcting pseudofermions distribution...\n");
  if(!_update_par.hasenbusch) {
    for (i=0;i<n_pf;++i) {
      spinor_field_copy_f(&pf[n_pf],&pf[i]);
      static_mass=_update_par.mass;
      D(&pf[i], &pf[n_pf]);
      spinor_field_g5_assign_f(&pf[i]);
    }
  } else {
    g5QMR_fltacc_par mpar;
    for (i=0;i<n_pf/2;++i) {
      /* S = | (a D +b) D^{-1} g5 psi |^2 */
      /* (a D +b) D^{-1} g5 psi = A */
      /* psi = g5 D (a D + b)^{-1} A */
      mpar.err2 = _update_par.n_MT_prec;
      mpar.max_iter = 0;
      mpar.err2_flt = _update_par.n_MT_prec_flt;
      mpar.max_iter_flt = 0;
      spinor_field_zero_f(&pf[n_pf]);
      static_mass=_update_par.mass+_update_par.hasen_dm;
      g5QMR_fltacc(&mpar, &D, &D_flt, &pf[i], &pf[n_pf]);
      static_mass=_update_par.mass;
      D(&pf[i], &pf[n_pf]);
      spinor_field_g5_assign_f(&pf[i]);
    }
    for (i=n_pf/2;i<n_pf;++i) {
      /* S = | Dt^{-1} g5 psi |^2 */
      /* Dt^{-1} g5 psi = A */
      /* psi = g5 Dt A */
      spinor_field_copy_f(&pf[n_pf],&pf[i]);
      static_mass=_update_par.mass+_update_par.hasen_dm;
      D(&pf[i], &pf[n_pf]);
      spinor_field_g5_assign_f(&pf[i]);
    }
  }
  
  /* integrate molecular dynamics */
  lprintf("HMC",30,"MD integration...\n");
  (*(integrator[0].integrator))(momenta,_update_par.tlen,&integrator[0]);
  
  
  /* project gauge field */
  project_gauge_field();
  represent_gauge_field();
  
  lprintf("HMC",30,"Computing new action density...\n");
  if(!_update_par.hasenbusch) {
    /* compute H2^{-1/2}*pf or H2^{-1}*pf */
    /* here we choose the first strategy which is more symmetric */
    /* for the HMC H2^-1/2 = H^-1 and we use MINRES for this inversion */
    g5QMR_fltacc_par mpar;
    for (i=0;i<n_pf;++i) {
      /* S = | D^{-1} g5 psi |^2 */
      mpar.err2 = _update_par.n_MT_prec;
      mpar.max_iter = 0;
      mpar.err2_flt = _update_par.n_MT_prec_flt;
      mpar.max_iter_flt = 0;
      spinor_field_g5_f(&pf[n_pf],&pf[i]);
      spinor_field_zero_f(&pf[i]);
      static_mass=_update_par.mass;
      g5QMR_fltacc(&mpar, &D, &D_flt, &pf[n_pf], &pf[i]);
    }
  } else {
    g5QMR_fltacc_par mpar;
    for (i=0;i<n_pf/2;++i) {
      /* S = | (D+b) D^{-1} g5 psi |^2 */
      mpar.err2 = _update_par.n_MT_prec;
      mpar.max_iter = 0;
      mpar.err2_flt = _update_par.n_MT_prec_flt;
      mpar.max_iter_flt = 0;
      spinor_field_g5_assign_f(&pf[i]);
      spinor_field_zero_f(&pf[n_pf]);
      static_mass=_update_par.mass;
      g5QMR_fltacc(&mpar, &D, &D_flt, &pf[i], &pf[n_pf]);
      static_mass=_update_par.mass+_update_par.hasen_dm;
      D(&pf[i], &pf[n_pf]);
    }
    for (i=n_pf/2;i<n_pf;++i) {
      /* S = | Dt^{-1} g5 psi |^2 */
      mpar.err2 = _update_par.h_MT_prec;
      mpar.max_iter = 0;
      mpar.err2_flt = _update_par.h_MT_prec_flt;
      mpar.max_iter_flt = 0;
      spinor_field_g5_f(&pf[n_pf],&pf[i]);
      spinor_field_zero_f(&pf[i]);
      static_mass=_update_par.mass+_update_par.hasen_dm;
      g5QMR_fltacc(&mpar, &D, &D_flt, &pf[n_pf], &pf[i]);
    }
  }
  
  /* compute new action */
  local_hmc_action(DELTA, &hmc_action_par, la, momenta, pf, pf);
  
  /* Metropolis test */
  deltaH=0.;
  _MASTER_FOR(la->type,i) {
    deltaH+=*_FIELD_AT(la,i);
  }
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


