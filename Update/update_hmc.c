/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
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
#include "rational_functions.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"

/* these are all alredy defined in update_rhmc.c */
/* State quantities for HMC */
/* suNg_av_field *momenta=NULL; */
/* spinor_field *pf=NULL; */
/* rhmc_par _update_par={0}; */
/* rational_app r_MD={0}; /\* used in the action MD evolution *\/ */
/* double minev, maxev; /\* min and max eigenvalue of H^2 *\/ */
/* END of State */
extern suNg_av_field *momenta;
extern spinor_field *pf;
extern rhmc_par _update_par;
extern rational_app r_MD; /* used in the action MD evolution */
extern spinor_operator H2;
extern spinor_operator H;
/* extern double minev, maxev; */ /* min and max eigenvalue of H^2 */


static short int init=0;

static suNg_field *u_gauge_old=NULL;
static scalar_field *la=NULL; /* local action field for Metropolis test */
static MINRES_par pfa;

void init_hmc(rhmc_par *par){
    
	if (init) {
        /* already initialized */
        lprintf("HMC",0,"WARNING: HMC already initialized!\nWARNNG: Ignoring call to init_hmc.\n");
		return;
    }
    
	lprintf("HMC",0,"Initializing...\n");
    
	/* copy input parameters into the internal variable and make some tests */
	_update_par=*par;
	if (_update_par.nf!=2*_update_par.n_pf) {
		lprintf("HMC",0,"The number of fermions is not twice the number of pseudofermions.\nTry with the RHMC algorithm.\n");
		error(1,1,"init_hmc","The HMC algorithm is not suitable for the parameters specified as input.\n");
	}
    
	lprintf("HMC",10,
			"Number of Flavors = %d\n"
			"beta = %.8f\n"
			"Mass = %.8f\n"
			"Metropolis test precision = %.8e\n"
			"RHMC force precision = %.8e\n"
			"Number of pseudofermions = %d\n"
			"MD trajectory length = %.8f\n"
			"MD steps = %d\n"
			"MD gauge substeps = %d\n"
			,_update_par.nf
			,_update_par.beta
			,_update_par.mass
			,_update_par.MT_prec
			,_update_par.force_prec
			,_update_par.n_pf
			,_update_par.MD_par->tlen
			,_update_par.MD_par->nsteps
			,_update_par.MD_par->gsteps
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
    
	/* set up rational approx needed for HMC */
	/* r_S = x^{-Nf/(4*NPf)} = x^-1/2 is used in the metropolis test */
	/* since H2^-1/2 = H^-1 we use the MINRES inverter */
	pfa.err2=_update_par.MT_prec;
	pfa.err2*=pfa.err2;
	pfa.max_iter=0;
	/* r_MD = x^{-Nf/(2*NPf)} = x^-1 is used in the molecular dynamics in the force*/
	r_MD.order=1;
	r_MD.n=-1;
	r_MD.d=1; 
	r_MD.rel_error=0.; /* not used */
	r_app_alloc(&r_MD);
	r_MD.a[0]=r_MD.b[0]=0.;
	r_MD.a[1]=1.;
	/* r_HB = x^{+Nf/(4*NPf)} = x^1/2 is used in the heat bath for pseudofermions */
	/* we don't need a rational approximation for this since H2^1/2 = H */
  
  /* set spinor_operator pointers H2 */
  H2.dbl=&H2_dbl;
  H2.flt=&H2_flt;
  H.dbl=&H_dbl;
  H.flt=&H_flt;
  
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
    
    r_app_free(&r_MD);
        
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
    
    /* find_spec_H2(&maxev,&minev, _update_par.mass); /\* find spectral interval of H^2 *\/ */
    
    /* generate new momenta and pseudofermions */
    lprintf("HMC",30,"Generating gaussian momenta and pseudofermions...\n");
    gaussian_momenta(momenta);
    for (i=0;i<_update_par.n_pf;++i)
        gaussian_spinor_field(&pf[i]);
    
    /* compute starting action */
    lprintf("HMC",30,"Computing action density...\n");
    local_hmc_action(NEW, la, momenta, pf, pf);
    
    /* compute H2^{1/2}*pf = H*pf */
    lprintf("HMC",30,"Correcting pseudofermions distribution...\n");
    for (i=0;i<_update_par.n_pf;++i) {
        spinor_field_copy_f(&pf[_update_par.n_pf],&pf[i]);
        H.dbl(&pf[i], &pf[_update_par.n_pf]);
    }
    
    /* integrate molecular dynamics */
    lprintf("HMC",30,"MD integration...\n");
    _update_par.integrator(momenta,_update_par.MD_par);
    
    /* project gauge field */
    project_gauge_field();
    represent_gauge_field();
    
    lprintf("HMC",30,"Computing new action density...\n");
    /* compute H2^{-1/2}*pf or H2^{-1}*pf */
    /* here we choose the first strategy which is more symmetric */
    /* for the HMC H2^-1/2 = H^-1 and we use MINRES for this inversion */
    for (i=0;i<_update_par.n_pf;++i) {
        spinor_field_copy_f(&pf[_update_par.n_pf],&pf[i]);
        MINRES(&pfa,H,&pf[_update_par.n_pf],&pf[i],0);
    }
    
    /* compute new action */
    local_hmc_action(DELTA, la, momenta, pf, pf);
    
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


