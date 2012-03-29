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
#include "rational_functions.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"

/* State quantities for RHMC */
suNg_av_field *momenta=NULL;
spinor_field *pf=NULL;
rhmc_par _update_par={0};
rational_app r_S={0};  /* used for computing the action S in the metropolis test */
rational_app r_MD={0}; /* used in the action MD evolution */
rational_app r_HB={0};  /* used in pseudofermions heatbath */
integrator_par *integrator = NULL;
double minev, maxev; /* min and max eigenvalue of H^2 */
/* END of State */


static short int init=0;


static suNg_field *u_gauge_old=NULL;
static scalar_field *la=NULL; /* local action field for Metropolis test */

/* this is the basic operator used in the update */
void H2(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
    g5Dphi_eopre_sq(_update_par.mass, out, in);
#else
    g5Dphi_sq(_update_par.mass, out, in);
#endif
}

void H(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
    g5Dphi_eopre(_update_par.mass, out, in);
#else
    g5Dphi(_update_par.mass, out, in);
#endif
}

void H_flt(spinor_field_flt *out, spinor_field_flt *in){
#ifdef UPDATE_EO
    g5Dphi_eopre_flt((float)(_update_par.mass), out, in);
#else
    g5Dphi_flt((float)(_update_par.mass), out, in);
#endif
}

static int gcd(int a, int b) {
    while (b!=0){ int t=b; b=a%t; a=t; }
    return a;
}

static void reduce_fraction(int *a, int *b){
    int f=gcd(abs(*a),abs(*b));
    if (*b!=0 && f!=1){ *a/=f; *b/=f; }
}

void init_rhmc(rhmc_par *par){
  
	if (init) {
	  /* already initialized */
	  lprintf("RHMC",0,"WARNING: RHMC already initialized!\nWARNNG: Ignoring call to init_hmc.\n");
		return;
	}
	
  lprintf("RHMC",0,"Initializing...\n");
  
  /* copy input parameters into the internal variable and make some tests */
  _update_par=*par;
  /* no tests yet... */
  
  lprintf("RHMC",10,
	  "Number of Flavors = %d\n"
	  "beta = %.8f\n"
	  "Mass = %.8f\n"
	  #ifdef ROTATED_SF
	  "rotatedSF zf = %.8f\n"
	  "rotatedSF ds = %.8f\n"
	  "rotatedSF sign = %d\n"
	  #endif /* ROTATED_SF */
	  "Metropolis test precision = %.8e\n"
	  "MD precision = %.8e\n"
	  "PF heat-bath precision = %.8e\n"
	  "RHMC force precision = %.8e\n"
	  "Number of pseudofermions = %d\n"
	  "MD trajectory length = %.8f\n"
	  "MD steps = %d\n"
	  "MD gauge substeps = %d\n"
	  ,_update_par.nf
	  ,_update_par.beta
	  ,_update_par.mass
	  #ifdef ROTATED_SF
	  ,_update_par.SF_zf
	  ,_update_par.SF_ds
	  ,_update_par.SF_sign
	  #endif /* ROTATED_SF */
	  ,_update_par.MT_prec
	  ,_update_par.MD_prec
	  ,_update_par.HB_prec
	  ,_update_par.force_prec
	  ,_update_par.n_pf
	  ,_update_par.tlen
	  ,_update_par.nsteps
	  ,_update_par.gsteps
	  );
	
	/* allocate space for the backup copy of gfield */
	if(u_gauge_old==NULL) u_gauge_old=alloc_gfield(&glattice);
	suNg_field_copy(u_gauge_old,u_gauge);
	
	/* allocate momenta */
	if(momenta==NULL) momenta = alloc_avfield(&glattice);
	
	/* allocate pseudofermions */
	if(pf==NULL) {
	  pf=alloc_spinor_field_f(_update_par.n_pf,
	    #ifdef UPDATE_EO
	    &glat_even /* even lattice for preconditioned dynamics */
	    #else
	    &glattice /* global lattice */
	    #endif 
	    );
	}
	
	/* allocate memory for the local action */
	if(la==NULL) la=alloc_sfield(&glattice);
	
	/* represent gauge field and find min and max eigenvalue of H^2 */
	represent_gauge_field();
	find_spec_H2(&maxev,&minev, par->mass); /* find spectral interval of H^2 */
	
	/* set up rational approx needed for RHMC */
	/* r_S = x^{-Nf/(4*NPf)} is used in the metropolis test */
	r_S.order=1;
	r_S.n=-_update_par.nf;
	r_S.d=4*_update_par.n_pf; 
	reduce_fraction(&r_S.n,&r_S.d);
	r_S.rel_error=_update_par.MT_prec;
	r_app_alloc(&r_S);
	r_app_set(&r_S,minev,maxev);
	/* r_D = x^{-Nf/(2*NPf)} is used in the molecula dynamics */
	r_MD.order=1;
	r_MD.n=-_update_par.nf;
	r_MD.d=2*_update_par.n_pf; 
	reduce_fraction(&r_MD.n,&r_MD.d);
	r_MD.rel_error=_update_par.MD_prec;
	r_app_alloc(&r_MD);
	r_app_set(&r_MD,minev,maxev);
	/* r_HB = x^{+Nf/(4*NPf)} is used in the heat bath for pseudofermions */
	r_HB.order=1;
	r_HB.n=_update_par.nf;
	r_HB.d=4*_update_par.n_pf; 
	reduce_fraction(&r_HB.n,&r_HB.d);
	r_HB.rel_error=_update_par.HB_prec;
	r_app_alloc(&r_HB);
	r_app_set(&r_HB,minev,maxev);



  /* integrator */
  integrator = (integrator_par*)malloc(sizeof(integrator_par)*3);

  integrator[0].level = 0;
  integrator[0].tlen = _update_par.tlen;
  integrator[0].nsteps = _update_par.nsteps;
  integrator[0].force = &force_rhmc;
  integrator[0].force_par = malloc(sizeof(force_rhmc_par));
  ((force_rhmc_par*)(integrator[0].force_par))->n_pf = _update_par.n_pf;
  ((force_rhmc_par*)(integrator[0].force_par))->pf = pf;
  ((force_rhmc_par*)(integrator[0].force_par))->ratio = &r_MD;
  ((force_rhmc_par*)(integrator[0].force_par))->inv_err2 = _update_par.force_prec;
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
  
  init_force_rhmc(r_MD.order+1);
	
	init = 1;
	
	lprintf("RHMC",0,"Initialization done.\n");
	
}

void free_rhmc(){
  
  if (!init) {
    /* not initialized */
    lprintf("RHMC",0,"WARNING: RHMC not initialized!\nWARNNG: Ignoring call to free_rhmc.\n");
		return;
	}
	
	/* free momenta */
	if(u_gauge_old!=NULL) free_gfield(u_gauge_old); u_gauge_old=NULL;
	if(momenta!=NULL) free_avfield(momenta); momenta=NULL;
	if(pf!=NULL) free_spinor_field(pf); pf=NULL;
	
	if(la!=NULL) free_sfield(la); la=NULL;
	
	r_app_free(&r_S);
	r_app_free(&r_MD);
	r_app_free(&r_HB);
  
  int i=0;
  while(1) {
    if(integrator[i].force_par != NULL) free(integrator[i].force_par);
    if(integrator[i].next == NULL) break;
    i++;
  }
  free(integrator);
  
  free_force_rhmc();
	
	init = 0;
	
	lprintf("RHMC",0,"Memory deallocated.\n");
	
}

int update_rhmc(){
    
    double deltaH;
    double oldmax,oldmin;
    _DECLARE_INT_ITERATOR(i);
    
    if (!init) {
        /* not initialized */
        lprintf("RHMC",0,"WARNING: RHMC not initialized!\nWARNNG: Ignoring call to update_rhmc.\n");
		return -1;
    }
    
    /* generate new momenta and pseudofermions */
    lprintf("RHMC",30,"Generating gaussian momenta and pseudofermions...\n");
    gaussian_momenta(momenta);
    for (i=0;i<_update_par.n_pf;++i)
        gaussian_spinor_field(&pf[i]);

#if defined(BASIC_SF) && !defined(NDEBUG)
    lprintf("MAIN",0,"SF_test_force_bcs(START): %1.8e\n",SF_test_force_bcs(momenta));
    for (i=0;i<_update_par.n_pf;++i)
    	lprintf("MAIN",0,"SF_test_spinor_bcs(START,%d): %1.8e\n",i,SF_test_spinor_bcs(&pf[i]));
#endif /* BASIC_SF && ! NDEBUG*/
      
    /* compute starting action */
    lprintf("RHMC",30,"Computing action density...\n");
    local_hmc_action(NEW, la, momenta, pf, pf);
    
    /* compute H2^{a/2}*pf */
    lprintf("RHMC",30,"Correcting pseudofermions distribution...\n");
    for (i=0;i<_update_par.n_pf;++i)
        rational_func(&r_HB, &H2, &pf[i], &pf[i]);
    
    /* integrate molecular dynamics */
    lprintf("RHMC",30,"MD integration...\n");
    (*(integrator[0].integrator))(momenta,&integrator[0]);
    
    /* project gauge field */
    project_gauge_field();
    represent_gauge_field();
    
    /* test min and max eigenvalue of H2 and update approx if necessary */
    /* now it just tests the approx !!! */
    oldmax = maxev; /* save old max */
    oldmin = minev; /* save old min */
    find_spec_H2(&maxev,&minev, _update_par.mass); /* find spectral interval of H^2 */
    r_app_set(&r_S,minev,maxev);
    r_app_set(&r_MD,minev,maxev);
    r_app_set(&r_HB,minev,maxev);
    
    lprintf("RHMC",30,"Computing new action density...\n");
    /* compute H2^{-a/2}*pf or H2^{-a}*pf */
    /* here we choose the first strategy which is more symmetric */
    for (i=0;i<_update_par.n_pf;++i)
        rational_func(&r_S, &H2, &pf[i], &pf[i]);
    
    /* compute new action */
    local_hmc_action(DELTA, la, momenta, pf, pf);
 
#if defined(BASIC_SF) && ! defined(NDEBUG)
    lprintf("MAIN",0,"SF_test_force_bcs(END): %1.8e\n",SF_test_force_bcs(momenta));
    for (i=0;i<_update_par.n_pf;++i)
    	lprintf("MAIN",0,"SF_test_spinor_bcs(END,%d): %1.8e\n",i,SF_test_spinor_bcs(&pf[i]));
#endif /* BASIC_SF && ! NDEBUG*/
    
    
    /* Metropolis test */
    deltaH=0.;
    _MASTER_FOR(la->type,i) {
        deltaH+=*_FIELD_AT(la,i);
    }
    global_sum(&deltaH, 1);
    lprintf("RHMC",10,"[DeltaS = %1.8e][exp(-DS) = %1.8e]\n",deltaH,exp(-deltaH));
    
    if(deltaH<0.) {
        suNg_field_copy(u_gauge_old,u_gauge);
    } else {
        double r;
        if (PID==0) { ranlxd(&r,1); if(r<exp(-deltaH)) r=1.0; else r=-1.0;}  /* make test on on PID 0 */
        bcast(&r,1);
        if(r>0.) {
            suNg_field_copy(u_gauge_old,u_gauge);
        } else {
            lprintf("RHMC",10,"Configuration rejected.\n");
            suNg_field_copy(u_gauge,u_gauge_old);
            start_gf_sendrecv(u_gauge); /* this may not be needed if we always guarantee that we copy also the buffers */
            represent_gauge_field();
            
            /* revert the approx to the old one */
            maxev=oldmax;
            minev=oldmin;
            r_app_set(&r_S,minev,maxev);
            r_app_set(&r_MD,minev,maxev);
            r_app_set(&r_HB,minev,maxev);
            
            return 0;
        }
    }
    
    lprintf("RHMC",10,"Configuration accepted.\n");
    
    return 1;
}



int update_rhmc_o(){
    
    double deltaH;
    double oldmax,oldmin;
    _DECLARE_INT_ITERATOR(i);
    
    if(!init)
        return -1;
    
    /* generate new momenta and pseudofermions */
    static unsigned int calln=0;
    
    if((calln++&1)==0){
        lprintf("RHMC",30,"Generating gaussian momenta and pseudofermions...\n");
        gaussian_momenta(momenta);
        for (i=0;i<_update_par.n_pf;++i)
            gaussian_spinor_field(&pf[i]);
    }
    else
        lprintf("RHMC",30,"NOT Generating momenta and pseudofermions...\n");
    /* compute starting action */
    lprintf("RHMC",30,"Computing action density...\n");
    local_hmc_action(NEW, la, momenta, pf, pf);
    
    /* compute H2^{a/2}*pf */
    lprintf("RHMC",30,"Correcting pseudofermions distribution...\n");
    for (i=0;i<_update_par.n_pf;++i)
        rational_func(&r_HB, &H2, &pf[i], &pf[i]);
    
    /* integrate molecular dynamics */
    lprintf("RHMC",30,"MD integration...\n");
    (*(integrator[0].integrator))(momenta,&integrator[0]);
    
    /* project gauge field */
    project_gauge_field();
    represent_gauge_field();
    
    /* test min and max eigenvalue of H2 and update approx if necessary */
    /* now it just tests the approx !!! */
    oldmax = maxev; /* save old max */
    oldmin = minev; /* save old min */
    find_spec_H2(&maxev,&minev, _update_par.mass); /* find spectral interval of H^2 */
    r_app_set(&r_S,minev,maxev);
    r_app_set(&r_MD,minev,maxev);
    r_app_set(&r_HB,minev,maxev);
    
    lprintf("RHMC",30,"Computing new action density...\n");
    /* compute H2^{-a/2}*pf or H2^{-a}*pf */
    /* here we choose the first strategy which is more symmetric */
    for (i=0;i<_update_par.n_pf;++i)
        rational_func(&r_S, &H2, &pf[i], &pf[i]);
    
    /* compute new action */
    local_hmc_action(DELTA, la, momenta, pf, pf);
    
    /* Metropolis test */
    deltaH=0.;
    _MASTER_FOR(la->type,i) {
        deltaH+=*_FIELD_AT(la,i);
    }
    global_sum(&deltaH, 1);
    lprintf("RHMC",10,"[DeltaS = %1.8e][exp(-DS) = %1.8e]\n",deltaH,exp(-deltaH));
    
    suNg_field_copy(u_gauge_old,u_gauge);
    
    return 1;
}





