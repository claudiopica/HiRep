/******************************************************************************
*
* File trunc_hairpin.c
*
* This program computes the hairpin term of the isosinglet mesonic correlators.
* It supports:
* - Calculation of low lying eigenvalues.
* - Dilution ("Dublin method", see Foley, Juge, O'Cais, Peardon, Ryan,
*              Skullerud, hep-lat/0505023).
* - Hopping parameter expansion.
* - Truncated solver method (see Sara Collins, Gunnar Bali, Andreas Schafer,
*                            hep-lat/07093217).
* The inverse of the dirac operator is split into the contribution of the
* low-lying eigenvalues and the rest. The former part is exactly computed.
* The latter can be expanded using the hopping expansion (perturbative terms
* are computed exactly), with the final remainder term estimated stochastically.
* The techniques of dilution and truncation are implemeted, in order to
* reduce the stochastic fluctuations.
*
* Parameters.
* n_eigenvalues : number of eigenvalues needed to compute the firts part
*                 of the quark propagator
* nevt :          nevt parameter of eva
*                 to estimate the second part of the quark propagator
* omega1 :        the absolute precision required to compute the eigenvalues
* omega2 :        the relative precision required to compute the eigenvalues
* kmax :          maximal degree of the Chebyshev polynomials used (for ev)
* imax :          maximal number of subspace iterations (for calculating ev)
*
* dilution :      dilution scheme used (0=none, 1=time, 2=time/spin, 3=exact)
*
* hopping_order : order to which hopping expansion is computed (-1 = disabled)
*
* n_truncation_steps :   number of iterations after which to stop inverter
* n_sources_truncation : number of sources for the truncated estimate
* n_sources_correction : number of sources for the correction term
* inverter_precision :   the precision to be used in the inversion routines
*
*
* Authors: Agostino Patella
*          Gregory Moraitis
*
******************************************************************************/


/******************************************************************************
* 
* H^{-1} = 
*
*     M    K^{2k} g5 e_i e_i^+
*   \sum ----------------------- +
*    k=0       (4+m)^{2k+1}
*
*     Nev      K^{2k+1}
*  + \sum  -------------- v_i [D^{-1} g5 v_i]^+ +
*     i=1   (4+m)^{2k+1}
*
*     1  dr*Nt     K^{2k+1}
*  + --- \sum  -------------- P[Nev] I(n,D,a g5 ł_i)/a ł_i^+ +
*     Nt  i=1   (4+m)^{2k+1}
*
*     1  dr*Nc     K^{2k+1}
*  + --- \sum  -------------- P[Nev] { D^{-1} g5 ß_i - I(n,D,a g5 ß_i)/a } ß_i^+
*     Nc  i=1   (4+m)^{2k+1}
*
*
******************************************************************************/

#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "error.h"
#include "logger.h"
#include "memory.h"
#include "random.h"
#include "update.h"
#include "communications.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>


#define QMR_INVERTER
#define EO_PRE

static int loglevel = 0;
static int init_flag = 0;

static ata_qprop_pars pars;
static double shift[256];
static int n_dilution_slices;

/*
static const float omega1 = 1e-16;
static const float omega2 = 1e-10;
static const float acc = 1e-16;
*/

static spinor_field **ev;       /* [m][a], m<n_masses, a<n_eigenvalues */

static spinor_field* max_H2_ev_ws;
static int max_H2_ev(double *max);

static spinor_field* compute_evs_ws;
static void compute_evs();

static void project_on_higher_evs(spinor_field *sp, int m);

static spinor_field* ev_propagator_ws;
static void ev_propagator(complex** prop);

static spinor_field* stoc_propagator_ws;
static void stoc_propagator(complex** prop, int n_src, int mode);

static spinor_field* hopping_propagator_ws;
static void hopping_propagator(complex** prop);

static spinor_field* hopping_remainder_ws;
static void hopping_remainder(spinor_field *out, spinor_field *in, int m);

static void create_diluted_source(spinor_field *source, int di, int dilution);

static spinor_field* create_sinks_QMR_ws;
enum {INVERSION, TRUNCATION, CORRECTION};
static void create_sinks_QMR(spinor_field *source, spinor_field *sink, int mode);

static void add_source_sink_contraction(complex *out, spinor_field *source, spinor_field *sink, double z);


#ifdef QMR_INVERTER
static spinor_field *QMR2_source;
static spinor_field *QMR2_sinks, *QMR2_sinks_trunc;
#endif /* QMR_INVERTER */


static void locH2(spinor_field *out, spinor_field *in);

#ifdef QMR_INVERTER
static void QMR_init();
static void D_qmr(spinor_field *out, spinor_field *in);
#ifdef EO_PRE
static void D_qmr_eo(spinor_field *out, spinor_field *in);
static void D_qmr_oe(spinor_field *out, spinor_field *in);
#endif /* EO_PRE */
#endif /* QMR_INVERTER */



/*******************************************************************************
* prop[xp][m][t*16+SPIN_2D_INDEX(a,b)] = 
* 
*    1   
*   ---  sum    tr_c [D^{-1} g5]_{a,b}(t,x,y,z;t,x,y,z)
*    V3  x,y,z
*
* tr_c = trace over the color indices
* xp = 0, ... , n_points (n_points independent stochastic estimates)
*******************************************************************************/

void traced_ata_qprop(complex*** prop, int n_points) {
/*
	prop[x][m][t*16+i]
	x in [0, n_points-1]
	m in [0, n_masses-1]
	t in [0 , T-1]
	i in [0, 15]
*/

	int m, xp;

#ifdef QMR_INVERTER
        if(pars.n_sources_truncation > 0 || pars.dilution == EXACT) QMR_init();
#endif /* QMR_INVERTER */

	for(xp = 0; xp < n_points; xp++)
	for(m = 0; m < pars.n_masses; m++)
	  memset(prop[xp][m], '\0', sizeof(complex)*GLB_T*16);

	/* If performing exact calculation (dilution==3) just invert fully
           using point sources. Otherwise, use ev, hopexp, truncation, etc. */
	if (pars.dilution == EXACT) {
	  compute_evs();
	  ev_propagator(prop[0]);
	  if(pars.hopping_order>=0) hopping_propagator(prop[0]);
	  stoc_propagator(prop[0],1,INVERSION);
	  for(xp = 1; xp < n_points; xp++)
	    for(m = 0; m < pars.n_masses; m++)
	      memcpy(prop[xp][m], prop[0][m], sizeof(complex)*GLB_T*16);

	} else {

	  compute_evs();
	  ev_propagator(prop[0]);

	  if(pars.hopping_order>=0) hopping_propagator(prop[0]);
	  
	  for(xp = 1; xp < n_points; xp++)
	    for(m = 0; m < pars.n_masses; m++) {
	      memcpy(prop[xp][m], prop[0][m], sizeof(complex)*GLB_T*16);
	    }
	
	for(xp = 0; xp < n_points; xp++) {
	  if(pars.n_truncation_steps>0) {
	    stoc_propagator(prop[xp],pars.n_sources_truncation,TRUNCATION);
	    stoc_propagator(prop[xp],pars.n_sources_correction,CORRECTION);
	  } else
	    stoc_propagator(prop[xp],pars.n_sources_truncation,INVERSION);
	    }
	}
}


static double hmass=0.;
static void locH2(spinor_field *out, spinor_field *in){
  g5Dphi_sq(hmass, out, in);
}
#ifdef QMR_INVERTER
static void D_qmr(spinor_field *out, spinor_field *in){
	Dphi(hmass,out,in);
}
#ifdef EO_PRE
static void D_qmr_eo(spinor_field *out, spinor_field *in){
  Dphi_eopre(hmass,out,in);
}
static void D_qmr_oe(spinor_field *out, spinor_field *in){
  Dphi_oepre(hmass,out,in);
}
#endif /* EO_PRE */
#endif /* QMR_INVERTER */


/******************************************************************************
*
* required workspace: 3
*
******************************************************************************/
static int max_H2_ev(double *max) {
  double norm, oldmax, dt;
  spinor_field *s1,*s2,*s3;
  int count;

  s1=max_H2_ev_ws;
  s2=s1+1;
  s3=s2+1;

  gaussian_spinor_field(s1);
  norm=sqrt(spinor_field_sqnorm_f(s1));
  spinor_field_mul_f(s1,1./norm,s1);
  norm=1.;

  dt=1.;

  locH2(s3,s1);

  count=1;
  do {
    ++count;
    spinor_field_mul_f(s1,dt,s3);
    norm=sqrt(spinor_field_sqnorm_f(s1));
    spinor_field_mul_f(s1,1./norm,s1);

    oldmax=*max;
    locH2(s3,s1);
    *max=spinor_field_prod_re_f(s1,s3);
  } while (fabs((*max-oldmax)/(*max))>1.e-3);

  *max*=1.1; /* do not know exact bound */
  
  return count;
}



/******************************************************************************
*
* required workspace: nevt-n_eigenvalues+2
*
******************************************************************************/

/* NOT WORKING FOR n_eigenvalues=1, but OK otherwise. */

static void compute_evs() {
  int m, p;
  spinor_field *ev_mask, *eva_ws;
  double *d;
  int status;
  double ubnd=0.;

#ifndef NDEBUG
  spinor_field *test = alloc_spinor_field_f(1,&glattice);
#endif /* NDEBUG */

  if(pars.n_eigenvalues<=0) return;
	
  ev_mask=malloc(sizeof(spinor_field)*pars.eva_nevt);
  for(p=pars.n_eigenvalues; p<pars.eva_nevt; p++)
    ev_mask[p]=compute_evs_ws[p-pars.n_eigenvalues];
  eva_ws=compute_evs_ws+pars.eva_nevt-pars.n_eigenvalues;
  d=malloc(sizeof(double)*pars.eva_nevt);
  
  for(m=0; m<pars.n_masses; m++) {
    for(p=0; p<pars.n_eigenvalues; p++)
      ev_mask[p]=ev[m][p];
      
    int init=0;
    if(m>0) {
      init=2;
      for(p=0; p<pars.n_eigenvalues; p++)
	spinor_field_copy_f(ev[m]+p,ev[m-1]+p);
    }
    hmass=pars.mass[m];
    max_H2_ev(&ubnd);
    eva(pars.n_eigenvalues,pars.eva_nevt,init,pars.eva_kmax,pars.eva_imax,ubnd,pars.eva_omega1,pars.eva_omega2,locH2,eva_ws,ev_mask,d,&status);

#ifndef NDEBUG
    int q;
    double dtmp[2];
    for(p=0; p<pars.n_eigenvalues; p++) {
      spinor_field_zero_f(test);
      locH2(test,&ev[m][p]);
      spinor_field_mul_add_assign_f(test, -d[p], &ev[m][p]);
      lprintf("COMPUTE_EVS",0,"Computing |H2*v_%d - lambda_%d*v_%d|^2 for mass=%f (must be about 0) = %e\n",p,p,p,pars.mass[m],spinor_field_sqnorm_f(test));
    }
    for(p=0; p<pars.n_eigenvalues; p++) {
      for(q=0; q<pars.n_eigenvalues; q++) {
	dtmp[0]=spinor_field_prod_re_f(&ev[m][p],&ev[m][q]);
	dtmp[1]=spinor_field_prod_im_f(&ev[m][p],&ev[m][q]);
	lprintf("COMPUTE_EVS",0,"Ortonormality test [%d,%d] = ( %e , %e )\n",p,q,dtmp[0],dtmp[1]);
      }
    }
#endif /* NDEBUG */

  }


#ifndef NDEBUG
	  free_spinor_field_f(test);
#endif /* NDEBUG */



  afree(ev_mask);
  afree(d);
}


/******************************************************************************
* 
* sp -> P[Nev] sp = 
*
*        Nev
*  sp - \sum  { v_i^+ sp } v_i +
*        i=1
*
******************************************************************************/
static void project_on_higher_evs(spinor_field *sp, int m) {
	int a;
	complex alpha;
	
	for(a=0; a<pars.n_eigenvalues; a++) {
		alpha.re = -spinor_field_prod_re_f(&ev[m][a],sp);
		alpha.im = -spinor_field_prod_im_f(&ev[m][a],sp);
		spinor_field_mulc_add_assign_f(sp,alpha,&ev[m][a]);
	}
}


/******************************************************************************
* 
* prop[m] += 
*
*   Nev      K^{2k+1}
*  \sum  -------------- v_i [D^{-1} g5 v_i]^+ +
*   i=1   (4+m)^{2k+1}
*
* required workspace: 3
*
******************************************************************************/
static void ev_propagator(complex** prop) {
  mshift_par QMR_par;
  int m,p;
  int cgiter=0;
  double shift[1]={0.};
  spinor_field *source, *sink, *tmp=NULL;
  
  source=ev_propagator_ws;
  sink=source+1;
  tmp=sink+1;
  
  /* set up inverters parameters */
  QMR_par.n = 1;
  QMR_par.shift = shift;
  QMR_par.err2 = pars.inverter_precision;
  QMR_par.max_iter = 0;
  
  for(m = 0; m < pars.n_masses; m++)
    for(p = 0; p < pars.n_eigenvalues; p++) {
      hmass=pars.mass[m];
      spinor_field_g5_f(tmp, ev[m]+p);
      cgiter+=g5QMR_mshift(&QMR_par, &D_qmr, tmp, source);
      if(pars.hopping_order<0) {
	add_source_sink_contraction(prop[m], source, ev[m]+p, 1.0f);
      } else {
	hopping_remainder(sink, ev[m]+p, m);
    	add_source_sink_contraction(prop[m], source, sink, 1.0f);
      }
    }
}


/******************************************************************************
* 
* IF mode==TRUNCATION
*
* prop[m] += 
*
*   1  dr*N      K^{2k+1}
*  --- \sum  -------------- P[Nev] I(n,D,a g5 ß_i)/a ß_i^+
*   N   i=1   (4+m)^{2k+1}
* 
* IF mode==CORRECTION
*
* prop[m] += 
*
*   1  dr*N      K^{2k+1}
*  --- \sum  -------------- P[Nev] { D^{-1} g5 ß_i - I(n,D,a g5 ß_i)/a } ß_i^+
*   N   i=1   (4+m)^{2k+1}
* 
* IF mode==INVERSION
*
* prop[m] += 
*
*   1  dr*N      K^{2k+1}
*  --- \sum  -------------- P[Nev] D^{-1} g5 ß_i ß_i^+
*   N   i=1   (4+m)^{2k+1}
*
* WHERE ß_i are diluted sources
*
* required workspace: 2+n_masses
*
******************************************************************************/
static void stoc_propagator(complex** prop, int n_src, int mode) {
  int r, di, m;
 spinor_field *source,*sinks,*tmp;

  source=stoc_propagator_ws;
  tmp=source+1;
  sinks=tmp+1;

  for(r = 0; r < n_src; r++)
  for(di = 0; di < n_dilution_slices; di++) {
    create_diluted_source(source,di,pars.dilution);
    create_sinks_QMR(source, sinks, mode);
    
    for(m = 0; m < pars.n_masses; m++) {
      project_on_higher_evs(sinks+m,m);
      hopping_remainder(tmp,sinks+m,m);
      spinor_field_copy_f(sinks+m,tmp);
      add_source_sink_contraction(prop[m], source, sinks+m, 1.0f/n_src);
    }
  }

}


/******************************************************************************
* 
* prop[m] += 
*
*     M    K^{2k} g5 e_i e_i^+
*   \sum -----------------------
*    k=0       (4+m)^{2k+1}
*
* required workspace: 3
*
* NB: Routine correct only for global L,T >= 4
* 
******************************************************************************/
static void hopping_propagator(complex** prop) {
  if (pars.hopping_order<0) return;
  
  int i, m;
  spinor_field *source, *sink, *tmp;

  source=hopping_propagator_ws;
  tmp=source+1;
  sink=tmp+1;
  
  /* 0th order */
  for(m = 0; m < pars.n_masses; m++)
  for(i = 0; i < GLB_T; i++) {
    prop[m][16*i].re    +=  NF*(1.0/(4.0+pars.mass[m]));
    prop[m][16*i+5].re  +=  NF*(1.0/(4.0+pars.mass[m]));
    prop[m][16*i+10].re += -NF*(1.0/(4.0+pars.mass[m]));
    prop[m][16*i+15].re += -NF*(1.0/(4.0+pars.mass[m]));
  }

  /* 1st, 2nd, 3rd order are zero, start calculating if order>=3
     NB. Only correct if L, T > 2                                  */
  if (pars.hopping_order>3) {

    int di=0;

    /* Generate exact sources (not noisy!) */

    for(di = 0; di < GLB_X*GLB_Y*GLB_Z*GLB_T*sizeof(suNf_spinor)/sizeof(complex); di++) {
    
      create_diluted_source(source,di,EXACT);
      for(m = 0; m < pars.n_masses; m++) {
        spinor_field_g5_f(tmp, source);
        spinor_field_mul_f(sink,1.0/(4.0+pars.mass[m]),tmp);
	for(i=0; i<pars.hopping_order/2; i++) {

          Dphi_(tmp, sink);
          Dphi_(sink, tmp);
          spinor_field_mul_f(sink,1.0/((4.0+pars.mass[m])*(4.0+pars.mass[m])),sink);

          if (i!=0) add_source_sink_contraction(prop[m], source, sink, 1.0);
	  }
	/*
	  for(i=0; i<order; i++) {
	    
	    Dphi(-4.0, tmp, sources);
	    
	    spinor_field_mul_f(noisy_sinks+m,1.0/(4.0+mass[m]),noisy_sinks+m);
	    spinor_field_copy_f(tmp_sink+m, noisy_sinks+m);

	    spinor_field_g5_f(noisy_sources+m, noisy_sources+m);
	      
	    if ((i%2==1) && (i!=1)) add_source_sink_contraction(prop[m], noisy_sources+m, noisy_sinks+m, 1.0);
	    }*/
      }
      
    }
    
  }
}


static double ipow(double x, int e) {
  if(e==0) return 1.0;
  else if(e<0) return ipow(1./x,-e);
  return x*ipow(x,e-1);
}

/******************************************************************************
* 
*            K^{2k+1}
*  out =  -------------- in
*          (4+m)^{2k+1}
*
* required workspace: 1
*
******************************************************************************/
static void hopping_remainder(spinor_field *out, spinor_field *in, int m) {
  int i;
  spinor_field *tmp=hopping_remainder_ws;
  
  spinor_field_copy_f(out,in);
  if(pars.hopping_order<0) return;
  
	for(i = 0; i<pars.hopping_order/2+1; i++) {
	  Dphi_(tmp,out);
	  Dphi_(out,tmp);
	}
  spinor_field_mul_f(out, ipow(-1.0/(4.0+pars.mass[m]), 2*(pars.hopping_order/2+1)), out);
}


/* Not fully tested yet with dilution!=3 !=0 */
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
/******************************************************************************
* 
* IF mode==TRUNCATION
*
* sink[m] = I(n,D,a g5 source)/a
* 
* IF mode==CORRECTION
*
* sink[m] = D^{-1} g5 source - I(n,D,a g5 source)/a
* 
* IF mode==INVERSION
*
* sink[m] = D^{-1} g5 source
*
*
* required workspace: 2+n_masses
*
******************************************************************************/
static void create_sinks_QMR(spinor_field *source, spinor_field *sink, int mode) {
	mshift_par QMR_par;
	int m;
	int cgiter=0;
#ifdef EO_PRE
  int i;
#endif /* EO_PRE */

	spinor_field *sinktmp, *sink_trunc;
#ifndef NDEBUG
  spinor_field *test = alloc_spinor_field_f(1,&glattice);
#endif /* NDEBUG */

	/* set up inverters parameters */
	QMR_par.n = pars.n_masses;
	QMR_par.shift = shift;
	QMR_par.err2 = pars.inverter_precision;
	QMR_par.max_iter = (mode==TRUNCATION)?pars.n_truncation_steps:0;

	sink_trunc=create_sinks_QMR_ws;
	sinktmp=sink_trunc+pars.n_masses;

	for(m = 0; m < pars.n_masses; m++) {
	  spinor_field_zero_f(sink + m);
	  spinor_field_zero_f(sink_trunc + m);
	}

	spinor_field_add_assign_f(source, QMR2_source);
	
	spinor_field_g5_f(source,source);

#ifdef EO_PRE

	spinor_field *b, b_even, b_odd, sinktmp_even, sinktmp_odd, source_even, source_odd, *sink_even, *sink_odd, *sink_trunc_even, *sink_trunc_odd;
  double* shift_eo;

	/* Prepare even/odd masks */	
	b=sinktmp+1;
	b_even=*b;
	b_even.type=&glat_even;
	/* b_even.ptr=b->ptr+glat_even.master_shift; */
	b_odd=*b;
	b_odd.type=&glat_odd;
	b_odd.ptr=b->ptr+glat_odd.master_shift;

	sinktmp_even=*sinktmp;
	sinktmp_even.type=&glat_even;
	/* sinktmp_even.ptr=sinktmp->ptr+glat_even.master_shift; */
	sinktmp_odd=*sinktmp;
	sinktmp_odd.type=&glat_odd;
	sinktmp_odd.ptr=sinktmp->ptr+glat_odd.master_shift;


	source_even=*source;
	source_even.type=&glat_even;
	/* source_even.ptr=source->ptr+glat_even.master_shift; */
	source_odd=*source;
	source_odd.type=&glat_odd;
	source_odd.ptr=source->ptr+glat_odd.master_shift;
	
	sink_even=(spinor_field*)malloc(sizeof(spinor_field)*pars.n_masses);
	sink_odd=(spinor_field*)malloc(sizeof(spinor_field)*pars.n_masses);
	sink_trunc_even=(spinor_field*)malloc(sizeof(spinor_field)*pars.n_masses);
	sink_trunc_odd=(spinor_field*)malloc(sizeof(spinor_field)*pars.n_masses);

	for(m = 0; m < pars.n_masses; m++) {
	  *(sink_even+m)=*(sink+m);
	  (sink_even+m)->type=&glat_even;
	  /*(sink_even+m)->ptr=(sink+m)->ptr+glat_even.master_shift; */

	  *(sink_odd+m)=*(sink+m);
	  (sink_odd+m)->type=&glat_odd;
	  (sink_odd+m)->ptr=(sink+m)->ptr+glat_odd.master_shift;

	  *(sink_trunc_even+m)=*(sink_trunc+m);
	  (sink_trunc_even+m)->type=&glat_even;
	  /*(sink_trunc_even+m)->ptr=(sink_trunc+m)->ptr+glat_even.master_shift; */

	  *(sink_trunc_odd+m)=*(sink_trunc+m);
	  (sink_trunc_odd+m)->type=&glat_odd;
	  (sink_trunc_odd+m)->ptr=(sink_trunc+m)->ptr+glat_odd.master_shift;

	}

  /* Start preconditioning & inversion */

	if (pars.n_masses==1) {

	  /* Non-multishift EO preconditioning, following DeGrand & DeTar p. 174-175.
             Requires only one non-trivial inversion. */
	  
	  Dphi_(&b_even, &source_odd);
          spinor_field_mul_f(&b_even, -1.0/(4.0+pars.mass[0]), &b_even);
          spinor_field_add_assign_f(&b_even, &source_even);
	  /* NOT OK -> spinor_field_copy_f(&b_odd, &source_odd); */
	  /* SLOW?  -> spinor_field_mul_f(&b_odd, 1.0, &source_odd); */
	  spinor_field_zero_f(&b_odd); /* Is this better? */
	  spinor_field_add_assign_f(&b_odd, &source_odd);

    cgiter+=g5QMR_mshift_trunc(&QMR_par, pars.n_truncation_steps, &D_qmr_eo, &b_even, sink_trunc_even, sink_even);

    spinor_field_mul_f(sink_even,(4.0+pars.mass[0]),sink_even);
    spinor_field_mul_f(sink_odd,1.0/(4.+pars.mass[0]),sink_odd);
    Dphi_(sink_odd, sink_even);
    spinor_field_sub_assign_f(sink_odd, &b_odd);
    spinor_field_mul_f(sink_odd, -1.0/(4.0+pars.mass[0]), sink_odd);

	  if(mode==CORRECTION) {
	    spinor_field_mul_f(sink_trunc_even,(4.0+pars.mass[0]),sink_trunc_even);
	    spinor_field_mul_f(sink_trunc_odd,1.0/(4.+pars.mass[0]),sink_trunc_odd);
	    Dphi_(sink_trunc_odd, sink_trunc_even);
	    spinor_field_sub_assign_f(sink_trunc_odd, &b_odd);
	    spinor_field_mul_f(sink_trunc_odd, -1.0/(4.0+pars.mass[0]), sink_trunc_odd);
	  }
	  
	}
	
	else {

	  /* Multishift preconditioning scheme, requires two non-trivial inversions. */
	  
	  shift_eo=(double*)malloc(sizeof(double)*pars.n_masses);
	  for(i=0; i<pars.n_masses; i++)
	    shift_eo[i]=(4.+pars.mass[0])*(4.+pars.mass[0])-(4.+pars.mass[i])*(4.+pars.mass[i]);
	  QMR_par.shift=shift_eo;
	  
	  cgiter+=g5QMR_mshift_trunc(&QMR_par, pars.n_truncation_steps, &D_qmr_eo, &source_even, sink_trunc_even, sink_even);
	  cgiter+=g5QMR_mshift_trunc(&QMR_par, pars.n_truncation_steps, &D_qmr_oe, &source_odd, sink_trunc_odd, sink_odd);
	  
    free(shift_eo);

    for(m = 0; m < pars.n_masses; m++) {
      spinor_field_copy_f(b,sink+m);
      spinor_field_mul_f(sink_even+m, 4.0+pars.mass[m], &b_even);
      Dphi_(&sinktmp_even, &b_odd);
      spinor_field_sub_assign_f(sink_even+m, &sinktmp_even);
      spinor_field_mul_f(sink_odd+m, 4.0+pars.mass[m], &b_odd);
      Dphi_(&sinktmp_odd, &b_even);
      spinor_field_sub_assign_f(sink_odd+m, &sinktmp_odd);

	    if(mode==CORRECTION) {
	      spinor_field_copy_f(b,sink_trunc+m);
	      spinor_field_mul_f(sink_trunc_even+m, 4.0+pars.mass[m], &b_even);
	      Dphi_(&sinktmp_even, &b_odd);
	      spinor_field_sub_assign_f(sink_trunc_even+m, &sinktmp_even);
	      spinor_field_mul_f(sink_trunc_odd+m, 4.0+pars.mass[m], &b_odd);
	      Dphi_(&sinktmp_odd, &b_even);
	      spinor_field_sub_assign_f(sink_trunc_odd+m, &sinktmp_odd);
	    }
    }
	}

	free(sink_even);
	free(sink_odd);
	free(sink_trunc_even);
	free(sink_trunc_odd);
	    
#else /* EO_PRE */
	cgiter+=g5QMR_mshift_trunc(&QMR_par, pars.n_truncation_steps, &D_qmr, source, sink_trunc, sink);
#endif /* EO_PRE */

#ifndef NDEBUG
  if(mode!=TRUNCATION) {
	  for(m = 0; m < pars.n_masses; m++) {
      Dphi(pars.mass[m],test,sink+m);
      spinor_field_sub_assign_f(test,source);
	    lprintf("GET_SINKS_QMR",0,"Invesion test (1) for mass=%f (must be about 0) = %e\n",pars.mass[m],spinor_field_sqnorm_f(test));
	  }
  }
#endif /* NDEBUG */

	spinor_field_g5_f(source,source);

	spinor_field_sub_assign_f(source, QMR2_source);

	for(m = 0; m < pars.n_masses; m++) {

	  if (mode==CORRECTION) {
	    spinor_field_sub_assign_f(sink + m, QMR2_sinks + m);
#ifndef NDEBUG
      g5Dphi(pars.mass[m],test,sink+m);
      spinor_field_sub_assign_f(test,source);
  	  lprintf("GET_SINKS_QMR",0,"Invesion test (2) for mass=%f (must be about 0) = %e\n",pars.mass[m],spinor_field_sqnorm_f(test));
#endif /* NDEBUG */
	    spinor_field_sub_assign_f(sink_trunc + m, QMR2_sinks_trunc + m);
	    spinor_field_sub_assign_f(sink + m, sink_trunc + m); 
	  }
	  else if (mode==INVERSION) {
	    spinor_field_sub_assign_f(sink + m, QMR2_sinks + m);
#ifndef NDEBUG
      g5Dphi(pars.mass[m],test,sink+m);
      spinor_field_sub_assign_f(test,source);
  	  lprintf("GET_SINKS_QMR",0,"Invesion test (2) for mass=%f (must be about 0) = %e\n",pars.mass[m],spinor_field_sqnorm_f(test));
#endif /* NDEBUG */
	  }
	  else if (mode==TRUNCATION) {
	    spinor_field_sub_assign_f(sink + m, QMR2_sinks_trunc + m);
	  }

	}

#ifndef NDEBUG
  free_spinor_field_f(test);
#endif /* NDEBUG */
	
	lprintf("GET_SINKS_QMR",loglevel+1,"QMR MVM = %d\n",cgiter);

}
#endif /* QMR_INVERTER */


static void add_source_sink_contraction(complex *out, spinor_field *source, spinor_field *sink, double z) {
  int i, j, t, x, index;
  suNf_vector *eta, *csi;
  complex tmp;
#ifndef NDEBUG
  complex trace, g5trace;
  double norm;
#endif /* NDEBUG */

#ifndef NDEBUG
  trace.re = trace.im = 0.;
  g5trace.re = g5trace.im = 0.;
#endif /* NDEBUG */

  int point[4];
  complex out_tmp[GLB_T*16];

  for(i = 0; i < GLB_T*16; i++)
    out_tmp[i].re = out_tmp[i].im = 0.0;

  for(t = 0; t < GLB_T; t++) {
    for(x = 0; x < GLB_X*GLB_Y*GLB_Z; x++) {
      point[0]=t;
      point[1]=(x/(GLB_X*GLB_Y)) % GLB_Z;
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
	    out_tmp[t*16+SPIN_2D_INDEX(i,j)].re += tmp.re*z/(GLB_X*GLB_Y*GLB_Z);
	    out_tmp[t*16+SPIN_2D_INDEX(i,j)].im += tmp.im*z/(GLB_X*GLB_Y*GLB_Z);
#ifndef NDEBUG
            if(i==j) {
              trace.re += tmp.re;
              trace.im += tmp.im;
              g5trace.re += (i==0 || i==1) ? tmp.re : -tmp.re;
              g5trace.im += (i==0 || i==1) ? tmp.im : -tmp.im;
            }
#endif /* NDEBUG */
          }
        }
      }
    }
  }

  global_sum((double*)out_tmp,2*GLB_T*16);

  for(i = 0; i < GLB_T*16; i++) {
    out[i].re += out_tmp[i].re;
    out[i].im += out_tmp[i].im;
  }

  lprintf("ADD_SOURCE_SINK_CONTRACTION",loglevel+2,"Written in %p\n",out);

#ifndef NDEBUG
  global_sum(&trace.re,1);
  global_sum(&trace.im,1);
  global_sum(&g5trace.re,1);
  global_sum(&g5trace.im,1);

  trace.re -= spinor_field_prod_re_f(source,sink);
  trace.im -= spinor_field_prod_im_f(source,sink);
  g5trace.re -= spinor_field_g5_prod_re_f(source,sink);
  g5trace.im -= spinor_field_g5_prod_im_f(source,sink);
  norm = sqrt( trace.re*trace.re + trace.im*trace.im );
  lprintf("ADD_SOURCE_SINK_CONTRACTION",0,"Testing trace (must be about 0) = %e\n",norm);
  norm = sqrt( g5trace.re*g5trace.re + g5trace.im*g5trace.im );
  lprintf("ADD_SOURCE_SINK_CONTRACTION",0,"Testing g5trace (must be about 0) = %e\n",norm);
#endif /* NDEBUG */

}





void ata_qprop_init(ata_qprop_pars *p) {
  int m;
  
  if(init_flag != 0) return;

  /* static parameters */
  pars=*p;
  
  error(pars.n_masses<1 || pars.n_masses>256,1,"ata_qprop_init [hairpin.c]", "Bad choice for n_masses");
  hmass = p->mass[0]; /* we can put any number for the index! */
  for(m = 0; m < pars.n_masses; ++m) {
    pars.mass[m] = p->mass[m];
    shift[m] = hmass - p->mass[m];
  }
  
  error(pars.n_eigenvalues<0,1,"ata_qprop_init [hairpin.c]", "Bad choice for n_eigenvalues");
  error(pars.eva_nevt<pars.n_eigenvalues,1,"ata_qprop_init [hairpin.c]", "Bad choice for nevt");
  
  if(pars.dilution == NO_DILUTION)
    n_dilution_slices = 1;
  else if(pars.dilution == TIME_DILUTION)
    n_dilution_slices = GLB_T;
  else if(pars.dilution == TIME_SPIN_DILUTION)
    n_dilution_slices = 4*GLB_T;
  else if(pars.dilution == EXACT)
    n_dilution_slices = GLB_X*GLB_Y*GLB_Z*GLB_T*sizeof(suNf_spinor)/sizeof(complex);
  else
    error(0==0,1,"ata_qprop_init [hairpin.c]", "Bad choice for dilution");

  error(pars.n_truncation_steps<0,1,"ata_qprop_init [hairpin.c]", "Bad choice for n_truncation_steps");

  /* Use hopping_order=-1 to disable hopping expansion (hopping_order=0 is 0th order) */
  error(pars.hopping_order<-1,1,"ata_qprop_init [hairpin.c]", "Bad choice for hopping_order");

  error(pars.n_sources_truncation<0,1,"ata_qprop_init [hairpin.c]", "Bad choice for n_sources_truncation");
  error(pars.n_sources_correction<0,1,"ata_qprop_init [hairpin.c]", "Bad choice for n_sources_correction");

  lprintf("ATA_QPROP",loglevel+1,"Number of masses = %d\n",pars.n_masses);
  for(m = 0;m < pars.n_masses; m++)
    lprintf("ATA_QPROP",loglevel+1,"Mass[%d] = %f\n",m,pars.mass[m]);
    
  lprintf("ATA_QPROP",loglevel+1,"Number of eigenvalues = %d (%d)\n",pars.n_eigenvalues,pars.eva_nevt);
  lprintf("ATA_QPROP",loglevel+1,"Eva parameters (omega1,omega2, imax,kmax) = %e,%e,%d,%d\n",
          pars.eva_omega1,pars.eva_omega2,pars.eva_imax,pars.eva_kmax);
  
  lprintf("ATA_QPROP",loglevel+1,"Order of hopping expansion = %d\n",pars.hopping_order);
  
  lprintf("ATA_QPROP",loglevel+1,"Iterations after which to truncate = %d\n",pars.n_truncation_steps);
  lprintf("ATA_QPROP",loglevel+1,"Number of global noisy sources for (truncated) inversion = %d\n",pars.n_sources_truncation);
  lprintf("ATA_QPROP",loglevel+1,"Number of global noisy sources for correction to truncation = %d\n",pars.n_sources_correction);
  lprintf("ATA_QPROP",loglevel+1,"Inverter precision = %e\n",pars.inverter_precision);
  lprintf("ATA_QPROP",loglevel+1,"Number of dilution slices = %d\n",n_dilution_slices);
  lprintf("ATA_QPROP",loglevel+1,"Dilution level = %d\n",pars.dilution);

  
	/* eigenvector related stuff */
	if(pars.n_eigenvalues != 0) {
    ev = (spinor_field**)malloc(sizeof(spinor_field*)*pars.n_masses);
    ev[0] = alloc_spinor_field_f(pars.n_masses*pars.n_eigenvalues,&glattice);
    for(m = 0; m < pars.n_masses; m++)
      ev[m] = ev[0]+m*pars.n_eigenvalues;
    compute_evs_ws = alloc_spinor_field_f(pars.eva_nevt-pars.n_eigenvalues+2,&glattice);
    max_H2_ev_ws = alloc_spinor_field_f(3,&glattice);
    ev_propagator_ws = alloc_spinor_field_f(3,&glattice);
	} else {
    ev = NULL;
    compute_evs_ws = NULL;
    max_H2_ev_ws = NULL;
    ev_propagator_ws = NULL;
  }

	/* hopping expansion related stuff */
	if(pars.hopping_order >= 0) {
    hopping_propagator_ws = alloc_spinor_field_f(3,&glattice);
    hopping_remainder_ws = alloc_spinor_field_f(1,&glattice);
	} else {
    hopping_propagator_ws = NULL;
    hopping_remainder_ws = NULL;
	}

	/* noisy sources related stuff */
	if(pars.n_sources_truncation > 0 || pars.dilution == EXACT) {
	  stoc_propagator_ws = alloc_spinor_field_f(2+pars.n_masses,&glattice);
#ifdef QMR_INVERTER
	  create_sinks_QMR_ws = alloc_spinor_field_f(2+pars.n_masses,&glattice);
	  QMR2_source = alloc_spinor_field_f(1,&glattice);
	  QMR2_sinks = alloc_spinor_field_f(pars.n_masses,&glattice);
	  QMR2_sinks_trunc = alloc_spinor_field_f(pars.n_masses,&glattice);
#endif /* QMR_INVERTER */
	}else {
	  stoc_propagator_ws = NULL;
#ifdef QMR_INVERTER
	  create_sinks_QMR_ws = NULL;
	  QMR2_source = NULL;
	  QMR2_sinks = NULL;
	  QMR2_sinks_trunc = NULL;
#endif /* QMR_INVERTER */
	}

  init_flag = 1;
}


void ata_qprop_free() {
	if(init_flag != 1) return;
	
	if(pars.n_eigenvalues != 0) {
	  free_spinor_field_f(ev[0]);
    afree(ev);
    free_spinor_field_f(compute_evs_ws);
    free_spinor_field_f(max_H2_ev_ws);
    free_spinor_field_f(ev_propagator_ws);
	}
	
	if(pars.hopping_order >= 0) {
    free_spinor_field_f(hopping_propagator_ws);
    free_spinor_field_f(hopping_remainder_ws);
	}

	if(pars.n_sources_truncation > 0 || pars.dilution == EXACT) {
	  free_spinor_field_f(stoc_propagator_ws);
#ifdef QMR_INVERTER
	  free_spinor_field_f(create_sinks_QMR_ws);
	  
	  free_spinor_field_f(QMR2_source);
	  free_spinor_field_f(QMR2_sinks);
	  free_spinor_field_f(QMR2_sinks_trunc);
#endif /* QMR_INVERTER */
	}
	init_flag = 0;
}


#ifdef QMR_INVERTER
static void QMR_init() {
	mshift_par QMR_par;
	int m;
	int cgiter=0;

#ifndef NDEBUG
  spinor_field *test = alloc_spinor_field_f(1,&glattice);
#endif /* NDEBUG */
	
  gaussian_spinor_field(QMR2_source);

	/* set up inverters parameters */
	QMR_par.n = pars.n_masses;
	QMR_par.shift = shift;
	QMR_par.err2 = pars.inverter_precision;
	QMR_par.max_iter = 0;

	for(m = 0; m < pars.n_masses; m++) {
		spinor_field_zero_f(QMR2_sinks+m);
		spinor_field_zero_f(QMR2_sinks_trunc+m);
	}

	spinor_field_g5_f(QMR2_source, QMR2_source);
	cgiter+=g5QMR_mshift_trunc(&QMR_par, pars.n_truncation_steps, &D_qmr, QMR2_source, QMR2_sinks_trunc, QMR2_sinks);
	if (pars.n_truncation_steps==0) {
	  for(m = 0; m < pars.n_masses; m++) {
	    spinor_field_copy_f(QMR2_sinks_trunc+m,QMR2_sinks+m);
	  }
	}

#ifndef NDEBUG
  for(m = 0; m < pars.n_masses; m++) {
    Dphi(pars.mass[m],test,QMR2_sinks+m);
    spinor_field_sub_assign_f(test,QMR2_source);
    lprintf("QMR_INIT",0,"Invesion test (1) for mass=%f (must be about 0) = %e\n",pars.mass[m],spinor_field_sqnorm_f(test));
  }
#endif /* NDEBUG */

	spinor_field_g5_f(QMR2_source, QMR2_source);

#ifndef NDEBUG
  for(m = 0; m < pars.n_masses; m++) {
    g5Dphi(pars.mass[m],test,QMR2_sinks+m);
    spinor_field_sub_assign_f(test,QMR2_source);
    lprintf("QMR_INIT",0,"Invesion test (2) for mass=%f (must be about 0) = %e\n",pars.mass[m],spinor_field_sqnorm_f(test));
  }
#endif /* NDEBUG */

#ifndef NDEBUG
  free_spinor_field_f(test);
#endif /* NDEBUG */


	lprintf("QMR_INIT",loglevel+1,"QMR MVM = %d\n",cgiter);	
}
#endif /* QMR_INVERTER */
