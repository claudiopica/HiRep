/***************************************************************************\
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen, Jarno Rantaharju            *
*                                                                           *
*                                                                           *
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
#include "random.h"
#include "communications.h"
#include "ranlux.h"

#define PI 3.141592653589793238462643383279502884197

//Helps QMR solver find more accurate solutions
#undef GAUSSIAN_NOISE 
//#define GAUSSIAN_NOISE 


static double hmass_pre;

static void D_pre(spinor_field *out, spinor_field *in){
#ifdef WITH_CLOVER
  Cphi_eopre(hmass_pre,out,in);
#else
  Dphi_eopre(hmass_pre,out,in);
#endif
}

static void H_pre(spinor_field *out, spinor_field *in){
  g5Dphi_eopre(hmass_pre,out,in);
}

static void H2_pre(spinor_field *out, spinor_field *in){
  g5Dphi_eopre_sq(hmass_pre,out,in);
}

//using g5 D g5, not the most efficient but not really crucial
static void Ddag_pre(spinor_field *out, spinor_field *in, spinor_field *ttmp){ 
      spinor_field_copy_f( ttmp, in ); 
      spinor_field_g5_assign_f( ttmp ); 
      H_pre(out, ttmp);
}


static int init=0;
static int init_odd=0;
static int init_eig=0;
static int neigs=0;

static mshift_par QMR_par;
static double *shift;
static double *mass;
#ifdef GAUSSIAN_NOISE
static spinor_field *QMR_noise;
static spinor_field *QMR_resdn;
#endif
static spinor_field *resd;
static spinor_field *tmp;
static spinor_field *tmp_odd;
// EVA parameters 
double *eva_val;
static spinor_field *eva_vec;
static spinor_field *tmp_sf;

enum {_g5QMR=0, _MINRES, _CG, _CG_4F};

/* Initialises the propagator, nm is the number of masses for multimass solver, 
   m is the array of masses, and acc is the inverter accuracy
*/

static void init_eva(int nevt){
  if(init_eig == 0){
    eva_val=malloc(sizeof(double)*nevt);
    eva_vec=alloc_spinor_field_f(nevt+1,&glat_even);
    tmp_sf=eva_vec+nevt;
    init_eig = 1;
  }
}

void init_propagator_eo(int nm, double *m, double acc){
  int i;
#ifdef GAUSSIAN_NOISE
  int cgiter=0;
  double norm;
#endif

  if(init==0) {
    shift=(double*)malloc(sizeof(double)*(nm));
    mass=(double*)malloc(sizeof(double)*(nm));
    hmass_pre=m[0]; /* we can put any number here!!! */
    for(i=0;i<nm;++i){
      mass[i]=m[i];
      shift[i]=(4.+hmass_pre)*(4.+hmass_pre)-(4.+m[i])*(4.+m[i]);
    }
    QMR_par.n = nm;
    QMR_par.shift = shift;
    QMR_par.err2 = .5*acc;
    QMR_par.max_iter = 0;
 
    resd=alloc_spinor_field_f(QMR_par.n,&glat_even);
    tmp=alloc_spinor_field_f(1,&glat_even);

#ifdef GAUSSIAN_NOISE
    QMR_noise=alloc_spinor_field_f(nm+1,&glat_even);
    QMR_resdn=QMR_noise+1;
    /* noisy background */
    gaussian_spinor_field(QMR_noise);
    norm=sqrt(spinor_field_sqnorm_f(QMR_noise));
    spinor_field_mul_f(QMR_noise,1./norm,QMR_noise);
#endif
    init=1;
  }
#ifdef GAUSSIAN_NOISE
  /* invert noise */
  for(i=0;i<QMR_par.n;++i) spinor_field_zero_f(&QMR_resdn[i]);
  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, QMR_noise, QMR_resdn);
  lprintf("Z2SEMWALL",10,"QMR_eo MVM = %d\n",cgiter);
#endif
}


void free_propagator_eo() {
  error(init==0,1,"calc_prop.c","propagator not initialized!");

  free_spinor_field_f(tmp);
  free_spinor_field_f(resd);


  free(shift);
  free(mass);
  
#ifdef GAUSSIAN_NOISE
  free_spinor_field_f(QMR_noise);
#endif
  if (init_odd){
    free_spinor_field_f(tmp_odd);
    init_odd = 0;
  }
  if (init_eig){
    free(eva_val);
    free_spinor_field_f(eva_vec);
    init_eig = 0;
  }
  init=0;
}




/***************************************************************************\

psi = D^{-1} eta

\***************************************************************************/


static void calc_propagator_eo_core(spinor_field *psi, spinor_field *eta, int solver) {
  spinor_field qprop_mask;
  int i, cgiter=0;
  error(init==0,1,"calc_prop.c","z2semwall method not initialized!");

  /* add source */
#ifdef GAUSSIAN_NOISE
  spinor_field_add_f(tmp,eta,QMR_noise);
#else
  spinor_field_copy_f(tmp,eta);
#endif

//if the solution vector is empty use zero guess
  if( spinor_field_sqnorm_f(psi) < 1e-28 ){
    for(i=0;i<QMR_par.n;++i){ 
      spinor_field_zero_f(&resd[i]); 
    }
  } 
  else {
    for(i=0;i<QMR_par.n;++i){ 
      psi[i].type = &glat_even; 
      spinor_field_mul_f(&resd[i],1/(4.+mass[i]),&psi[i]);
      psi[i].type = &glattice; 
    }
  }

  if(solver == _CG){
    qprop_mask.type=&glat_even;
    spinor_field_copy_f( &qprop_mask, tmp ); 
    spinor_field_g5_assign_f( &qprop_mask ); 
    H_pre(tmp, &qprop_mask);
    cgiter+=cg_mshift(&QMR_par, &H2_pre, tmp, resd);
  } else if(solver == _MINRES){
    spinor_field_g5_f(tmp,tmp);
    cgiter+=MINRES_mshift(&QMR_par, &H_pre, tmp, resd);
  } else {
    cgiter+=g5QMR_mshift(&QMR_par, &D_pre, tmp, resd);
  }
  
  for(i=0;i<QMR_par.n;++i){
#ifdef GAUSSIAN_NOISE
    spinor_field_sub_assign_f(&resd[i],&QMR_resdn[i]);
#endif
    /* compute solution */
    qprop_mask=psi[i];
    qprop_mask.type=&glat_even;
    spinor_field_mul_f(&qprop_mask,(4.+mass[i]),&resd[i]);
    qprop_mask.type=&glat_odd;
    qprop_mask.ptr=psi[i].ptr+glat_odd.master_shift;
    spinor_field_zero_f(&qprop_mask);
    Dphi_(&qprop_mask,&resd[i]);
    spinor_field_minus_f(&qprop_mask,&qprop_mask);
    if(i&1) ++cgiter; /* count only half of calls. works because the number of sources is even */
  }
  start_sf_sendrecv(psi);
  complete_sf_sendrecv(psi);
  lprintf("CALC_PROP",10,"QMR_eo MVM = %d\n",cgiter);
}








static void calc_propagator_core(spinor_field *psi, spinor_field *eta, int solver) {

  start_sf_sendrecv(eta);
  complete_sf_sendrecv(eta);

  spinor_field qprop_mask_eta;
  spinor_field qprop_mask_psi;
  int cgiter=0;
  if (init_odd==0){
    tmp_odd = alloc_spinor_field_f(1,&glat_odd);
    init_odd=1;
  }
  error(init==0,1,"calc_prop.c","calc_propagator_core method not initialized!");

  /* Construct source
     eta_even' = eta_even - D_eo D_oo^-1 eta_odd
   */
  qprop_mask_eta=*eta;
  qprop_mask_eta.type=&glat_odd;
  qprop_mask_eta.ptr=eta->ptr+glat_odd.master_shift; 
  spinor_field_mul_f(tmp_odd,(1./( 4.+ mass[0] )),&qprop_mask_eta);
  Dphi_(tmp,tmp_odd);
  qprop_mask_eta.type=&glat_even;
  qprop_mask_eta.ptr=eta->ptr;
  spinor_field_sub_f(tmp,&qprop_mask_eta,tmp);
#ifdef GAUSSIAN_NOISE
  spinor_field_add_assign_f(tmp,QMR_noise);
#endif
  //spinor_field_sub_f(resd,&qprop_mask_eta,tmp);

//if the solution vector is empty use zero guess
if( spinor_field_sqnorm_f(psi) < 1e-28 ){
  spinor_field_zero_f(resd); 
} else {
	  psi[0].type = &glat_even; 
          spinor_field_mul_f(resd,1/(4.+mass[0]),psi);
	  psi[0].type = &glattice; 
}

#ifdef GAUSSIAN_NOISE
  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, tmp, resd);
#else
  if(solver == _CG){
    qprop_mask_eta.type=&glat_even;
    spinor_field_copy_f( &qprop_mask_eta, tmp ); 
    spinor_field_g5_assign_f( &qprop_mask_eta ); 
    H_pre(tmp, &qprop_mask_eta);
    cgiter+=cg_mshift(&QMR_par, &H2_pre, tmp, resd);
  } else if(solver == _MINRES){
    spinor_field_g5_f(tmp,tmp);
    cgiter+=MINRES_mshift(&QMR_par, &H_pre, tmp, resd);
  } else {
    cgiter+=g5QMR_mshift(&QMR_par, &D_pre, tmp, resd);
  }

  //spinor_field_g5_f(tmp,tmp);
  //cgiter+=MINRES_mshift(&QMR_par, &H_pre, tmp, resd); 

  //cgiter+=g5QMR_mshift(&QMR_par, &D_pre, tmp, resd);

  //cg stuff
  /*qprop_mask_eta.type=&glat_even;
  spinor_field_copy_f( &qprop_mask_eta, tmp ); 
  spinor_field_g5_assign_f( &qprop_mask_eta ); 
  H_pre(tmp, &qprop_mask_eta);
  cgiter+=cg_mshift(&QMR_par, &H2_pre, tmp, resd);*/
#endif

#ifdef GAUSSIAN_NOISE
  spinor_field_sub_assign_f(resd,QMR_resdn);
#endif
  /* compute solution 
     psi_even = D_ee*resd_e
     psi_odd = D_oo^-1*eta_odd-D_oe resd_e
  */

  qprop_mask_psi=*psi;
  qprop_mask_psi.type=&glat_even;
  spinor_field_mul_f(&qprop_mask_psi,(4.+mass[0]),resd);

  qprop_mask_psi.type=&glat_odd;
  qprop_mask_psi.ptr=psi->ptr+glat_odd.master_shift; 
  Dphi_(&qprop_mask_psi,resd);
  
  spinor_field_sub_f(&qprop_mask_psi,tmp_odd,&qprop_mask_psi);

  ++cgiter; /* One whole call*/
  lprintf("CALC_PROP_CORE",10,"QMR_eo MVM = %d\n",cgiter);


   start_sf_sendrecv(psi);
   complete_sf_sendrecv(psi);

}

static void calc_propagator_clover(spinor_field *dptr, spinor_field *sptr)
{
	static spinor_field *etmp, *otmp, *stmp;
	static int init = 0;

	// Inverter
	mshift_par mpar;
	double tmp;

	mpar.err2 = QMR_par.err2;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &tmp;
	mpar.shift[0] = 0;

	// Allocate temporary fields
	if(init == 0)
	{
		etmp = alloc_spinor_field_f(1, &glat_even);
		otmp = alloc_spinor_field_f(1, &glat_odd);
		stmp = alloc_spinor_field_f(1, &glattice);
		init = 1;
	}

	// Destination even/odd
	spinor_field dptr_e, dptr_o;
	dptr_e = *dptr;
	dptr_e.type = &glat_even;
	dptr_o = *dptr;
	dptr_o.ptr += glat_odd.master_shift;
	dptr_o.type = &glat_odd;

	// Source even/odd
	spinor_field sptr_e, sptr_o;
	sptr_e = *stmp;
	sptr_e.type = &glat_even;
	sptr_o = *stmp;
	sptr_o.ptr += glat_odd.master_shift;
	sptr_o.type = &glat_odd;

	// Handle source
	if(sptr->type == &glat_even)
	{
		spinor_field_zero_f(stmp);
		spinor_field_copy_f(&sptr_e, sptr);
	}
	else
	{
		spinor_field_copy_f(stmp, sptr);
	}

#ifdef WITH_CLOVER

	// etmp = sptr_e - D_eo D_oo^-1 sptr_o
	Cphi_diag_inv(hmass_pre, otmp, &sptr_o);
	Dphi_(etmp, otmp);
	spinor_field_minus_f(etmp, etmp);
	spinor_field_add_assign_f(etmp, &sptr_e);

	// Call inverter
	g5QMR_mshift(&mpar, D_pre, etmp, &dptr_e);

	// dptr_o = D_oo^-1 ( sptr_o - D_oe dptr_e )
	Dphi_(&dptr_o, &dptr_e);
	spinor_field_minus_f(&dptr_o, &dptr_o);
	spinor_field_add_assign_f(&dptr_o, &sptr_o);
	Cphi_diag_inv(hmass_pre, &dptr_o, &dptr_o);

#else

	// etmp = D_ee sptr_e - D_eo sptr_o
	Dphi_(etmp, &sptr_o);
	spinor_field_minus_f(etmp, etmp);
	spinor_field_mul_add_assign_f(etmp, 4.+hmass_pre, &sptr_e);

	// Call inverter
	g5QMR_mshift(&mpar, D_pre, etmp, &dptr_e);

	// dptr_o = D_oo^-1 ( sptr_o - D_oe dptr_e )
	Dphi_(&dptr_o, &dptr_e);
	spinor_field_minus_f(&dptr_o, &dptr_o);
	spinor_field_add_assign_f(&dptr_o, &sptr_o);
	spinor_field_mul_f(&dptr_o, 1./(4.+hmass_pre), &dptr_o);

#endif
}

void calc_propagator(spinor_field *psi, spinor_field* eta, int ndilute){
	int beta,i,n_masses;
	double *m;
	m = mass;
	n_masses=QMR_par.n;
	QMR_par.n=1;
	for (beta=0;beta<ndilute;++beta){
		for (i=0;i<n_masses;++i){
			lprintf("CALC_PROPAGATOR",10,"n masses=%d, mass = %g\n",n_masses, mass[0]);
			hmass_pre = mass[0];
#ifdef WITH_CLOVER
			calc_propagator_clover(&psi[beta*n_masses+i],&eta[beta]);
#else
			calc_propagator_core(&psi[beta*n_masses+i],&eta[beta],_g5QMR);
#endif
			mass++;
		}
		mass = m;
	}
	QMR_par.n = n_masses;
	hmass_pre = mass[0];
}

void calc_propagator_eo(spinor_field *psi, spinor_field *eta, int ndilute) {
#ifdef WITH_CLOVER
	calc_propagator(psi, eta, ndilute);
#else
  int beta;
  lprintf("CALC_PROPAGATOR_EO",20,"Calculating EO propagator with ndilute: %d\n",ndilute);
  for (beta=0;beta<ndilute;++beta){
    calc_propagator_eo_core(&psi[beta*QMR_par.n],&eta[beta],_g5QMR);
  }
#endif
}

/*Different source for each mass. Needed in sequential propagators
  with multiple masses */
void calc_propagator_multisource(spinor_field *psi, spinor_field* eta, int ndilute){
  int beta,i,n_masses;
  double *m;
  m = mass;
  n_masses=QMR_par.n;
  QMR_par.n=1;
  for (i=0;i<n_masses;++i){
    hmass_pre = mass[0];
    for (beta=0;beta<ndilute;++beta){
#ifdef WITH_CLOVER
		calc_propagator_clover(&psi[beta*n_masses+i],&eta[beta*n_masses+i]);
#else
      calc_propagator_core(&psi[beta*n_masses+i],&eta[beta*n_masses+i],_g5QMR);
#endif
      mass++;
    }
  }
  QMR_par.n = n_masses;
  mass = m;
  hmass_pre = mass[0];
}

void eig_init(int nev, int nevt, int kmax, int maxiter, double lbnd, double omega1, double omega2){

 if(!init_eig){ init_eva(nevt); }
 hmass_pre = mass[0];

 double max, mupp;
 int status,ie,n;
 mupp=fabs(hmass_pre+4)+4;
 mupp*=mupp;

 //Eigen Stuff
 int MVM=0; // counter for matrix-vector multiplications 

 max_H(&H2_pre, &glat_even, &max);
 //lprintf("MAIN",0,"MAXCHECK: cnfg=%e  uppbound=%e diff=%e %s\n",max,mupp,mupp-max,(mupp-max)<0?"[FAILED]":"[OK]");
 max=1.1*max;

 ie=eva_tuned(nev,nevt,0,kmax,maxiter,lbnd,max,omega1,omega2,&H2_pre,eva_vec,eva_val,&status);
 MVM+=status;
 while (ie!=0) { // if failed restart EVA 
  lprintf("MAIN",0,"Restarting EVA!\n");
  ie=eva_tuned(nev,nevt,2,kmax,maxiter,lbnd,max,omega1,omega2,&H2_pre,eva_vec,eva_val,&status);
  MVM+=status;
 }
 lprintf("MAIN",0,"EVA MVM = %d\n",MVM);
 neigs = nev;

 for (n=0;n<nev;++n) {
   H2_pre(tmp_sf,&eva_vec[n]);
   lprintf("RESULT",0,"Eig %d = %.15e %.15e\n",n,eva_val[n],
   spinor_field_prod_re_f(tmp_sf,&eva_vec[n])/spinor_field_sqnorm_f(&eva_vec[n]));
 } 

}

void copy_evec( int n, spinor_field* psi1, double *eval ){

	spinor_field_copy_f(psi1, &eva_vec[n]);
	*eval = eva_val[n];

}

void calc_deflated_propagator(spinor_field *psi, spinor_field* eta, int ndilute, int Nuse){
  int beta,i,n_masses,n;
  double *m;
  m = mass;
  n_masses=QMR_par.n;
  QMR_par.n=1;
  if(Nuse < 0 || Nuse > neigs){ Nuse = neigs; }
  lprintf("CALC_DEFLATED_PROPAGATOR",10,"n masses=%d, mass = %g, neigs = %d\n",n_masses, mass[0],Nuse);
  for (i=0;i<n_masses;++i){
    hmass_pre = mass[0];
    for (beta=0;beta<ndilute;++beta){
	psi[beta*n_masses+i].type = &glat_even; //even guy
	eta[beta].type = &glat_even; //even guy
        spinor_field_zero_f(&psi[beta*n_masses+i]);
        Ddag_pre(tmp, &eta[beta], tmp_sf);
        spinor_field_mul_f(tmp,(4.+m[0]),tmp);
	for (n=0;n<Nuse;++n) {
	  complex p = spinor_field_prod_f(&eva_vec[n],tmp);
	  _complex_mulr( p, ( 1./eva_val[n] ), p );
	  spinor_field_mulc_add_assign_f(&psi[beta*n_masses+i],p,&eva_vec[n]);
	}
        calc_propagator_core(&psi[beta*n_masses+i],&eta[beta],_MINRES);
    }
    mass++;
  }
  QMR_par.n = n_masses;
  mass = m;
  hmass_pre = mass[0];
}

#undef GAUSSIAN_NOISE
