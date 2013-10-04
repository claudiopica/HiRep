/***************************************************************************\
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen                              *
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
#include "gamma_spinor.h"

#define PI 3.141592653589793238462643383279502884197

//Helps QMR solver find more accurate solutions
#undef GAUSSIAN_NOISE 


/* Random timeslize not previously chosen */
static int random_tau(){
  static int* slices=NULL;
  if (slices == NULL) slices = (int*) malloc(GLB_T*sizeof(int));
  static int counter = 0;
  int itmp,tau,i;
  double ran;

  if (counter == 0){
    for (i=0;i<GLB_T;++i){
      slices[i]=i;
    }
    counter=GLB_T;
  }
  do{
    ranlxd(&ran,1);
    itmp=(int)(ran*counter);    
  } while(itmp==counter);
  counter--;
  tau = slices[itmp];
  slices[itmp]=slices[counter];
  slices[counter]=tau;
  bcast_int(&tau,1);
  return tau;
}

/***************************************************************************\

	Sources: 
		point_source: 			
						source[spin](t,x) = \delta_{a color} \delta_{s, spin} \delta( (t,x) - (tau,0) )
		diluted_source_equal_eo:		
						\xi(x) = Z(2) x Z(2)  -  NF color vector at x
						eta(t,x) = \delta(t - tau) \xi(x)
						source[spin](x) = \delta_{s spin} eta(t,x)  -  x even
		diluted_source_equal:		
						\xi(x) = Z(2) x Z(2)  -  NF color vector at x
						eta(t,x) = \delta(t - tau) \xi(x)
						source[spin](t,x) = \delta_{s spin} eta(t,x)  -  x even & odd
		noise_source_equal_eo:
						\xi(x) = Z(2) x Z(2)  -  NF color vector at x
						eta(t,x) = \xi(t,x)
						source[spin](t,x) = \delta_{s spin} eta(t,x)  -  x even
		gauge_fixed_wall_source:
						source[spin](t,x) = \delta_{a color} \delta_{s spin} \delta(t - tau) 1
		sequential_source:
						source[spin](tf,ti,x) = \gamma_5 S(x,tf; 0,ti)
		gauge_fixed_momentum_source:
						source[spin](t,x) = \delta_{a color} \delta_{s spin} e^{ i p_\mu x_\mu }
\***************************************************************************/
void create_point_source(spinor_field *source,int tau, int color) {
  int beta;
  for (beta=0;beta<4;++beta){
    spinor_field_zero_f(&source[beta]);
  }
  if(COORD[0]==tau/T && COORD[1]==0 && COORD[2]==0 && COORD[3]==0) {
    int ix=ipt(tau,0,0,0);
    for (beta=0;beta<4;++beta){
      _FIELD_AT(&source[beta],ix)->c[beta].c[color].re = 1.;
    }
  }
}

/* Creates four Z2xZ2 noise sources localised on time slice tau. The noise 
   vectors are equal in each source but placed at a different spin. */

int create_diluted_source_equal_eo(spinor_field *source) {
  int c[4];
  suNf_vector *v1,*v2;
  int i;
  int tau = random_tau();
  for (i=0;i<4;++i){
    spinor_field_zero_f(&source[i]);
  }
  
  if(COORD[0]==tau/T) {// Check that tau is in this thread.
    c[0]=tau%T;
    for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	  if(((tau+zerocoord[1]+c[1]+zerocoord[2]+c[2]+zerocoord[3]+c[3])&1)==0){
	    v1 = &((_FIELD_AT(&source[0],ipt(c[0],c[1],c[2],c[3])))->c[0]);
	    ranz2((double*)(v1),sizeof(suNf_vector)/sizeof(double)); // Make new sources
	    for (i=1;i<4;++i){
	      v2 = &((_FIELD_AT(&source[i],ipt(c[0],c[1],c[2],c[3])))->c[i]); //Copy previous index.
	      *v2 = *v1;
	    }
	  }
	}
  }
  return tau;
}

void create_diluted_source_equal_atau_eo(spinor_field *source, int tau){
  int c[4];
  suNf_vector *v1,*v2;
  int i;
  //int tau = random_tau();
  for (i=0;i<4;++i){
    spinor_field_zero_f(&source[i]);
  }
  
  if(COORD[0]==tau/T) {// Check that tau is in this thread.
    c[0]=tau%T;
    for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	  if(((tau+zerocoord[1]+c[1]+zerocoord[2]+c[2]+zerocoord[3]+c[3])&1)==0){
	    v1 = &((_FIELD_AT(&source[0],ipt(c[0],c[1],c[2],c[3])))->c[0]);
	    ranz2((double*)(v1),sizeof(suNf_vector)/sizeof(double)); // Make new sources
	    for (i=1;i<4;++i){
	      v2 = &((_FIELD_AT(&source[i],ipt(c[0],c[1],c[2],c[3])))->c[i]); //Copy previous index.
	      *v2 = *v1;
	    }
	  }
	}
  }
  //return tau;
}

int create_diluted_source_equal(spinor_field *source) {
  int c[4];
  suNf_vector *v1,*v2;
  int i;
  int tau = random_tau();
  for (i=0;i<4;++i){
    spinor_field_zero_f(&source[i]);
  }
  
  if(COORD[0]==tau/T) {// Check that tau is in this thread.
    c[0]=tau%T;
    for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	  v1 = &((_FIELD_AT(&source[0],ipt(c[0],c[1],c[2],c[3])))->c[0]);
	  ranz2((double*)(v1),sizeof(suNf_vector)/sizeof(double)); // Make new sources
	  for (i=1;i<4;++i){
	    v2 = &((_FIELD_AT(&source[i],ipt(c[0],c[1],c[2],c[3])))->c[i]); //Copy previous index.
	    *v2 = *v1;
	  }
	}
  }
  return tau;
}

/* Creates four Z2xZ2 noise sources NOT localised on time slice but spread over
all timeslices. The noise vectors are equal in each source but placed at a different spin. */

void create_noise_source_equal_eo(spinor_field *source) {
  int c[4];
  suNf_vector *v1,*v2;
  int i;
  //int tau = random_tau();
  for (i=0;i<4;++i){
    spinor_field_zero_f(&source[i]);
  }
  
  //if(COORD[0]==tau/T) {// Check that tau is in this thread.
  //  c[0]=tau%T;
    for(c[0]=0; c[0]<T; c[0]++) for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	  if(((zerocoord[0]+c[0]+zerocoord[1]+c[1]+zerocoord[2]+c[2]+zerocoord[3]+c[3])&1)==0){
	    v1 = &((_FIELD_AT(&source[0],ipt(c[0],c[1],c[2],c[3])))->c[0]);
	    ranz2((double*)(v1),sizeof(suNf_vector)/sizeof(double)); // Make new sources
	    for (i=1;i<4;++i){
	      v2 = &((_FIELD_AT(&source[i],ipt(c[0],c[1],c[2],c[3])))->c[i]); //Copy previous index.
	      *v2 = *v1;
	    }
	  }
	}
  //}
  //return tau;
}

//create a wall source at timeslice tau, all parity sites.
void create_gauge_fixed_wall_source(spinor_field *source, int tau, int color) {
  int c[4];
  int beta;

  for (beta=0;beta<4;++beta){
    spinor_field_zero_f(&source[beta]);
  }

  if(COORD[0]==tau/T) {// Check that tau is in this thread.
    c[0]=tau%T;
    for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	  for (beta=0;beta<4;++beta){
	      _FIELD_AT(&source[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].re = 1.;
	  }
    }
  }

}

void create_sequential_source(spinor_field *source, int tf, spinor_field* prop){
  int c[4];
  int beta;
  suNf_spinor sp;
  for (beta=0;beta<4;++beta){
    spinor_field_zero_f(&source[beta]);
  }

  if(COORD[0]==tf/T) {// Check that tf is in this thread.
    c[0]=tf%T;
    for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	  for (beta=0;beta<4;++beta){
	    _spinor_g5_f( sp , *_FIELD_AT(&prop[beta], ipt(c[0],c[1],c[2],c[3]) ) );
	    _spinor_plus_f( *_FIELD_AT(&source[beta], ipt(c[0],c[1],c[2],c[3])) , sp);
	  }
	}
  }
}

//create a e^ipx source
void create_gauge_fixed_momentum_source(spinor_field *source, int pt, int px, int py, int pz, int color) {
  int c[4];
  int beta;
  double pdotx;

  for (beta=0;beta<4;++beta){
    spinor_field_zero_f(&source[beta]);
  }

  for(c[0]=0; c[0]<T; c[0]++) for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	  pdotx = 2.*PI*( c[0]*pt/GLB_T + c[1]*px/GLB_X + c[2]*py/GLB_Y + c[3]*pz/GLB_Z );
	  for (beta=0;beta<4;++beta){
	     _FIELD_AT(&source[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].re = cos(pdotx);
	     _FIELD_AT(&source[beta], ipt(c[0],c[1],c[2],c[3]) )->c[beta].c[color].im = sin(pdotx);
	  }
  }
}

static double hmass_pre;

static void D_pre(spinor_field *out, spinor_field *in){
  Dphi_eopre(hmass_pre,out,in);
}

static void H_pre(spinor_field *out, spinor_field *in){
  g5Dphi_eopre(hmass_pre,out,in);
}

static int init=0;
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
static int init_odd=0;
/* Initialises the propagator, nm is the number of masses for multimass solver, 
   m is the array of masses, and acc is the inverter accuracy
*/



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
  init=0;
}



/***************************************************************************\

psi = D^{-1} eta

\***************************************************************************/


static void calc_propagator_eo_core(spinor_field *psi, spinor_field *eta) {
  spinor_field qprop_mask;
  int i, cgiter=0;
  error(init==0,1,"calc_prop.c","z2semwall method not initialized!");

  /* add source */
#ifdef GAUSSIAN_NOISE
  spinor_field_add_f(tmp,eta,QMR_noise);
#else
  spinor_field_copy_f(tmp,eta);
#endif

  /* invert source */
  for(i=0;i<QMR_par.n;++i){
    spinor_field_zero_f(&resd[i]);
  }

  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, tmp, resd);

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
    Dphi_(&qprop_mask,&resd[i]);
    spinor_field_minus_f(&qprop_mask,&qprop_mask);
    if(i&1) ++cgiter; /* count only half of calls. works because the number of sources is even */
  }
  lprintf("CALC_PROP",10,"QMR_eo MVM = %d\n",cgiter);
}

static void calc_propagator_core(spinor_field *psi, spinor_field *eta) {
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
  /* invert source */
  spinor_field_zero_f(resd);

#ifdef GAUSSIAN_NOISE
  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, tmp, resd);
#else
  spinor_field_g5_f(tmp,tmp);
  cgiter+=MINRES_mshift(&QMR_par, &H_pre, tmp, resd);  
  //  cgiter+=HBiCGstab_mshift(&QMR_par, &H_pre, tmp, resd);  
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
}


void calc_propagator_eo(spinor_field *psi, spinor_field *eta, int ndilute) {
  int beta;
  for (beta=0;beta<ndilute;++beta){
    calc_propagator_eo_core(&psi[beta*QMR_par.n],&eta[beta]);
  }
}

void calc_propagator(spinor_field *psi, spinor_field* eta, int ndilute){
  int beta,i,n_masses;
  double *m;
  m = mass;
  n_masses=QMR_par.n;
  QMR_par.n=1;
  lprintf("CALC_PROPAGATOR",10,"n masses=%d, mass = %g",n_masses, mass[0]);
  for (i=0;i<n_masses;++i){
    hmass_pre = mass[0];
    for (beta=0;beta<ndilute;++beta){
      calc_propagator_core(&psi[beta*n_masses+i],&eta[beta]);
    }
    mass++;
  }
  QMR_par.n = n_masses;
  mass = m;
  hmass_pre = mass[0];
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
      calc_propagator_core(&psi[beta*n_masses+i],&eta[beta*n_masses+i]);
      mass++;
    }
  }
  QMR_par.n = n_masses;
  mass = m;
  hmass_pre = mass[0];
}


#undef GAUSSIAN_NOISE
