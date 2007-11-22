#include "update.h"
#include "inverters.h"
#include "linear_algebra.h"
#include "dirac.h"
#include "suN.h"
#include "random.h"
#include "memory.h"
#include "global.h"
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <logger.h>

/* Neuberger bound. See hep-lat/9911004 - Phys.Rev.D61:085015,2000. */
/* find a lower bound for the minimal eigenvalue of H(-1); */
/*
static void min_bound_H2(double *min) {
	int ix,iy,iz,mu,nu;
	suNg *v1,*v2,*v3,*v4,w1,w2,w3;
	double eps,max;
	const double c=(2.+sqrt(2.));

	*min = 0.;
	for (mu=1;mu<4;++mu) {
		for (nu=0;nu<mu;++nu) {
			max=0.;
			for (ix=0; ix<VOLUME; ++ix) {
				iy=iup[ix][mu];
				iz=iup[ix][nu];

				v1=pu_gauge(ix,mu);
				v2=pu_gauge(iy,nu);
				v3=pu_gauge(iz,mu);
				v4=pu_gauge(ix,nu);

				_suNg_times_suNg(w1,(*v1),(*v2));
				_suNg_times_suNg(w2,(*v4),(*v3));
				_suNg_times_suNg_dagger(w3,w1,w2);

				eps=sqrt(_suNg_sqnorm_m1(w3));
				if(max<eps) max=eps;
			}
			*min+=max;
		}
	}

	*min=1.-c*(*min);
	if ((*min)<0.) *min=0.;

}
*/

/* use power method to find max eigvalue */
int max_H2(double *max, double mass) {
  double norm, oldmax, dt;
  int count, ix, sx;
  spinor_field *s1,*s2,*s3;
	unsigned int len;

	get_spinor_len(&len);
  s1=alloc_spinor_field_f(3);
  s2=s1+1;
  s3=s2+1;

  /* random spinor */
   FOR_LOCAL_SC(ix,sx) {
      suNf_spinor* sptr=_SPINOR_AT(s1,sx);
      ranlxd((double*)sptr,sizeof(suNf_spinor)/sizeof(double));
   }
  norm=sqrt(spinor_field_sqnorm_f(s1));
  spinor_field_mul_f(s1,1./norm,s1);
  norm=1.;

  dt=1.;

  g5Dphi(mass,s2,s1);
  g5Dphi(mass,s3,s2);

  count=1;
  do {
    /* multiply vector by H2 */
    ++count;

    spinor_field_mul_f(s1,dt,s3);

    /*
    g5Dphi(mass,s2,s3);
    g5Dphi(mass,s3,s2);
    
    spinor_field_mul_add_assign_f(s1,dt*dt*0.5,s3);

    g5Dphi(mass,s2,s3);
    g5Dphi(mass,s3,s2);
    
    spinor_field_mul_add_assign_f(s1,dt*dt*dt/6.,s3);
    */

    norm=sqrt(spinor_field_sqnorm_f(s1));
    spinor_field_mul_f(s1,1./norm,s1);


    oldmax=*max;
    g5Dphi(mass,s2,s1);
    g5Dphi(mass,s3,s2);
    *max=spinor_field_prod_re_f(s1,s3);
    
    /* printf("Iter %d: %4.5e\n",count,fabs(oldnorm-norm)); */
  } while (fabs((*max-oldmax)/(*max))>1.e-3);

  *max*=1.1; /* do not know exact bound */

  lprintf("MaxH2",10,"Max_eig = %1.8e [MVM = %d]\n",*max,count); 

  free_field(s1);
  
	return count;
}

static spinor_field *ev;

void find_spec_H2(double *max, double *min, double mass) {
	/* EVA parameters */
	const int nevt=5; /* use 5-dim space */
	const int nev=1; /* require only the smallest to be accurate */
	const int kmax=200; /* max degree of polynomial */
	const int maxiter=20; /* max number of subiterations */
	static double *d1;
	const double omega1=1.e-4; /* absolute precision */
	const double omega2=1.e-1; /* relative precision */
	int status,ie;
	suNf_spinor **ws;
	/* END of EVA parameters */
	int MVM=0; /* counter for matrix-vector multiplications */
	unsigned int len;


	MVM+=max_H2(max, mass);

	get_spinor_len(&len);

	d1=malloc(sizeof(*d1)*nevt);
	ev=alloc_spinor_field_f((nevt+2));
	ws=ev+nevt;

	ie=eva(len,nev,nevt,0,kmax,maxiter,*max,omega1,omega2,&H2,ws,ev,d1,&status);
	MVM+=status;
	while (ie!=0) { /* if failed restart EVA */
		ie=eva(len,nev,nevt,2,kmax,maxiter,*max,omega1,omega2,&H2,ws,ev,d1,&status);
		MVM+=status;
	}

	*min=d1[0]*(1-omega2)-omega1;

	lprintf("SPECLIMITS",0,"Range = [%1.8e,%1.8e] [MVM = %d]\n",*min,*max,MVM);

	free(d1);
	free_spinor_field(ev);

	return;
}


