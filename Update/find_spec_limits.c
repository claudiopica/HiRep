#include "update.h"
#include "inverters.h"
#include "linear_algebra.h"
#include "dirac.h"
#include "suN.h"
#include "random.h"
#include "global.h"
#include <math.h>
#include <malloc.h>
#include <stdio.h>

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

/* use power method to find min eigvalue */
int max_H2(double *max, double mass) {
  double norm, oldmax, dt;
  int count;
  suNf_spinor *s1,*s2,*s3;
  s1=(suNf_spinor*)malloc(sizeof(suNf_spinor)*3*VOLUME);
  s2=s1+VOLUME;
  s3=s2+VOLUME;

  /* random spinor */
  ranlxs((float*)s1,(sizeof(suNf_spinor)/sizeof(float))*VOLUME);
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

  printf("[H2 Max_eig = %1.8e][H2 Max_count = %d]\n",*max,count); 

  free(s1);
  
	return count;
}

int min_H2(double *min, double max, double mass) {
  double norm, oldmin;
  int count;
  suNf_spinor *s1,*s2,*s3;
  s1=(suNf_spinor*)malloc(sizeof(suNf_spinor)*4*VOLUME);
  s2=s1+VOLUME;
  s3=s2+VOLUME;  

  /* random spinor */
  ranlxs((float*)s1,(sizeof(suNf_spinor)/sizeof(float))*VOLUME);
  norm=sqrt(spinor_field_sqnorm_f(s1));
  spinor_field_mul_f(s1,1./norm,s1);
  
  if(max==0.) {
    max=4.+fabs(4.+mass);
    max*=max;
  }
  
  g5Dphi(mass,s2,s1);
  g5Dphi(mass,s3,s2);
    
  count=1;
  do {
    /* multiply vector by max-H2 */
    ++count;

    spinor_field_mul_f(s1,max,s1);
    spinor_field_sub_assign_f(s1,s3);

    norm=sqrt(spinor_field_sqnorm_f(s1));
    spinor_field_mul_f(s1,1./norm,s1);
    
    /* check the eigen value */
    g5Dphi(mass,s2,s1);
    g5Dphi(mass,s3,s2);

    oldmin=*min;
    *min = spinor_field_prod_re_f(s1,s3);

    /* printf("Iter %d: %4.5e\n",count,fabs(oldnorm-norm)); */
  } while (fabs((oldmin-*min)/(*min))>1.e-4);

  *min*=0.95; /* do not know exact bound */

  printf("[H2 Min_eig = %1.8e][H2 Min_count = %d]\n",*min,count); 

  free(s1);  

	return count;
}

void new_min_H2(double *min, double mass) {
  int count;
  double lambda, a, b, g, delta, norm, oldl;
  suNf_spinor *x,*z,*s1,*s2;
  x=(suNf_spinor*)malloc(sizeof(suNf_spinor)*4*VOLUME);
  z=x+VOLUME;
  s1=z+VOLUME;  
  s2=s1+VOLUME;  

  /* random spinor */
  ranlxs((float*)x,(sizeof(suNf_spinor)/sizeof(float))*VOLUME);
  norm=sqrt(spinor_field_sqnorm_f(x));
  spinor_field_mul_f(x,1./norm,x);

  g5Dphi(mass,s2,x);
  g5Dphi(mass,s1,s2);
  
  lambda=spinor_field_prod_re_f(x,s1);

  spinor_field_mul_f(z,lambda,x);
  spinor_field_sub_assign_f(z,s1);

  delta=sqrt(spinor_field_sqnorm_f(z));

  count=1;
  do{
    g5Dphi(mass,s2,s1);
    g5Dphi(mass,s1,s2);
    a=spinor_field_prod_re_f(x,s1)-lambda*lambda;
    g5Dphi(mass,s2,s1);
    g5Dphi(mass,s1,s2);
    b=spinor_field_prod_re_f(x,s1)-(3.*a+lambda*lambda)*lambda;
    g=b/(2.*a);
    g=(sqrt(g*g+a)-g)/a;
    spinor_field_mul_add_assign_f(x,g,z);
    norm=sqrt(spinor_field_sqnorm_f(x));
    spinor_field_mul_f(x,1./norm,x);
    
    g5Dphi(mass,s2,x);
    g5Dphi(mass,s1,s2);
    oldl=lambda;
    lambda=spinor_field_prod_re_f(x,s1);

    spinor_field_mul_f(z,lambda,x);
    spinor_field_sub_assign_f(z,s1);

    delta=sqrt(spinor_field_sqnorm_f(z));
    count+=3;
    
  } while (fabs((oldl-lambda)/lambda)>1.e-4);

  printf("NEW [H2 eig] %d Min: %e\n",count, lambda); 
  *min=0.95*lambda;

  free(x);
  
}

/* use power method to find min eigvalue */
void test_min_H2(double *max, double mass) {
  double norm, oldmax, dt, m;
  int count;
  suNf_spinor *s1,*s2,*s3, *s4;
  s1=(suNf_spinor*)malloc(sizeof(suNf_spinor)*4*VOLUME);
  s2=s1+VOLUME;
  s3=s2+VOLUME;
  s4=s3+VOLUME;

  /* random spinor */
  ranlxs((float*)s1,(sizeof(suNf_spinor)/sizeof(float))*VOLUME);
  norm=sqrt(spinor_field_sqnorm_f(s1));
  spinor_field_mul_f(s1,1./norm,s1);
  norm=1.;

  dt=1.;

  m=64.;

  g5Dphi(mass,s2,s1);
  g5Dphi(mass,s3,s2);
  spinor_field_mul_add_assign_f(s3,-m,s1);
  g5Dphi(mass,s2,s3);
  g5Dphi(mass,s1,s2);
  spinor_field_mul_add_assign_f(s1,-m,s3);
  g5Dphi(mass,s2,s1);
  g5Dphi(mass,s3,s2);
  spinor_field_mul_add_assign_f(s3,-m,s1);

  count=3;
  do {
    /* multiply vector by H2 */
    count+=3;

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
    spinor_field_mul_add_assign_f(s3,-m,s1);
    g5Dphi(mass,s2,s3);
    g5Dphi(mass,s4,s2);
    spinor_field_mul_add_assign_f(s4,-m,s3);
    g5Dphi(mass,s2,s4);
    g5Dphi(mass,s3,s2);
    spinor_field_mul_add_assign_f(s3,-m,s4);

    *max=spinor_field_prod_re_f(s1,s3);

    /* printf("Iter %d: %4.5e\n",count,fabs(oldnorm-norm)); */
  } while (fabs((*max-oldmax)/(*max))>1.e-4);

  *max = m - pow(fabs(*max),1./3.);

  printf("[H2 eig] %d Test Min: %e\n",count, *max); 

  free(s1);
  
}


static suNf_spinor **ev;

void find_spec_H2(double *max, double *min, double mass) {
	const int nevt=5;
	const int nev=1;
	const int kmax=200;
	const int maxiter=20;
	static float *d1;
	const float omega1=1.e-4; /* absolute precision */
	const float omega2=1.e-1; /* relative precision */
	int status,ie;
	suNf_spinor **ws;

	max_H2(max, mass);

	d1=malloc(sizeof(*d1)*nevt);
	ev=malloc(sizeof(*ev)*(nevt+2));
	ev[0]=malloc(sizeof(**ev)*(nevt+2)*VOLUME);
	for (status=1;status<nevt+2;++status) {
		ev[status]=ev[status-1]+VOLUME;
	}
	ws=ev+nevt;

	ie=eva(VOLUME,nev,nevt,0,kmax,maxiter,1,*max,omega1,omega2,&H2,ws,ev,d1,&status);
	printf("[EVA] status: %d, ret: %d\n",status, ie);

	*min=d1[0]*(1-omega2)-omega1;


	free(d1);
	free(ev[0]);
	free(ev);

	return;
}


