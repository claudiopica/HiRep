/* arXiv:1006.4518 [hep-lat] */

#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "suN_repr_func.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "update.h"
#include "utils.h"
#include "random.h"
#include "communications.h"
#include "wilsonflow.h"
#include <math.h>

#define REUNIT 10

void unit_gauge(suNg_field *gauge){
  int mu;
    int x, y, z, t, ix; for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ ix=ipt(t,x,y,z);
	    for (mu=0; mu<4; ++mu) {
		_suNg_unit(*_4FIELD_AT(gauge,ix,mu));
	    }
    }
  start_gf_sendrecv(gauge);
  complete_gf_sendrecv(gauge);
}


suNg_field* g;
static int init = 0;

static void init_g(){
	if(!init){
		g=alloc_gfield(&glattice);
		unit_gauge(g);
	}
}
static void free_g(){
	if(init){
		free_gfield(g);
	}
}

static void gUgmu(suNg_field *gauge){
     suNg w1, w2;
     int t, x, y, z, ix, mu;
     for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ ix=ipt(t,x,y,z);
	for (mu=0; mu<4; ++mu) {

		_suNg_times_suNg(w1,(*_4FIELD_AT(g,ix,0)),(*_4FIELD_AT(gauge,ix,mu)));
		_suNg_times_suNg_dagger(w2,w1,(*_4FIELD_AT(g,iup(ix,mu),0))); 
		*_4FIELD_AT(gauge,ix,mu)=w2;
	
	} 
    }
    start_gf_sendrecv(gauge);
    complete_gf_sendrecv(gauge);
}



void random_gauge_transform(suNg_field *gauge){

	init_g();
	suNg w1;
	int x, y, z, t, ix;

	for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ ix=ipt(t,x,y,z);
		random_suNg(&w1);
		*_4FIELD_AT(g,ix,0)=w1;
	}

	start_gf_sendrecv(g);
	complete_gf_sendrecv(g);
	gUgmu( gauge );
}

double calc_plaq(suNg_field* V){
  
int mu, nu;
int iy,iz;
suNg *v1,*v2,*v3,*v4,w1,w2,w3;
double pl, E = 0;

int t, x, y, z, ix; for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ix=ipt(t,x,y,z);
	for(mu=0;mu<4;mu++) for(nu=mu+1;nu<4;nu++) {
	
	  iy=iup(ix,mu);
	  iz=iup(ix,nu);

	  v1=_4FIELD_AT(V,ix,mu);
	  v2=_4FIELD_AT(V,iy,nu);
	  v3=_4FIELD_AT(V,iz,mu);
	  v4=_4FIELD_AT(V,ix,nu);

	  _suNg_times_suNg(w1,(*v1),(*v2));
	  _suNg_times_suNg(w2,(*v4),(*v3));
	  _suNg_times_suNg_dagger(w3,w1,w2);      
	      
	  _suNg_trace_re(pl,w3);

	#ifdef PLAQ_WEIGHTS
	  if(plaq_weight!=NULL) pl*=plaq_weight[ix*16+mu*4+nu];
	#endif
	  E += pl;
	}
}

  global_sum(&E, 1);
  return E/(6.*NG)/GLB_VOLUME;
}

void reunit(suNg_field *fixed_gauge){

suNg *u;
int mu;
int t, x, y, z, ix; for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ix=ipt(t,x,y,z);
	for(mu=0;mu<4;mu++){ 
		u=_4FIELD_AT(fixed_gauge,ix,mu);		//u_mu(x)
		project_to_suNg(u);
	}
}
start_gf_sendrecv(fixed_gauge);
complete_gf_sendrecv(fixed_gauge);

}

double gaugefix_action(int fix_dir, suNg_field *gauge ){
  suNg *u1;
  double action = 0.0;
  double retr;
  int mu;

  int t, x, y, z, ix; for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ix=ipt(t,x,y,z);
    for(mu=0;mu<4;mu++){
	if(mu != fix_dir){
		u1 = _4FIELD_AT(gauge,ix,mu);
		_suNg_trace_re(retr,(*u1) );
		action += retr;
	}
    }
  }
  int ndir = 0; for(mu=0;mu<4;mu++){ if(mu != fix_dir){ndir++;} }
  global_sum(&action,1);
  action *= 1./(NG*ndir*GLB_T*GLB_X*GLB_Y*GLB_Z);
  return action;
}

void su2_hit(int fix_dir, int parity, double overrelax, suNg_field *fixed_gauge, int c )
{
  suNg *u1, *u2, v1, v2;
  int p;   
  double r[3], a[3], asq;
#ifndef GAUGE_SON
  double r0, a0, a0sq, ax, ar, axdr;
#endif
  double rnorm;
  int t,x,y,z,idx,mu,i,j;

  //Get the SU(2) index of the SU(N) matrix. Copied from chroma.
  int i1, i2;
  int found = 0;
  int del_i = 0;
  int index = -1;

  while ( del_i < (NG-1) && found == 0 )
    {
      del_i++;
      for ( i1 = 0; i1 < (NG-del_i); i1++ )
	{
	  index++;
	  if ( index == c )
	    {
	      found = 1;
	      break;
	    }
	}
    }
  i2 = i1 + del_i;

  int i11 = i1*NG + i1;
  int i12 = i1*NG + i2;
  int i21 = i2*NG + i1;
  int i22 = i2*NG + i2;

  for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++){ idx=ipt(t,x,y,z);
	  p = (zerocoord[0] + t + zerocoord[1] + x + zerocoord[2] + y + zerocoord[3] + z)%2; 
	  if(p == parity){

	    _suNg_zero(v1);
	    for(mu=0;mu<4;mu++){ 
	      if(mu != fix_dir){
		u1=_4FIELD_AT(fixed_gauge,idx,mu);		//u_mu(x)
		u2=_4FIELD_AT(fixed_gauge,idn(idx,mu),mu);	//u_mu(x - mu)

		//v1 = u1 + u2^dag. v1 is not in SU(NG). 
		for(i=0;i<NG;i++){ for(j=0;j<NG;j++){
#ifdef GAUGE_SON
		    v1.c[i*NG + j]+= u1->c[i*NG+j]+u2->c[j*NG+i];
#else
		    _complex_add_star_assign(v1.c[i*NG + j], (u1->c[i*NG+j]), (u2->c[j*NG+i]) );
#endif
		  }}
	      }}
#ifdef GAUGE_SON
	    //Extract gauge transformation
	    r[0] = v1.c[i11] + v1.c[i22];
	    r[1] = v1.c[i12] - v1.c[i21];
	    rnorm=1./sqrt(r[0]*r[0]+r[1]*r[1]);
	    if(1./rnorm > 1e-16){
	      a[0] = r[0]*rnorm;
	      a[1] = r[1]*rnorm;
	    } else {
	      a[0] = 1;
	      a[1] = 0;
	    }
	    //Overrelax
	    asq = ( acos(a[0]) )*overrelax;
	    a[0] = cos(asq);
	    a[1] = sin(asq);
	    if(r[1] < 0) a[1] *= -1; 

	    _suNg_unit(v2);
	    v2.c[i11]=v2.c[i22]=(a[0]);
	    v2.c[i21]=(a[1]);
	    v2.c[i12]=-v2.c[i21];
#else		
	    //Extract gauge transform
	    r0 =   v1.c[i11].re + v1.c[i22].re;
	    r[0] = v1.c[i12].im + v1.c[i21].im; 
	    r[1] = v1.c[i12].re - v1.c[i21].re; 
	    r[2] = v1.c[i11].im - v1.c[i22].im; 
	    rnorm = 1./sqrt(r0*r0 + r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	    if(1./rnorm > 1e-16){
	      a0 = r0*rnorm;
	      a[0] = -r[0]*rnorm;
	      a[1] = -r[1]*rnorm;
	      a[2] = -r[2]*rnorm;
	    } else {
	      a0 = 1;
	      a[0] = 0;
	      a[1] = 0;
	      a[2] = 0;
	    }
	    //Overrelaxation
	    asq = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
	    a0sq = a0*a0;
	    ax = (overrelax*a0sq + asq)/(a0sq + asq);
	    ar = sqrt((double)(a0sq + ax*ax*asq));
	    axdr = ax/ar;
	    a0 = a0/ar; a[0] = a[0]*axdr; a[1] = a[1]*axdr; a[2] = a[2]*axdr;

	    _suNg_unit(v2);
	    v2.c[i11].re = a0; 	v2.c[i11].im = a[2]; 
	    v2.c[i12].re = a[1]; 	v2.c[i12].im = a[0]; 
	    v2.c[i21].re =-a[1]; 	v2.c[i21].im = a[0]; 
	    v2.c[i22].re = a0; 	v2.c[i22].im =-a[2]; 
#endif		

	    //_suNg_times_suNg( (*_4FIELD_AT(g,idx,0)), (*_4FIELD_AT(g,idx,0)), v2);
	    *_4FIELD_AT(g,idx,0) = v2;
	  } 
	}
     start_gf_sendrecv(g);
     complete_gf_sendrecv(g);
} 



double gaugefixstep(int fix_dir,double overrelax, suNg_field *fixed_gauge )
{
	int c, parity;
	for(parity=0; parity<2; parity++){
		//Loop over SU(2) subgroups
		for(c=0; c< NG*(NG-1)/2; ++c){
			unit_gauge(g);
			su2_hit(fix_dir, parity, overrelax, fixed_gauge, c); //<-- Should Fix This Maybe?
			gUgmu(fixed_gauge);
		}
	}
	return gaugefix_action(fix_dir, fixed_gauge );
} 

double gaugefix(int fix_dir,double overrelax,int max_it,
	      double fix_tol, suNg_field *fixed_gauge )
{
  lprintf("GAUGE FIXING",30,"Reunitarizing every %d steps.\n",REUNIT);
  int it;
  double new_act, old_act=0., diff_act;
  init_g();
  for (it=0; it < max_it; it++) {
 	//Do gauge fixing  
    	new_act = gaugefixstep(fix_dir, overrelax, fixed_gauge); //returns new action
	double p2 = calc_plaq(fixed_gauge);
        lprintf("GAUGE FIXING",20,"Iteration %d action %1.12f plaq %1.6f\n",it, new_act, p2);

	if(it != 0) {
	  diff_act = new_act - old_act;
	  if (fabs(diff_act) < fix_tol) break;		//Stop when the change in action is very small
	}
    	old_act = new_act;

	// Reunitarize periodically
	if((it % REUNIT) == (REUNIT - 1)){ reunit(fixed_gauge); }
  }
  // Reunitarize at the end
  if((it % REUNIT) != 0){ reunit(fixed_gauge); }

  start_gf_sendrecv(fixed_gauge);
  complete_gf_sendrecv(fixed_gauge);

  new_act = gaugefix_action(fix_dir, fixed_gauge);
  free_g();
  return new_act;
}


