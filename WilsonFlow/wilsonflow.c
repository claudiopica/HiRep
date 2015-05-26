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
#include "communications.h"
#include "wilsonflow.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

/*
  #define EXP_CHECK
  #define PLAQ_CHECK
*/

/*
 * d/dt V = Z(V) V
 * S_W = 1/g^2 \sum_{p oriented} Re tr ( 1 - V(p) )
 * d_L f(V) = T^a d/ds f(e^{sT^a} V)
 * T_a^dag = -T_a      tr T_a T_b = -delta_{ab}/2
 * Z(V) = -g^2 d_L S_W = \sum_{p oriented} d_L Re tr ( V(p) )
 * Z(V) = d_L Re ( tr V s + tr V^dag s^dag ) =
 *      = d_L ( tr V s + tr V^dag s^dag ) =
 *      = T^a tr T^a ( V s - s^dag V^dag )
 *      = - 1/2 ( V s - s^dag V^dag - 1/N tr (V s - s^dag V^dag) )
 */
static void Zeta(suNg_field *Z, const suNg_field* U, const double alpha){

  error(Z->type!=&glattice,1,"wilson_flow.c","'Z' in Zeta must be defined on the whole lattice");
  error(U->type!=&glattice,1,"wilson_flow.c","'U' in Zeta must be defined on the whole lattice");


#ifdef PLAQ_CHECK
  double plaq=0.;
  _MASTER_FOR_SUM(&glattice,i,plaq) {
#else
  _MASTER_FOR(&glattice,i) {
#endif
    suNg staple, tmp1, tmp2;
    for (int mu=0; mu<4; ++mu) {
      _suNg_zero(staple);
      for(int nu=(mu+1)%4; nu!=mu; nu=(nu+1)%4) {
        suNg *u1=_4FIELD_AT(U,iup(i,mu),nu);
        suNg *u2=_4FIELD_AT(U,iup(i,nu),mu);
        suNg *u3=_4FIELD_AT(U,i,nu);
        _suNg_times_suNg_dagger(tmp1,*u1,*u2);
        _suNg_times_suNg_dagger(tmp2,tmp1,*u3);
#ifdef PLAQ_WEIGHTS
        if(plaq_weight!=NULL) {
          _suNg_mul(tmp2,plaq_weight[i*16+nu*4+mu],tmp2);
        }
#endif
        _suNg_add_assign(staple,tmp2);
        
        int j=idn(i,nu);
        u1=_4FIELD_AT(U,iup(j,mu),nu);
        u2=_4FIELD_AT(U,j,mu);
        u3=_4FIELD_AT(U,j,nu);
        _suNg_times_suNg(tmp1,*u2,*u1);
        _suNg_dagger_times_suNg(tmp2,tmp1,*u3);

#ifdef PLAQ_WEIGHTS
        if(plaq_weight!=NULL) {
          _suNg_mul(tmp2,plaq_weight[j*16+nu*4+mu],tmp2);
        }
#endif
        _suNg_add_assign(staple,tmp2);
      }
      
      _suNg_times_suNg(tmp1,*_4FIELD_AT(U,i,mu),staple);
      
#ifdef PLAQ_CHECK
      double retr;
      _suNg_trace_re(retr,tmp1);
      plaq += retr/(NG*6);
#endif
      
      _suNg_dagger(tmp2,tmp1);
      _suNg_sub_assign(tmp1,tmp2);
#ifndef GAUGE_SON
#ifndef WITH_QUATERNIONS
      double imtr;
      _suNg_trace_im(imtr,tmp1);
      imtr = imtr/NG;
      for(int k=0;k<NG*NG;k+=NG+1) {
        tmp1.c[k].im -= imtr;
      }
#endif
/*
      for(int k=0; k<NG*NG; k++) {
      	_complex_mulr_assign(_4FIELD_AT(Z,i,mu)->c[k],-alpha/2.,tmp1.c[k]);
      }
*/
      _suNg_mul(tmp1,-alpha/2.,tmp1);
      _suNg_add_assign(*_4FIELD_AT(Z,i,mu),tmp1);
#else
/*
      for(int k=0; k<NG*NG; k++) {
      	_4FIELD_AT(Z,i,mu)->c[k]=-alpha/2.*tmp1.c[k];
      }  
*/
      _suNg_mul(tmp1,-alpha/2.,tmp1);
      _suNg_add_assign(*_4FIELD_AT(Z,i,mu),tmp1);
#endif
    }

  } //Master_for
 
#ifdef PLAQ_CHECK
  plaq /= (4*GLB_T*GLB_X*GLB_Y*GLB_Z);
  global_sum(&plaq,1);
  lprintf("WILSONFLOW",0,"ZETA Plaquette check %f %e\n",plaq,24.*(1.-plaq));
#endif
}


#ifdef EXP_CHECK
static void WF_Exp_check(suNg *u, suNg *X) {
  suNg Xk, tmp;
  _suNg_unit(*u);
  _suNg_unit(Xk);
  
  int k=1;
  double error;
  while(1) {
    _suNg_times_suNg(tmp,Xk,*X);
    _suNg_mul(Xk,1./k,tmp);
    k++;
    _suNg_add_assign(*u,Xk);    

    _suNg_sqnorm(error,Xk);
    if(error<1e-28) break;
  }
  
}
#endif

#if NG==2 && !defined(WITH_QUATERNIONS)

/*
 *  u = exp(X)
 *
 * I AM ASSUMING
 * X^dag = -X
 * tr X = 0
 */
static void WF_Exp(suNg *u, suNg *X) {
  suNg_algebra_vector h,v;

  h.c[0] = X->c[1].im;
  h.c[1] = X->c[1].re;
  h.c[2] = X->c[0].im;
  
  double z=sqrt(h.c[0]*h.c[0]+h.c[1]*h.c[1]+h.c[2]*h.c[2]);
  double s=1.;
  if(z>1e-16) s=sin(z)/z;
  double c=cos(z);
  v.c[0]=h.c[0]*s;
  v.c[1]=h.c[1]*s;
  v.c[2]=h.c[2]*s;

  u->c[0].re = c;       u->c[0].im = v.c[2];
  u->c[1].re = v.c[1];  u->c[1].im = v.c[0];
  u->c[2].re = -v.c[1]; u->c[2].im = v.c[0];
  u->c[3].re = c;       u->c[3].im = -v.c[2];
  
#ifdef EXP_CHECK
  suNg w;
  double error;
  WF_Exp_check(&w,X);
  _suNg_sub_assign(w,*u);
  _suNg_sqnorm(error,w);
  lprintf("WILSONFLOW",0,"WF EXP CHECK %e\n",sqrt(error));
#endif
}

/*#elif NG==3*/


/* static void WF_Exp(suNg *u, suNg *X) { */
/*   suNg X2; */
/*   complex c[3], s[3], tmp; */
/*   double alpha, beta; */
/*   double norm, error; */
/*   int n; */
  
  
/* /\* X2 = X.X *\/ */
/*   _suNg_times_suNg(X2,*X,*X); */
  
/* /\* alpha = Im det(X) *\/ */
/*   #define _X_(a,b) (X->c[a+3*b]) */
/*   #define ImProd(a,b,c) (a.re*(b.re*c.im+b.im*c.re)+a.im*(b.re*c.re-b.im*c.im)) */
/*   alpha= */
/*     +ImProd(_X_(0,0),_X_(1,1),_X_(2,2)) */
/*     +ImProd(_X_(0,1),_X_(1,2),_X_(2,0)) */
/*     +ImProd(_X_(0,2),_X_(1,0),_X_(2,1)) */
/*     -ImProd(_X_(0,2),_X_(1,1),_X_(2,0)) */
/*     -ImProd(_X_(0,1),_X_(1,0),_X_(2,2)) */
/*     -ImProd(_X_(0,0),_X_(1,2),_X_(2,1)); */
/*   #undef _X_ */
/*   #undef ImProd */
  
/* /\* beta = tr (X2) / 2 *\/ */
/* /\* norm = sqrt( |tr (X2)| ) *\/ */
/*   #define _X2_(a,b) (X2.c[a+3*b]) */
/*   beta=(_X2_(0,0).re+_X2_(1,1).re+_X2_(2,2).re)/2.; */
/*   norm=sqrt(-_X2_(0,0).re-_X2_(1,1).re-_X2_(2,2).re); */
/*   #undef _X2_ */
  
/*   s[0].re = 1.;          s[0].im = alpha/6.; */
/*   s[1].re = 1.+beta/6.;  s[1].im = 0.; */
/*   s[2].re = .5;          s[2].im = 0.; */
  
/*   n=3; */
/*   c[0].re = 0.;          c[0].im = alpha/6.; */
/*   c[1].re = beta/6.;     c[1].im = 0.; */
/*   c[2].re = 0.;          c[2].im = 0.; */
  
/*   /\* error = |X|^{n+1}/(n+1)! exp(|X|) *\/   */
/*   error= exp(norm)*norm*norm*norm*norm/24.; */
/* #error The error must be rechecked!!! */

/*   /\* */
/*   c[0][n] = i*c[2][n-1]*alpha/n */
/*   c[1][n] = (c[0][n-1]+c[2][n-1]*beta)/n */
/*   c[2][n] = c[1][n-1]/n */
/*   *\/ */
/*   while(1) { */
/*     n++; */
/*     tmp=c[1]; */
/*     c[1].re=(c[0].re+c[2].re*beta)/n; */
/*     c[1].im=(c[0].im+c[2].im*beta)/n; */
/*     c[0].re=-c[2].im*alpha/n; */
/*     c[0].im=c[2].re*alpha/n; */
/*     c[2].re=tmp.re/n; */
/*     c[2].im=tmp.im/n; */
    
/*     s[0].re+=c[0].re; s[0].im+=c[0].im; */
/*     s[1].re+=c[1].re; s[1].im+=c[1].im; */
/*     s[2].re+=c[2].re; s[2].im+=c[2].im; */

/*     error *= norm/(n+1); */
/*     if(error < 1.e-20) break; */
/*   } */
  
/*   _suNg_zero(*u); */
/*   u->c[0].re=s[0].re; u->c[0].im=s[0].im; */
/*   u->c[4].re=s[0].re; u->c[4].im=s[0].im; */
/*   u->c[8].re=s[0].re; u->c[8].im=s[0].im; */
/*   for(int i=0; i<9; i++) { */
/*     _complex_mul_assign(u->c[i],s[1],X->c[i]); */
/*     _complex_mul_assign(u->c[i],s[2],X2.c[i]); */
/*   } */
  
/* #ifdef EXP_CHECK */
/*   suNg v; */
/*   WF_Exp_check(&v,X); */
/*   _suNg_sub_assign(v,*u); */
/*   _suNg_sqnorm(error,v); */
/*   lprintf("WILSONFLOW",0,"WF EXP CHECK %e\n",sqrt(error)); */
/* #endif */
  
/* } */
  
#else

static void WF_Exp(suNg *u, suNg *X) {
  suNg Xk, tmp;
  _suNg_unit(*u);
  _suNg_unit(Xk);
  
  int k=1;
  double error;
  while(1) {
    _suNg_times_suNg(tmp,Xk,*X);
    _suNg_mul(Xk,1./k,tmp);
    k++;
    _suNg_add_assign(*u,Xk);    

    _suNg_sqnorm(error,Xk);
    if(error<1e-28) break;
  }
  
}

#endif




static suNg_field* ws_gf=NULL;
static suNg_field* ws_gf_tmp=NULL;
static suNg_field* Vprime=NULL;
static suNg_field* u_gauge_backup=NULL;



void WF_initialize() {
  if(ws_gf==NULL) 
	{		
		ws_gf=alloc_gfield(&glattice);
		ws_gf_tmp=alloc_gfield(&glattice);
		Vprime=alloc_gfield(&glattice);  	
		u_gauge_backup=alloc_gfield(&glattice);  	
}
}


void WF_free() {
  if(ws_gf!=NULL)
  {
		free_gfield(ws_gf);
		free_gfield(ws_gf_tmp);
		free_gfield(Vprime);
		free_gfield(u_gauge_backup);
	}
}


void WilsonFlow1(suNg_field* V, const double epsilon) {

  _MASTER_FOR(&glattice,ix) {
    for (int mu=0; mu<4; ++mu) {
      _suNg_zero(*_4FIELD_AT(ws_gf,ix,mu));
    }
  }

  Zeta(ws_gf,V,epsilon);
  
  _MASTER_FOR(&glattice,ix) {
    suNg utmp[2];
    for (int mu=0; mu<4; ++mu) {
      WF_Exp(&utmp[0],_4FIELD_AT(ws_gf,ix,mu));
      _suNg_times_suNg(utmp[1],utmp[0],*_4FIELD_AT(V,ix,mu));
      *_4FIELD_AT(V,ix,mu)=utmp[1];
    }
  }
  
  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

}


//define distance between complex matrices
double max_distance(suNg_field* V, suNg_field* Vprime){

double d,tmp;
suNg diff;
_suNg_zero(diff);



tmp=0;
_MASTER_FOR(&glattice,ix) {
	d=0.;

	for (int mu=0; mu<4; ++mu) {
	 _suNg_mul(diff,1.,*_4FIELD_AT(V,ix,mu));
	 _suNg_sub_assign(diff,*_4FIELD_AT(Vprime,ix,mu));
	 _suNg_sqnorm(d,diff);
	 if (d > tmp) tmp=d;

	}

}
global_max(&tmp,1);


return tmp/(double)NG;
}

// following 1301.4388
double WilsonFlow3_adaptative(suNg_field* V, double epsilon,double delta) {
 
suNg_field_copy(u_gauge_backup,u_gauge);
double varepsilon,d;
 
  _MASTER_FOR(&glattice,ix) {
    for (int mu=0; mu<4; ++mu) {
      _suNg_zero(*_4FIELD_AT(ws_gf,ix,mu));
      _suNg_zero(*_4FIELD_AT(ws_gf_tmp,ix,mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif


  Zeta(ws_gf,V,epsilon/4.); //ws_gf = Z0/4
  
  _MASTER_FOR(&glattice,ix) {
    suNg utmp[2];
    for (int mu=0; mu<4; ++mu) {
      WF_Exp(&utmp[0],_4FIELD_AT(ws_gf,ix,mu));
      _suNg_times_suNg(utmp[1],utmp[0],*_4FIELD_AT(V,ix,mu));
      *_4FIELD_AT(V,ix,mu)=utmp[1]; // V = exp(Z0/4) W0
      _suNg_mul(*_4FIELD_AT(ws_gf_tmp,ix,mu),-4.,*_4FIELD_AT(ws_gf,ix,mu)); //ws_gf_tmp = -Z0
      _suNg_mul(*_4FIELD_AT(ws_gf,ix,mu),-17./9.,*_4FIELD_AT(ws_gf,ix,mu)); //ws_gf =  -17*Z0/36 
    }
  }
  

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  Zeta(ws_gf,V,8.*epsilon/9.); // ws_gf = 8 Z1 /9 - 17 Z0/36
  Zeta(ws_gf_tmp,V,2.*epsilon); // ws_gf_tmp = 2 Z1 - Z0
  
  _MASTER_FOR(&glattice,ix) {
    suNg utmp[4];
    for (int mu=0; mu<4; ++mu) {
      WF_Exp(&utmp[0],_4FIELD_AT(ws_gf,ix,mu));
      WF_Exp(&utmp[2],_4FIELD_AT(ws_gf_tmp,ix,mu));
      _suNg_times_suNg(utmp[1],utmp[0],*_4FIELD_AT(V,ix,mu)); // utmp[1] = exp(8 Z1/9 - 17 Z0/36) W1
      _suNg_times_suNg(utmp[3],utmp[2],*_4FIELD_AT(V,ix,mu)); // utmp[4] = exp( Z1 -  Z0) W1
      *_4FIELD_AT(V,ix,mu)=utmp[1];
      *_4FIELD_AT(Vprime,ix,mu)=utmp[3];
      _suNg_mul(*_4FIELD_AT(ws_gf,ix,mu),-1.,*_4FIELD_AT(ws_gf,ix,mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);

  start_gf_sendrecv(Vprime);
  complete_gf_sendrecv(Vprime);

#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  
  Zeta(ws_gf,V,3.*epsilon/4.);
  
  _MASTER_FOR(&glattice,ix) {
    suNg utmp[2];
    for (int mu=0; mu<4; ++mu) {
      WF_Exp(&utmp[0],_4FIELD_AT(ws_gf,ix,mu));
      _suNg_times_suNg(utmp[1],utmp[0],*_4FIELD_AT(V,ix,mu));
      *_4FIELD_AT(V,ix,mu)=utmp[1];
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif


// now need to get the maximum of the distance 
	d = max_distance(V,Vprime);

// compute new value of epsilon	
	varepsilon = 0.95*epsilon*pow(delta/d,1./3.);

	if (varepsilon < 0.1) epsilon = varepsilon ;
	else epsilon = 0.1;

	if (d > delta ) 
	{
		suNg_field_copy(u_gauge,u_gauge_backup);
    epsilon = -1.;
		lprintf("WARNING",0,"d > delta ! Epsilon is set to -1 in order to repeat the calculation \n");
	}

	return epsilon;
}




void WilsonFlow3(suNg_field* V, const double epsilon) {

  
  _MASTER_FOR(&glattice,ix) {
    for (int mu=0; mu<4; ++mu) {
      _suNg_zero(*_4FIELD_AT(ws_gf,ix,mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif


  Zeta(ws_gf,V,epsilon/4.);
  
  _MASTER_FOR(&glattice,ix) {
    suNg utmp[2];
    for (int mu=0; mu<4; ++mu) {
      WF_Exp(&utmp[0],_4FIELD_AT(ws_gf,ix,mu));
      _suNg_times_suNg(utmp[1],utmp[0],*_4FIELD_AT(V,ix,mu));
      *_4FIELD_AT(V,ix,mu)=utmp[1];
      _suNg_mul(*_4FIELD_AT(ws_gf,ix,mu),-17./9.,*_4FIELD_AT(ws_gf,ix,mu));
    }
  }
  

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  Zeta(ws_gf,V,8.*epsilon/9.);
  
  _MASTER_FOR(&glattice,ix) {
    suNg utmp[2];
    for (int mu=0; mu<4; ++mu) {
      WF_Exp(&utmp[0],_4FIELD_AT(ws_gf,ix,mu));
      _suNg_times_suNg(utmp[1],utmp[0],*_4FIELD_AT(V,ix,mu));
      *_4FIELD_AT(V,ix,mu)=utmp[1];
      _suNg_mul(*_4FIELD_AT(ws_gf,ix,mu),-1.,*_4FIELD_AT(ws_gf,ix,mu));
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

  
  Zeta(ws_gf,V,3.*epsilon/4.);
  
  _MASTER_FOR(&glattice,ix) {
    suNg utmp[2];
    for (int mu=0; mu<4; ++mu) {
      WF_Exp(&utmp[0],_4FIELD_AT(ws_gf,ix,mu));
      _suNg_times_suNg(utmp[1],utmp[0],*_4FIELD_AT(V,ix,mu));
      *_4FIELD_AT(V,ix,mu)=utmp[1];
    }
  }

  start_gf_sendrecv(V);
  complete_gf_sendrecv(V);
#ifdef ROTATED_SF
  apply_BCs_on_fundamental_gauge_field();
#endif

}


static void WF_plaq(double *ret,suNg_field* V,int ix,int mu,int nu)
{
  int iy,iz;
  suNg *v1,*v2,*v3,*v4,w1,w2,w3;

  iy=iup(ix,mu);
  iz=iup(ix,nu);

  v1=_4FIELD_AT(V,ix,mu);
  v2=_4FIELD_AT(V,iy,nu);
  v3=_4FIELD_AT(V,iz,mu);
  v4=_4FIELD_AT(V,ix,nu);

  _suNg_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg(w2,(*v4),(*v3));
  _suNg_times_suNg_dagger(w3,w1,w2);      
      
  _suNg_trace_re(*ret,w3);

#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL)
    *ret *= plaq_weight[ix*16+nu*4+mu];
#endif


}


double WF_E(suNg_field* V) {
  double E=0.;

  _MASTER_FOR_SUM(&glattice,ix,E) {
    double p;
    for(int mu=0;mu<4;mu++) for(int nu=mu+1;nu<4;nu++) {
      WF_plaq(&p,V, ix, mu, nu);
      E += ((double)NG)-p;
    }
  }
  
  E *= 2./((double)GLB_VOLUME);
  
  global_sum(&E,1);
  
  return E;
}

void WF_E_T(double* E, suNg_field* V) {
  int t,x,y,z,ix;
  int mu,nu;
  double p;

  for(t=0;t<2*GLB_T;t++) E[t]=0.;

  for (t=0; t<T; t++){
    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++)
					      for(mu=0;mu<4;mu++) for(nu=mu+1;nu<4;nu++) {
						  ix=ipt(t,x,y,z);
						  WF_plaq(&p,V, ix, mu, nu);
						  if(mu==0) E[2*((zerocoord[0]+t+GLB_T)%GLB_T)] += p;
						  else E[2*((zerocoord[0]+t+GLB_T)%GLB_T)+1] += p;
						}
    E[2*((zerocoord[0]+t+GLB_T)%GLB_T)] = - E[2*((zerocoord[0]+t+GLB_T)%GLB_T)]/(3.*GLB_VOL3);
    E[2*((zerocoord[0]+t+GLB_T)%GLB_T)+1] = - E[2*((zerocoord[0]+t+GLB_T)%GLB_T)+1]/(3.*GLB_VOL3);
  }

  global_sum(E,2*GLB_T);

  for(t=0;t<2*GLB_T;t++) E[t] += NG;
  E[3]=0.0;
}


/* This gives F_{\mu\nu}^A */
static void WF_clover_F(suNg_algebra_vector *F, suNg_field* V, int ix, int mu, int nu) {
  int iy,iz,iw;
  suNg *v1,*v2,*v3,*v4,w1,w2,w3;

  _suNg_unit(w3);
  _suNg_mul(w3,-4.,w3);

  iy=iup(ix,mu);
  iz=iup(ix,nu);
  
  v1=_4FIELD_AT(V,ix,mu);
  v2=_4FIELD_AT(V,iy,nu);
  v3=_4FIELD_AT(V,iz,mu);
  v4=_4FIELD_AT(V,ix,nu);
  
  _suNg_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg_dagger(w2,w1,(*v3));
  _suNg_times_suNg_dagger(w1,w2,(*v4));
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(w1,plaq_weight[ix*16+nu*4+mu],w1);
  }
#endif
  _suNg_add_assign(w3,w1);
  
  iy=idn(ix,mu);
  iz=iup(iy,nu);
  
  v1=_4FIELD_AT(V,ix,nu);
  v2=_4FIELD_AT(V,iz,mu);
  v3=_4FIELD_AT(V,iy,nu);
  v4=_4FIELD_AT(V,iy,mu);
  
  _suNg_times_suNg_dagger(w1,(*v1),(*v2));
  _suNg_times_suNg_dagger(w2,w1,(*v3));
  _suNg_times_suNg(w1,w2,(*v4));
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(w1,plaq_weight[iy*16+nu*4+mu],w1);
  }
#endif
  _suNg_add_assign(w3,w1);
  
  iy=idn(ix,mu);
  iz=idn(iy,nu);
  iw=idn(ix,nu);
  
  v1=_4FIELD_AT(V,iy,mu);
  v2=_4FIELD_AT(V,iz,nu);
  v3=_4FIELD_AT(V,iz,mu);
  v4=_4FIELD_AT(V,iw,nu);
  
  _suNg_times_suNg(w1,(*v2),(*v1));
  _suNg_dagger_times_suNg(w2,w1,(*v3));
  _suNg_times_suNg(w1,w2,(*v4));
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(w1,plaq_weight[iz*16+nu*4+mu],w1);
  }
#endif
  _suNg_add_assign(w3,w1);

  iy=idn(ix,nu);
  iz=iup(iy,mu);
  
  v1=_4FIELD_AT(V,iy,nu);
  v2=_4FIELD_AT(V,iy,mu);
  v3=_4FIELD_AT(V,iz,nu);
  v4=_4FIELD_AT(V,ix,mu);
  
  _suNg_dagger_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg(w2,w1,(*v3));
  _suNg_times_suNg_dagger(w1,w2,(*v4));
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(w1,plaq_weight[iy*16+nu*4+mu],w1);
  }
#endif	
  _suNg_add_assign(w3,w1);
  
  _fund_algebra_project(*F,w3);
  
  _algebra_vector_mul_g(*F,1/4.,*F);
}


double WF_Esym(suNg_field* V) {
  double E=0.;
  
  _MASTER_FOR_SUM(&glattice,ix,E) {
    suNg_algebra_vector clover;
    double p;
    for(int mu=0;mu<4;mu++) for(int nu=mu+1;nu<4;nu++) {
      WF_clover_F(&clover,V,ix,mu,nu);
      _algebra_vector_sqnorm_g(p,clover);
      E += p;
    }
  }
  
  E *= _FUND_NORM2/((double)GLB_VOLUME);
  
  global_sum(&E,1);
  return E;
}



void WF_Esym_T(double* E, suNg_field* V) {
  int t,x,y,z,ix;
  int mu,nu;
  suNg_algebra_vector clover;
  double p;
  
  for(t=0;t<2*GLB_T;t++) E[t]=0.;

  for (t=0; t<T; t++){
    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++)
					      for(mu=0;mu<4;mu++) for(nu=mu+1;nu<4;nu++) {
						  ix=ipt(t,x,y,z);
						  WF_clover_F(&clover,V, ix, mu, nu);
						  _algebra_vector_sqnorm_g(p,clover);
						  if(mu==0) E[2*((zerocoord[0]+t+GLB_T)%GLB_T)] += p;
						  else E[2*((zerocoord[0]+t+GLB_T)%GLB_T)+1] += p;
						}
    E[2*((zerocoord[0]+t+GLB_T)%GLB_T)] *= _FUND_NORM2/(6.*GLB_VOL3);
    E[2*((zerocoord[0]+t+GLB_T)%GLB_T)+1] *= _FUND_NORM2/(6.*GLB_VOL3);
  }

  global_sum(E,2*GLB_T);
}


/*
  q = 1/(16 \pi^2) \epsilon_{\mu\nu\rho\sigma} \tr F_{\mu\nu} F_{\rho\sigma}
*/
double WF_topo(suNg_field* V) {
  double TC=0.;
  
  _MASTER_FOR_SUM(&glattice,ix,TC) {
    suNg_algebra_vector F1, F2;
    WF_clover_F(&F1,V,ix,1,2);
    WF_clover_F(&F2,V,ix,0,3);
    for(int i=0;i<NG*NG-1;i++) TC += F1.c[i]*F2.c[i];
    
    WF_clover_F(&F1,V,ix,1,3);
    WF_clover_F(&F2,V,ix,0,2);
    for(int i=0;i<NG*NG-1;i++) TC -= F1.c[i]*F2.c[i];
    
    WF_clover_F(&F1,V,ix,0,1);
    WF_clover_F(&F2,V,ix,2,3);
    for(int i=0;i<NG*NG-1;i++) TC += F1.c[i]*F2.c[i];
  }
  
  TC *= _FUND_NORM2/(4.*M_PI*M_PI);

  global_sum(&TC,1);
  
  return TC;
}

