#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "update.h"
#include "utils.h"
#include <math.h>


#ifdef WITH_MPI
#error Please compile without MPI!
#endif


#if NG==2
static void project_on_suN(suNg *A) {
  project_cooling_to_suNg(A,A,0);
}
#else
#error ERROR: The cooling parameter must be chosen!!!
#endif

static int index2[4][4]={
{-1, 0, 1, 2},
{ 0,-1, 3, 4},
{ 1, 3,-1, 5},
{ 2, 4, 5,-1}
};


void HYP_smearing(suNg_field* out, suNg_field* in, double weight[3]) {
  _DECLARE_INT_ITERATOR(ix);
  int mu,nu,rho,eta;
  int iy;
  int i;
  suNg tmp[3];
  suNg_field* Vbar[6];
  suNg_field* Vtilde[4];

  error(out->type!=&glattice,1,"HYP_smearing.c","'out' in HYP_core must be defined on the whole lattice");
  error(in->type!=&glattice,1,"HYP_smearing.c","'in' in HYP_core must be defined on the whole lattice");

  lprintf("HYP",30,"HYP smearing with weights %f %f %f\n",weight[0],weight[1],weight[2]);

  for(i=0;i<6;i++) Vbar[i]=alloc_gfield(&glattice);
  for(i=0;i<4;i++) Vtilde[i]=alloc_gfield(&glattice);


  for(nu=0;nu<4;nu++)
  for(rho=nu+1;rho<4;rho++) {
    i=index2[nu][rho];
    
    _MASTER_FOR(&glattice,ix) {
      for(mu=0;mu<4;mu++) {
        if(mu==nu || mu==rho) continue;

        _suNg_zero(tmp[0]);
        for(eta=0;eta<4;eta++){
          if(eta==mu || eta==nu || eta==rho) continue;
          
          iy=iup(ix,eta);
          _suNg_times_suNg(tmp[1],*_4FIELD_AT(in,ix,eta),*_4FIELD_AT(in,iy,mu));
          iy=iup(ix,mu);
          _suNg_times_suNg_dagger(tmp[2],tmp[1],*_4FIELD_AT(in,iy,eta));
          _suNg_add_assign(tmp[0],tmp[2]);
          
          iy=idn(ix,eta);
          _suNg_dagger(tmp[2],*_4FIELD_AT(in,iy,eta));
          _suNg_times_suNg(tmp[1],tmp[2],*_4FIELD_AT(in,iy,mu));
          iy=iup(iy,mu);
          _suNg_times_suNg(tmp[2],tmp[1],*_4FIELD_AT(in,iy,eta));
          _suNg_add_assign(tmp[0],tmp[2]);
        }

        *_4FIELD_AT(Vbar[i],ix,mu)=*_4FIELD_AT(in,ix,mu);
        if(mu!=nu && mu!=rho) {
          _suNg_mul(*_4FIELD_AT(Vbar[i],ix,mu),(1.-weight[2]),*_4FIELD_AT(Vbar[i],ix,mu));
          _suNg_mul(tmp[0],weight[2]/2.,tmp[0]);
          _suNg_add_assign(*_4FIELD_AT(Vbar[i],ix,mu),tmp[0]);
          
          project_on_suN(_4FIELD_AT(Vbar[i],ix,mu));
        }
          
      }
    }
    
  }


  for(nu=0;nu<4;nu++) {

    _MASTER_FOR(&glattice,ix) {
      for(mu=0;mu<4;mu++) {
        if(mu==nu) continue;

        _suNg_zero(tmp[0]);
        for(rho=0;rho<4;rho++){
          if(rho==mu || rho==nu) continue;
          
          iy=iup(ix,rho);
          _suNg_times_suNg(tmp[1],*_4FIELD_AT(Vbar[index2[nu][mu]],ix,rho),*_4FIELD_AT(Vbar[index2[rho][nu]],iy,mu));
          iy=iup(ix,mu);
          _suNg_times_suNg_dagger(tmp[2],tmp[1],*_4FIELD_AT(Vbar[index2[nu][mu]],iy,rho));
          _suNg_add_assign(tmp[0],tmp[2]);
          
          iy=idn(ix,rho);
          _suNg_dagger(tmp[2],*_4FIELD_AT(Vbar[index2[nu][mu]],iy,rho));
          _suNg_times_suNg(tmp[1],tmp[2],*_4FIELD_AT(Vbar[index2[rho][nu]],iy,mu));
          iy=iup(iy,mu);
          _suNg_times_suNg(tmp[2],tmp[1],*_4FIELD_AT(Vbar[index2[nu][mu]],iy,rho));
          _suNg_add_assign(tmp[0],tmp[2]);
        }

        *_4FIELD_AT(Vtilde[nu],ix,mu)=*_4FIELD_AT(in,ix,mu);
        if(mu!=nu) {
          _suNg_mul(*_4FIELD_AT(Vtilde[nu],ix,mu),(1.-weight[1]),*_4FIELD_AT(Vtilde[nu],ix,mu));
          _suNg_mul(tmp[0],weight[1]/4.,tmp[0]);
          _suNg_add_assign(*_4FIELD_AT(Vtilde[nu],ix,mu),tmp[0]);
          
          project_on_suN(_4FIELD_AT(Vtilde[nu],ix,mu));
        }
          
      }
    }
  
  }


  _MASTER_FOR(&glattice,ix) {
    for(mu=0;mu<4;mu++) {

      _suNg_zero(tmp[0]);
      for(nu=0;nu<4;nu++){
        if(nu==mu) continue;
        
        iy=iup(ix,nu);
        _suNg_times_suNg(tmp[1],*_4FIELD_AT(Vtilde[mu],ix,nu),*_4FIELD_AT(Vtilde[nu],iy,mu));
        iy=iup(ix,mu);
        _suNg_times_suNg_dagger(tmp[2],tmp[1],*_4FIELD_AT(Vtilde[mu],iy,nu));
        _suNg_add_assign(tmp[0],tmp[2]);
        
        iy=idn(ix,nu);
        _suNg_dagger(tmp[2],*_4FIELD_AT(Vtilde[mu],iy,nu));
        _suNg_times_suNg(tmp[1],tmp[2],*_4FIELD_AT(Vtilde[nu],iy,mu));
        iy=iup(iy,mu);
        _suNg_times_suNg(tmp[2],tmp[1],*_4FIELD_AT(Vtilde[mu],iy,nu));
        _suNg_add_assign(tmp[0],tmp[2]);
      }

      *_4FIELD_AT(out,ix,mu)=*_4FIELD_AT(in,ix,mu);
      _suNg_mul(*_4FIELD_AT(out,ix,mu),(1.-weight[0]),*_4FIELD_AT(out,ix,mu));
      _suNg_mul(tmp[0],weight[0]/6.,tmp[0]);
      _suNg_add_assign(*_4FIELD_AT(out,ix,mu),tmp[0]);
      
      project_on_suN(_4FIELD_AT(out,ix,mu));
        
    }
  }


  for(i=0;i<6;i++) free_gfield(Vbar[i]);
  for(i=0;i<4;i++) free_gfield(Vtilde[i]);

}


double min_tplaq(suNg_field* g) {
  _DECLARE_INT_ITERATOR(ix);
  int mu;
  int iy,iz;
  double ret=2.,p;
  suNg *v1,*v2,*v3,*v4,w1,w2,w3;

  _MASTER_FOR(&glattice,ix) {
    iy=iup(ix,0);
    v1=_4FIELD_AT(g,ix,0);
    for(mu=1;mu<4;mu++) {
      iz=iup(ix,mu);
      v2=_4FIELD_AT(g,iy,mu);
      v3=_4FIELD_AT(g,iz,0);
      v4=_4FIELD_AT(g,ix,mu);
      _suNg_times_suNg(w1,(*v1),(*v2));
      _suNg_times_suNg(w2,(*v4),(*v3));
      _suNg_times_suNg_dagger(w3,w1,w2);      
      _suNg_trace_re(p,w3);
    }
    if(p<ret) ret=p;
  }
  
  return ret;
}


void HYP_span_parameters(double mtp[6859]) {
  int i,j,k;
  double w[3];
  suNg_field* sg;

  sg=alloc_gfield(&glattice);
  
  for(i=0;i<19;i++) {
    w[0]=.05*(i+1);
    for(j=0;j<19;j++) {
      w[1]=.05*(j+1);
      for(k=0;k<19;k++) {
        w[2]=.05*(k+1);
        
        HYP_smearing(sg,u_gauge,w);
        mtp[i+19*(j+19*k)]+=min_tplaq(sg);
      }
    }
  }
}


int HYP_best_parameters(double mtp[6859], double w[3]) {
  int i, ibest;
  double max;
  
  max=mtp[0];
  w[0]=w[1]=w[2]=.05;
  ibest=0;
  for(i=1;i<6859;i++) {
    if(mtp[i]>max) {
      max=mtp[i];
      w[0]=.05*(i%19 + 1);
      w[1]=.05*((i/19)%19 + 1);
      w[2]=.05*(i/361 + 1);
      ibest=i;
    }
  }
  
  return ibest;
}

