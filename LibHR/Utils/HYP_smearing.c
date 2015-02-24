/* hep-lat/0103029 */

#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "update.h"
#include "utils.h"
#include "communications.h"
#include <math.h>
#include <stdio.h>


#if NG==2
static void project_on_suN(suNg *A) {
  project_cooling_to_suNg(A,A,0);
}
#else
//#error ERROR: The cooling parameter must be chosen!!!
static void project_on_suN(suNg *A) {
  error(1,1,"project_on_suN",
      "Error function only defined for NG=2");
}

#endif


static int two_to_index[4][4];
static int index_to_two[12][2];
static int three_to_index[4][4][4];
static int index_to_three[24][3];


static void init_HYP_indices() {
  static int init=(1==0);
  if(init) return;
  
  int i;
  
  i=0;
  for(int mu=0;mu<4;mu++) for(int nu=0;nu<4;nu++) {
    if(mu==nu) {
      two_to_index[mu][nu]=-1;
    } else {
      two_to_index[mu][nu]=i;
      index_to_two[i][0]=mu;
      index_to_two[i][1]=nu;
      i++;
    }
  }
  
/*  printf("two indices: %d\n",i);*/

  i=0;
  for(int mu=0;mu<4;mu++) for(int nu=0;nu<4;nu++) for(int rho=nu+1;rho<4;rho++) {
    if(mu==nu || mu==rho) {
      three_to_index[mu][nu][rho]=-1;
      three_to_index[mu][rho][nu]=-1;
    } else {
      three_to_index[mu][nu][rho]=i;
      three_to_index[mu][rho][nu]=i;
      index_to_three[i][0]=mu;
      index_to_three[i][1]=nu;
      index_to_three[i][2]=rho;
      i++;
    }
  }

/*  printf("three indices: %d\n",i);*/
  
  init=(1==1);
}


/*
void spatialHYP_smearing(suNg_field* out, suNg_field* in, double weight[3]) {
  _DECLARE_INT_ITERATOR(ix);
  int mu,nu,rho,eta;
  int iy;
  int i;
  suNg tmp[3];
  suNg_field* Vbar[6];
  suNg_field* Vtilde[4];

  error(out->type!=&glattice,1,"spatialHYP_smearing.c","'out' in HYP_core must be defined on the whole lattice");
  error(in->type!=&glattice,1,"spatialHYP_smearing.c","'in' in HYP_core must be defined on the whole lattice");

  lprintf("HYP",30,"Spatial HYP smearing with weights %f %f %f\n",weight[0],weight[1],weight[2]);

  for(i=0;i<6;i++) Vbar[i]=alloc_gfield(&glattice);
  for(i=0;i<4;i++) Vtilde[i]=alloc_gfield(&glattice);


  for(nu=1;nu<4;nu++)
  for(rho=nu+1;rho<4;rho++) {
    i=index2[nu][rho];
    
    _MASTER_FOR(&glattice,ix) {
      for(mu=1;mu<4;mu++) {
        if(mu==nu || mu==rho) continue;

        _suNg_zero(tmp[0]);
        for(eta=1;eta<4;eta++){
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


  for(nu=1;nu<4;nu++) {

    _MASTER_FOR(&glattice,ix) {
      for(mu=1;mu<4;mu++) {
        if(mu==nu) continue;

        _suNg_zero(tmp[0]);
        for(rho=1;rho<4;rho++){
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

      *_4FIELD_AT(out,ix,mu)=*_4FIELD_AT(in,ix,mu);
      if(mu==0) continue;

      _suNg_zero(tmp[0]);
      for(nu=1;nu<4;nu++){
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

      _suNg_mul(*_4FIELD_AT(out,ix,mu),(1.-weight[0]),*_4FIELD_AT(out,ix,mu));
      _suNg_mul(tmp[0],weight[0]/6.,tmp[0]);
      _suNg_add_assign(*_4FIELD_AT(out,ix,mu),tmp[0]);
      
      project_on_suN(_4FIELD_AT(out,ix,mu));
        
    }
  }


  for(i=0;i<6;i++) free_gfield(Vbar[i]);
  for(i=0;i<4;i++) free_gfield(Vtilde[i]);

}
*/


void HYP_smearing(suNg_field* out, suNg_field* in, double weight[3]) {
  suNg_field* Vbar[24];
  suNg_field* Vtilde[12];

  error(out->type!=&glattice,1,"HYP_smearing.c","'out' in HYP_core must be defined on the whole lattice");
  error(in->type!=&glattice,1,"HYP_smearing.c","'in' in HYP_core must be defined on the whole lattice");

  lprintf("HYP",30,"HYP smearing with weights %f %f %f\n",weight[0],weight[1],weight[2]);

  for(int i=0;i<24;i++) Vbar[i]=alloc_gtransf(&glattice);
  for(int i=0;i<12;i++) Vtilde[i]=alloc_gtransf(&glattice);

  init_HYP_indices();
  
/*

Vbar{x,mu;nu,rho} = Proj[
(1+a3) U{x,mu} +
a3/2 sum_{\pm eta \neq rho,nu,mu}
U{x,eta} U{x+eta,mu} U{x+mu,eta}^\dag
]

*/

  start_gf_sendrecv(in);

  _PIECE_FOR(&glattice,ixp) {
    if(ixp==glattice.inner_master_pieces) {
      _OMP_PRAGMA( master )
      complete_gf_sendrecv(in);
      _OMP_PRAGMA( barrier )
    }
    suNg tmp[3];
    _SITE_FOR(&glattice,ixp,ix) {

      for(int i=0;i<24;i++) {
        int mu=index_to_three[i][0];
        int nu=index_to_three[i][1];
        int rho=index_to_three[i][2];
        int eta=(0+1+2+3)-mu-nu-rho;
        int iy;
        
        _suNg_zero(tmp[0]);
        
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
  
        *_FIELD_AT(Vbar[i],ix)=*_4FIELD_AT(in,ix,mu);
        _suNg_mul(*_FIELD_AT(Vbar[i],ix),(1.-weight[2]),*_FIELD_AT(Vbar[i],ix));
        _suNg_mul(tmp[0],weight[2]/2.,tmp[0]);
        _suNg_add_assign(*_FIELD_AT(Vbar[i],ix),tmp[0]);
        
        project_on_suN(_FIELD_AT(Vbar[i],ix));
        
      }
    
    } /* SITE_FOR */
  } /* PIECE FOR */  

/*

Vtilde{x,mu;nu} = Proj[
(1+a2) U{x,mu} +
a2/4 sum_{\pm rho \neq nu,mu}
Vbar{x,rho;nu,mu} Vbar{x+rho,mu;rho,nu} Vbar{x+mu,rho;nu,mu}^\dag
]

*/

  for(int i=0;i<24;i++) start_gt_sendrecv(Vbar[i]);

  _PIECE_FOR(&glattice,ixp) {
    if(ixp==glattice.inner_master_pieces) {
      _OMP_PRAGMA( master )
      for(int i=0;i<24;i++) complete_gt_sendrecv(Vbar[i]);
      _OMP_PRAGMA( barrier )
    }
    suNg tmp[3];
    _SITE_FOR(&glattice,ixp, ix) {

      for(int i=0;i<12;i++) {
        int mu=index_to_two[i][0];
        int nu=index_to_two[i][1];
        
        _suNg_zero(tmp[0]);
        for(int rho=0;rho<4;rho++){
          int iy;
          if(rho==mu || rho==nu) continue;
          
          iy=iup(ix,rho);
          _suNg_times_suNg(tmp[1],*_FIELD_AT(Vbar[three_to_index[rho][nu][mu]],ix),*_FIELD_AT(Vbar[three_to_index[mu][rho][nu]],iy));
          iy=iup(ix,mu);
          _suNg_times_suNg_dagger(tmp[2],tmp[1],*_FIELD_AT(Vbar[three_to_index[rho][nu][mu]],iy));
          _suNg_add_assign(tmp[0],tmp[2]);
          
          iy=idn(ix,rho);
          _suNg_dagger(tmp[2],*_FIELD_AT(Vbar[three_to_index[rho][nu][mu]],iy));
          _suNg_times_suNg(tmp[1],tmp[2],*_FIELD_AT(Vbar[three_to_index[mu][rho][nu]],iy));
          iy=iup(iy,mu);
          _suNg_times_suNg(tmp[2],tmp[1],*_FIELD_AT(Vbar[three_to_index[rho][nu][mu]],iy));
          _suNg_add_assign(tmp[0],tmp[2]);
        }
  
        *_FIELD_AT(Vtilde[i],ix)=*_4FIELD_AT(in,ix,mu);
        _suNg_mul(*_FIELD_AT(Vtilde[i],ix),(1.-weight[1]),*_FIELD_AT(Vtilde[i],ix));
        _suNg_mul(tmp[0],weight[1]/4.,tmp[0]);
        _suNg_add_assign(*_FIELD_AT(Vtilde[i],ix),tmp[0]);
          
        project_on_suN(_FIELD_AT(Vtilde[i],ix));
          
      }
    } /* SITE_FOR */
  } /* PIECE FOR */  

/*

V{x,mu} = Proj[
(1+a1) U{x,mu} +
a1/6 sum_{\pm nu \neq mu}
Vtilde{x,nu;mu} Vtilde{x+nu,mu;nu} Vtilde{x+mu,nu;mu}^\dag
]

*/

  for(int i=0;i<12;i++) start_gt_sendrecv(Vtilde[i]);

  _PIECE_FOR(&glattice,ixp) {
    if(ixp==glattice.inner_master_pieces) {
      _OMP_PRAGMA( master )
      for(int i=0;i<12;i++) complete_gt_sendrecv(Vtilde[i]);
      _OMP_PRAGMA( barrier )
    }
    suNg tmp[3];
    _SITE_FOR(&glattice,ixp,ix) {

      for(int mu=0;mu<4;mu++) {

        _suNg_zero(tmp[0]);
        for(int nu=0;nu<4;nu++){
          int iy;
          if(nu==mu) continue;
          
          iy=iup(ix,nu);
          _suNg_times_suNg(tmp[1],*_FIELD_AT(Vtilde[two_to_index[nu][mu]],ix),*_FIELD_AT(Vtilde[two_to_index[mu][nu]],iy));
          iy=iup(ix,mu);
          _suNg_times_suNg_dagger(tmp[2],tmp[1],*_FIELD_AT(Vtilde[two_to_index[nu][mu]],iy));
          _suNg_add_assign(tmp[0],tmp[2]);
          
          iy=idn(ix,nu);
          _suNg_dagger(tmp[2],*_FIELD_AT(Vtilde[two_to_index[nu][mu]],iy));
          _suNg_times_suNg(tmp[1],tmp[2],*_FIELD_AT(Vtilde[two_to_index[mu][nu]],iy));
          iy=iup(iy,mu);
          _suNg_times_suNg(tmp[2],tmp[1],*_FIELD_AT(Vtilde[two_to_index[nu][mu]],iy));
          _suNg_add_assign(tmp[0],tmp[2]);
        }
  
        *_4FIELD_AT(out,ix,mu)=*_4FIELD_AT(in,ix,mu);
        _suNg_mul(*_4FIELD_AT(out,ix,mu),(1.-weight[0]),*_4FIELD_AT(out,ix,mu));
        _suNg_mul(tmp[0],weight[0]/6.,tmp[0]);
        _suNg_add_assign(*_4FIELD_AT(out,ix,mu),tmp[0]);
        
        project_on_suN(_4FIELD_AT(out,ix,mu));
        
      }
    } /* SITE_FOR */
  } /* PIECE FOR */  


  for(int i=0;i<24;i++) free_gtransf(Vbar[i]);
  for(int i=0;i<12;i++) free_gtransf(Vtilde[i]);

}


double min_tplaq(suNg_field* g) {

  double ret=2.;

  _MASTER_FOR_MIN(&glattice,ix,ret) {
    suNg w1,w2,w3;
    double p;
    int iy=iup(ix,0);
    suNg *v1=_4FIELD_AT(g,ix,0);
    for(int mu=1;mu<4;mu++) {
      int iz=iup(ix,mu);
      suNg *v2=_4FIELD_AT(g,iy,mu);
      suNg *v3=_4FIELD_AT(g,iz,0);
      suNg *v4=_4FIELD_AT(g,ix,mu);
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

