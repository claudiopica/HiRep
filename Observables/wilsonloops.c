#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "update.h"


int multistep(int ix,int mu,int dx) {
  if(dx==0) return ix;
  else if(dx>0) return multistep(iup(ix,mu),mu,dx-1);
  else if(dx<0) return multistep(idn(ix,mu),mu,dx+1);
  return 0;
}

void wilsonloops(int mu, int x, suNg_field* g) {
  suNg_field* gfixed = alloc_gfield(&glattice);
  int i, y, nu;
  _DECLARE_INT_ITERATOR(ix);
  int ix1;
  suNg tmp, tmp1;
  int glb_size[4] = {GLB_T, GLB_X, GLB_Y, GLB_Z};
  complex wilson_fund;
  complex ctmp;
  double wilson_adj;

  error(g->type!=&glattice,1,"wilson.c","The gauge field must be defined on the whole lattice");
  error(x<=0,1,"wilson.c","The length must be > 0");
  error(mu<0 || mu>3,1,"wilson.c","Not defined direction");

  suNg_field_copy(gfixed,g);

  
  _MASTER_FOR(&glattice,ix) {
    ix1=ix;
    for(i=1;i<x;i++) {
      ix1=iup(ix1,mu);
      _suNg_times_suNg(tmp,*_4FIELD_AT(gfixed,ix,mu),*_4FIELD_AT(g,ix1,mu));
      *_4FIELD_AT(gfixed,ix,mu)=tmp;
    }
  }
  
  for(nu=0;nu<4;nu++) {
    if(nu==mu) continue;
    for(y=1;y<glb_size[nu];y++) {

      wilson_fund.re=wilson_fund.im=0.0;
      wilson_adj=0.0;
      _MASTER_FOR(&glattice,ix) {
        ix1=multistep(ix,nu,y);
        _suNg_times_suNg(tmp,*_4FIELD_AT(gfixed,ix,nu),*_4FIELD_AT(gfixed,ix1,mu));
        ix1=multistep(ix,mu,x);
        _suNg_times_suNg_dagger(tmp1,tmp,*_4FIELD_AT(gfixed,ix1,nu));
        _suNg_times_suNg_dagger(tmp,tmp1,*_4FIELD_AT(gfixed,ix,mu));
        _suNg_trace_re(ctmp.re,tmp);
        _suNg_trace_im(ctmp.im,tmp);
        _complex_add_assign(wilson_fund,ctmp);
        wilson_adj += _complex_prod_re(ctmp,ctmp)-1.;
      }
      _complex_mulr(wilson_fund,1./(GLB_T*GLB_X*GLB_Y*GLB_Z*NG),wilson_fund);
      wilson_adj /= GLB_T*GLB_X*GLB_Y*GLB_Z*(NG*NG-1);
      
      lprintf("WILSON LOOPS",0,"(mu,x,nu,y,f.re,f.im,adj) = %d %d %d %d %.8e %.8e %.8e\n",mu,x,nu,y,wilson_fund.re,wilson_fund.im,wilson_adj);
      
      if(y==glb_size[nu]-1) break;

      _MASTER_FOR(&glattice,ix) {
        ix1=multistep(ix,nu,y);
        _suNg_times_suNg(tmp,*_4FIELD_AT(gfixed,ix,nu),*_4FIELD_AT(g,ix1,nu));
        *_4FIELD_AT(gfixed,ix,nu)=tmp;
      }

    }
  }
  
  free_gfield(gfixed);
}

