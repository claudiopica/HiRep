/*******************************************************************************
*
* Wrapper functions for different type of measurements 
*
*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "spectrum.h"


void measure_spectrum_semwall(int nm, double* m, int nhits,int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glat_even);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  int tau,k;
  init_propagator_eo(nm, m, precision);
  for (k=0;k<nhits;++k){
    tau=create_diluted_source_equal_eo(source);
    calc_propagator_eo(prop,source,4);//4 for spin dilution
    measure_mesons(prop,source,nm,tau);
  }
  print_mesons(nhits*GLB_VOL3/2.,conf_num,nm,m,"DEFAULT_SEMWALL");
  free_propagator_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
}

void measure_spectrum_pt(int tau, int nm, double* m, int n_mom,int nhits,int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  init_propagator_eo(nm, m, precision);
  int k;
  lprintf("MAIN",0,"Point Source at (%d,0,0,0) \n",tau);
  for (k=0;k<NF;++k){
    create_point_source(source,tau,k);
    calc_propagator(prop,source,4);//4 for spin components
    if (n_mom>0){
      measure_point_mesons_momenta(prop, source, nm, tau, n_mom);
    }
    else{
      measure_mesons(prop, source, nm, tau);
    }
  }
  if (n_mom>0){
    print_mesons_momenta(conf_num,nm,m,n_mom,"DEFAULT_POINT");
  }
  else{
    print_mesons(1.,conf_num,nm,m,"DEFAULT_POINT");
  }
  free_propagator_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
}

static void flip_T_bc(int tau){
  int index;
  int ix,iy,iz;
  suNf *u;
  if (tau==0) tau+= GLB_T;
  if(COORD[0]==(tau-1)/T) {
    for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
      index=ipt(tau-1,ix,iy,iz);
      if(index!=-1) {
        u=pu_gauge_f(index,0);
        _suNf_minus(*u,*u);
      }
    }
  }
}


void measure_spectrum_semwall_ext(int nm, double* m, int nhits,int conf_num, double precision){
  int k,l,tau;
  spinor_field* source = alloc_spinor_field_f(4,&glat_even);
  spinor_field* prop_p = alloc_spinor_field_f(8*nm,&glattice);
  spinor_field* prop_a = prop_p+4*nm;
  init_propagator_eo(nm, m, precision);
  for (k=0;k<nhits;++k){	
    tau=create_diluted_source_equal_eo(source);
    calc_propagator_eo(prop_p,source,4);//4 for spin dilution
    flip_T_bc(tau);
    calc_propagator_eo(prop_a,source,4);//4 for spin dilution	
    flip_T_bc(tau);
    for (l=0;l<4*nm;++l) spinor_field_add_assign_f(&prop_p[l],&prop_a[l]);
    measure_mesons_ext(prop_p,source,nm,tau,1);
    for (l=0;l<4*nm;++l) spinor_field_mul_add_assign_f(&prop_p[l],-2.,&prop_a[l]);
    measure_mesons_ext(prop_p,source,nm,tau,0);
  }
  print_mesons_ext(4.*nhits*GLB_VOL3/2.,conf_num,nm,m,"EXTENDED_SEMWALL");
  free_propagator_eo();
  free_spinor_field_f(source);
  free_spinor_field_f(prop_p);
  free_spinor_field_f(prop_a);
}



void measure_spectrum_pt_ext(int tau, int nm, double* m, int n_mom,int nhits,int conf_num, double precision){
  int k,l;
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop_p = alloc_spinor_field_f(8*nm,&glattice);
  spinor_field* prop_a = prop_p+4*nm;
  init_propagator_eo(nm, m, precision);
  lprintf("MAIN",10,"Point Source at (%d,0,0,0) \n",tau);
  for (k=0;k<NF;++k){
    create_point_source(source,tau,k);
    calc_propagator(prop_p,source,4);//4 for spin components
    flip_T_bc(tau);
    calc_propagator(prop_a,source,4);//4 for spin components
    flip_T_bc(tau);
    for (l=0;l<4*nm;++l) spinor_field_add_assign_f(&prop_p[l],&prop_a[l]);
    if (n_mom>0){
      measure_point_mesons_momenta_ext(prop_p, source, nm, tau, n_mom,1);
      for (l=0;l<4*nm;++l) spinor_field_mul_add_assign_f(&prop_p[l],-2.,&prop_a[l]);
      measure_point_mesons_momenta_ext(prop_p, source, nm, tau, n_mom,0);	  
    }
    else{
      measure_mesons_ext(prop_p, source, nm, tau, 1);
      for (l=0;l<4*nm;++l) spinor_field_mul_add_assign_f(&prop_p[l],-2.,&prop_a[l]);
      measure_mesons_ext(prop_p, source, nm, tau, 0);
    }
  }
  if (n_mom>0){
    print_mesons_momenta_ext(conf_num,nm,m,n_mom,"EXTENDED_POINT");
  }
  else{
    print_mesons_ext(4.,conf_num,nm,m,"EXTENDED_POINT");
  }
  free_propagator_eo();
  free_spinor_field_f(source);
  free_spinor_field_f(prop_p);
  free_spinor_field_f(prop_a);
}

static void fix_T_bc(int tau){
  int index;
  int ix,iy,iz;
  suNf *u;
  if (tau<=0) tau+= GLB_T;
  if(COORD[0]==(tau-1)/T) {
    for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
      index=ipt(tau-1,ix,iy,iz);
      if(index!=-1) {
        u=pu_gauge_f(index,0);
        _suNf_zero(*u);
      }
    }
  }
}

void measure_spectrum_semwall_fixedbc(int dt, int nm, double* m, int nhits,int conf_num, double precision){
  int tau,k;
  spinor_field* source = alloc_spinor_field_f(4,&glat_even);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  suNf_field* u_gauge_old=alloc_gfield_f(&glattice);
  char label[100];
  suNf_field_copy(u_gauge_old,u_gauge_f);
  init_propagator_eo(nm, m, precision);
  for (k=0;k<nhits;++k){	
    tau=create_diluted_source_equal_eo(source);
    fix_T_bc(tau-dt);
    calc_propagator_eo(prop,source,4);//4 for spin dilution
    measure_mesons(prop,source,nm,tau);
    suNf_field_copy(u_gauge_f,u_gauge_old);
  }
  sprintf(label,"DIRICHLET_SEMWALL dt=%d",dt);
  print_mesons(nhits*GLB_VOL3/2.,conf_num,nm,m,label);
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
  free_gfield_f(u_gauge_old);
  free_propagator_eo(); 
}

void measure_spectrum_pt_fixedbc(int tau, int dt, int nm, double* m, int n_mom,int nhits,int conf_num, double precision){
  int k;
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  suNf_field* u_gauge_old=alloc_gfield_f(&glattice);
  char label[100];
  suNf_field_copy(u_gauge_old,u_gauge_f);
  init_propagator_eo(nm, m, precision);
  fix_T_bc(tau);
  lprintf("MAIN",0,"Point Source at (%d,0,0,0) \n",tau);
  for (k=0;k<NF;++k){
    create_point_source(source,tau,k);
    calc_propagator(prop,source,4);//4 for spin components
    if (n_mom>0){
      measure_point_mesons_momenta(prop, source, nm, tau, n_mom);
    }
    else{
      measure_mesons(prop, source, nm, tau);
    }
  }
  sprintf(label,"DIRICHLET_POINT dt=%d",dt);
  if (n_mom>0){
    print_mesons_momenta(conf_num,nm,m,n_mom,label);
  }
  else{
    print_mesons(1.,conf_num,nm,m,label);
  }
  suNf_field_copy(u_gauge_f,u_gauge_old);
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
  free_gfield_f(u_gauge_old);
  free_propagator_eo(); 
}
