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
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "propagator.h"
#include "gaugefix.h"

static void fix_T_bc(int tau){
  int index;
  int ix,iy,iz;
  suNf *u;
  if (--tau<0) tau+= GLB_T;
  lprintf("meson_measurements",15,"Setting Dirichlet boundary conidtion at global time slice %d, %d\n",tau,T_BORDER);
  if((zerocoord[0]-1<=tau && zerocoord[0]+T>tau) || (zerocoord[0]==0 && tau==GLB_T-1)) { 
    for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
	  if( ( (tau==zerocoord[0]-1) || (zerocoord[0]==0 && tau==GLB_T-1)) && (NP_T>1) ){
	    index=ipt_ext(0,ix,iy,iz);
	  }
	  else{
	    index=ipt_ext(T_BORDER+tau-zerocoord[0],ix,iy,iz); 
	  }
	  if(index!=-1) {
	    u=pu_gauge_f(index,0);
	    _suNf_zero(*u);
	  }
	}
  }
  lprintf("meson_measurements",50,"Boundaries set!\n");
}

static void flip_T_bc(int tau){
  int index;
  int ix,iy,iz;
  suNf *u;
  tau-=1;
  if (tau<0) tau+= GLB_T;
  lprintf("meson_measurements",15,"Flipping the boundary at global time slice %d\n",tau);
  fflush(stdout);
  if((zerocoord[0]-1<=tau && zerocoord[0]+T>tau) || (zerocoord[0]==0 && tau==GLB_T-1)) { 
    for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
	  if ((tau==zerocoord[0]-1) || (zerocoord[0]==0 && tau==GLB_T-1)){
	    index=ipt_ext(0,ix,iy,iz);
	  }
	  else{
	    index=ipt_ext(T_BORDER+tau-zerocoord[0],ix,iy,iz);
	  }
	  if(index!=-1) {
	    u=pu_gauge_f(index,0);
	    _suNf_minus(*u,*u);
	  }
	}
  }
  lprintf("meson_measurements",50,"Flipping DONE!\n");
}


void measure_spectrum_semwall(int nm, double* m, int nhits,int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice); //This isn't glat_even so that the odd sites will be set to zero explicitly
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  int tau,k,beta;
  init_propagator_eo(nm, m, precision);
  for (k=0;k<nhits;++k){
    tau=create_diluted_source_equal_eo(source);
    for(beta=0;beta<4;beta++) source[beta].type = &glat_even;
    calc_propagator_eo(prop,source,4);//4 for spin dilution
    for(beta=0;beta<4;beta++) source[beta].type = &glattice;
    measure_mesons(prop,source,nm,tau);
  }
  print_mesons(nhits*GLB_VOL3/2.,conf_num,nm,m,"DEFAULT_SEMWALL");
  free_propagator_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
}

void measure_spectrum_discon_semwall(int nm, double* m, int nhits,int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  int tau,k,beta;
  tau = 0;
  init_propagator_eo(nm, m, precision);
  for (k=0;k<nhits;++k){
      create_noise_source_equal_eo(source);
      for(beta=0;beta<4;beta++) source[beta].type = &glat_even;
      calc_propagator(prop,source,4);//4 for spin dilution
      for(beta=0;beta<4;beta++) source[beta].type = &glattice;
      measure_discon(prop,source,nm,tau);
  }
  print_mesons(nhits*GLB_VOL3/2.,conf_num,nm,m,"DISCON_SEMWALL");
  free_propagator_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
}

void measure_spectrum_discon_gfwall(int nm, double* m, int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  suNg_field* u_gauge_old=alloc_gfield(&glattice);
  int tau,k;
  tau = 0;

  suNg_field_copy(u_gauge_old,u_gauge);
  //Fix the Gauge
  double act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
	1.8,	//overrelax
	10000,	//maxit
	1e-12, //tolerance
	u_gauge //gauge
	);
  lprintf("GFWALL",0,"Gauge fixed action  %1.6f\n",act);
  double p2 = calc_plaq(u_gauge);
  lprintf("TEST",0,"fixed_gauge plaq %1.6f\n",p2);
    full_plaquette();
  represent_gauge_field();
  init_propagator_eo(nm, m, precision);

  for (tau=0;tau<GLB_T;++tau){
  for (k=0;k<NF;++k){
      create_gauge_fixed_wall_source(source, tau, k);
      calc_propagator(prop,source,4);//4 for spin dilution
      create_point_source(source,tau,k); //to get the contraction right
      measure_discon(prop,source,nm,tau);

  }}
  //This gets the norm of the 2pt wrong by a factor GLB_VOL3 but the norm of the disconnected right
  print_mesons(GLB_T,conf_num,nm,m,"DISCON_GFWALL"); 

  suNg_field_copy(u_gauge,u_gauge_old);
  represent_gauge_field();

  free_propagator_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
  free_gfield(u_gauge_old);
}

void measure_spectrum_discon_volume(int nm, double* m, int conf_num, double precision, int dil){

  //Spin diluted
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);

  init_propagator_eo(nm, m, precision);
  int p;
  for(p=0;p<dil;p++){
  	create_diluted_volume_source(source, p, dil);
  	calc_propagator(prop,source,4);//spin dilution
  	measure_discon(prop,source,nm,0);
  }

  print_mesons(1.,conf_num,nm,m,"DISCON_VOLUME"); 

  free_propagator_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);

}


void measure_spectrum_gfwall(int nm, double* m, int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  suNg_field* u_gauge_old=alloc_gfield(&glattice);

  int tau,k;
  tau = 0;
  suNg_field_copy(u_gauge_old,u_gauge);
  //Fix the Gauge
  double act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
	1.8,	//overrelax
	10000,	//maxit
	1e-12, //tolerance
	u_gauge //gauge
	);
  lprintf("GFWALL",0,"Gauge fixed action  %1.6f\n",act);
  double p2 = calc_plaq(u_gauge);
  lprintf("TEST",0,"fixed_gauge plaq %1.6f\n",p2);
    full_plaquette();
  represent_gauge_field();

  init_propagator_eo(nm, m, precision);
  for (k=0;k<NF;++k){
      create_gauge_fixed_wall_source(source, tau, k);
      calc_propagator(prop,source,4);//4 for spin dilution
      measure_mesons(prop, source, nm, tau);
  }
  print_mesons(GLB_VOL3,conf_num,nm,m,"DEFAULT_GFWALL");

  suNg_field_copy(u_gauge,u_gauge_old);
  represent_gauge_field();

  free_propagator_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
  free_gfield(u_gauge_old);
}


static void generate_random_point(int *pr) {

  double ran;
  ranlxd(&ran,1);
  pr[1] = (int)(ran*GLB_X);
  ranlxd(&ran,1);
  pr[2] = (int)(ran*GLB_Y);
  ranlxd(&ran,1);
  pr[3] = (int)(ran*GLB_Z);
  ranlxd(&ran,1);
  pr[0] = (int)(ran*GLB_T);

  bcast_int(pr,4);

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



void measure_spectrum_semwall_ext(int nm, double* m, int nhits,int conf_num, double precision){
  int k,l,tau;
  spinor_field* source = alloc_spinor_field_f(4,&glat_even);
  spinor_field* prop_p = alloc_spinor_field_f(8*nm,&glattice);
  spinor_field* prop_a = prop_p+4*nm;
  init_propagator_eo(nm, m, precision);
  for (k=0;k<nhits;++k){	
    tau=create_diluted_source_equal_eo(source);
    lprintf("MEASURE_SPECTRUM_SEMWALL_EXT",10,"SEM wall source (noise) at time slice %d.\n",tau);
    calc_propagator_eo(prop_p,source,4);//4 for spin dilution
    measure_mesons_ext(prop_p,source,nm,tau,0);
    flip_T_bc(tau);
    calc_propagator_eo(prop_a,source,4);//4 for spin dilution	
    flip_T_bc(tau);
    for (l=0;l<4*nm;++l) {
      spinor_field_add_assign_f(&prop_p[l],&prop_a[l]);
      spinor_field_mul_f(&prop_p[l],0.5,&prop_p[l]);
    }
    measure_mesons_ext(prop_p,source,nm,tau,1);
    for (l=0;l<4*nm;++l) spinor_field_sub_assign_f(&prop_p[l],&prop_a[l]);
    measure_mesons_ext(prop_p,source,nm,tau,2);
  }
  print_mesons_ext(1.*nhits*GLB_VOL3/2.,conf_num,nm,m,"EXTENDED_SEMWALL");
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
    if (n_mom>0){
      measure_point_mesons_momenta_ext(prop_p, source, nm, tau, n_mom,0);
    }
    else{
      measure_mesons_ext(prop_p,source,nm,tau,0);
    }
    flip_T_bc(tau);
    calc_propagator(prop_a,source,4);//4 for spin components
    flip_T_bc(tau);
    for (l=0;l<4*nm;++l){
      spinor_field_add_assign_f(&prop_p[l],&prop_a[l]);
      spinor_field_mul_f(&prop_p[l],0.5,&prop_p[l]);
    }
    if (n_mom>0){
      measure_point_mesons_momenta_ext(prop_p, source, nm, tau, n_mom,1);
      for (l=0;l<4*nm;++l) spinor_field_mul_add_assign_f(&prop_p[l],-1.,&prop_a[l]);
      measure_point_mesons_momenta_ext(prop_p, source, nm, tau, n_mom,2);	  
    }
    else{
      measure_mesons_ext(prop_p, source, nm, tau, 1);
      for (l=0;l<4*nm;++l) spinor_field_mul_add_assign_f(&prop_p[l],-1.,&prop_a[l]);
      measure_mesons_ext(prop_p, source, nm, tau, 2);
    }
  }
  if (n_mom>0){
    print_mesons_momenta_ext(conf_num,nm,m,n_mom,"EXTENDED_POINT");
  }
  else{
    print_mesons_ext(1.,conf_num,nm,m,"EXTENDED_POINT");
  }
  free_propagator_eo();
  free_spinor_field_f(source);
  free_spinor_field_f(prop_p);
  free_spinor_field_f(prop_a);
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
  fix_T_bc(tau-dt);//Apply fixed boundaryconditions by zeroing links at time slice tau to direction 0.
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

void measure_spectrum_gfwall_fixedbc(int dt, int nm, double* m, int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  suNg_field* u_gauge_old=alloc_gfield(&glattice);

  int tau,k;
  tau = 0;
  suNg_field_copy(u_gauge_old,u_gauge);

  //Fix the Gauge
  double act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
	1.8,	//overrelax
	10000,	//maxit
	1e-12, //tolerance
	u_gauge //gauge
	);
  lprintf("GFWALL",0,"Gauge fixed action  %1.6f\n",act);
  double p2 = calc_plaq(u_gauge);
  lprintf("TEST",0,"fixed_gauge plaq %1.6f\n",p2);
    full_plaquette();
  represent_gauge_field();

  fix_T_bc(tau-dt);//Apply fixed boundaryconditions by zeroing links at time slice tau to direction 0.

  init_propagator_eo(nm, m, precision);
  for (k=0;k<NF;++k){
      create_gauge_fixed_wall_source(source, tau, k);
      calc_propagator(prop,source,4);//4 for spin dilution
      measure_mesons(prop, source, nm, tau);
  }
  print_mesons(GLB_VOL3,conf_num,nm,m,"DIRICHLET_GFWALL");

  suNg_field_copy(u_gauge,u_gauge_old);
  represent_gauge_field();

  free_propagator_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
  free_gfield(u_gauge_old);
}

void measure_formfactor_pt(int ti, int tf, int nm, double* m, int n_mom, int conf_num, double precision){
  spinor_field* source;
  spinor_field* source_seq;
  spinor_field* prop_i;
  spinor_field* prop_seq;
  int k;

  source = alloc_spinor_field_f(4,&glattice);
  prop_i = alloc_spinor_field_f(4*NF,&glattice);
  prop_seq = alloc_spinor_field_f(4*NF,&glattice);
  source_seq = alloc_spinor_field_f(4*NF,&glattice);

  init_propagator_eo(1, m, precision);//1 for number of masses 

  for (k=0;k<NF;++k){
    create_point_source(source,ti,k);
    calc_propagator(prop_i+4*k,source,4);//4 for spin components
  }
    create_sequential_source(source_seq,tf,prop_i);
    calc_propagator(prop_seq,source_seq,4*NF);

  measure_formfactors(prop_seq, prop_i, source_seq, nm, ti, tf, n_mom); //eats two propagators
  print_formfactor(conf_num,nm,m,n_mom, "DEFAULT_FF_POINT",tf-ti);
  free_spinor_field_f(source);
  free_spinor_field_f(source_seq);
  free_spinor_field_f(prop_i);
  free_spinor_field_f(prop_seq);
  free_propagator_eo(); 
}

void measure_formfactor_fixed(int ti, int tf, int dt, int nm, double* m, int n_mom, int conf_num, double precision){
  spinor_field* source;
  spinor_field* prop_i;
  spinor_field* source_seq;
  spinor_field* prop_seq;

  int k;
  char label[100];
  suNf_field* u_gauge_old=alloc_gfield_f(&glattice);
  suNf_field_copy(u_gauge_old,u_gauge_f);//Save the gaugefield

  double p = avr_plaquette();
  lprintf("MESON_MEASUREMENTS",0,"<P> = %g\n", p);
  fix_T_bc(ti-dt);//Apply fixed boundaryconditions by zeroing links at time slice tau to direction 0.
  p = avr_plaquette();
  lprintf("MESON_MEASUREMENTS",0,"<P> = %g\n", p);

  source = alloc_spinor_field_f(4,&glattice);
  source_seq = alloc_spinor_field_f(4*NF,&glattice);
  prop_i = alloc_spinor_field_f(4*NF,&glattice);
  prop_seq = alloc_spinor_field_f(4*NF,&glattice);

  init_propagator_eo(1, m, precision);//1 for number of masses 
  for (k=0;k<NF;++k){
    create_point_source(source,ti,k);
    calc_propagator(prop_i+4*k,source,4);//4 for spin components
  }
  create_sequential_source(source_seq,tf,prop_i); //prop_i = S(x,0);
  calc_propagator(prop_seq,source_seq,4*NF); //prop_seq = S(y,x) S(x,0) delta(x, (2,0,0,0) )

  measure_formfactors(prop_seq, prop_i, source, nm, ti, tf, n_mom); //eats two propagators

  sprintf(label,"DIRICHLET_FF_POINT dt=%d",dt);
  print_formfactor(conf_num,nm,m,n_mom, label,tf-ti);
  suNf_field_copy(u_gauge_f,u_gauge_old);//Restore the gaugefield

  free_spinor_field_f(source);  
  free_spinor_field_f(source_seq);
  free_spinor_field_f(prop_i);
  free_spinor_field_f(prop_seq);

  free_gfield_f(u_gauge_old);
  free_propagator_eo(); 
}

