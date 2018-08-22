/*******************************************************************************
*                                                                              *
* Wrapper functions for different correlator measurements                      *
* for four fermion interaction                                                 *
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen, Jarno Rantaharju               *
*                                                                              *
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
#include "gaugefix.h"
#include "meson_observables.h"

//The flavors are not degenerate, the pi field has opposite sign
//We can to measure both separately and add together

//This changes the sign of the pi field
void flip_scalar_field(scalar_field *f);

void create_noise_source_eo(spinor_field *source);
void calc_propagator_hopping_series(spinor_field *psi, spinor_field *eta, int hopping, int ndilute);
void calc_propagator_hopping_oe(spinor_field *psi, spinor_field *eta, int hopping, int ndilute);
void create_noise_source_equal_oe(spinor_field *source);


//We often need corr_u-corr_d,
//This adds the minus-sign
static void flip_corrs(meson_observable* mo){
  meson_observable* motmp=mo;
  int i;
  while (motmp!=NULL){
    for(i=0; i<motmp->corr_size; i++){
      motmp->corr_re[i] = -motmp->corr_re[i];
      motmp->corr_im[i] = -motmp->corr_im[i];
    }
    motmp=motmp->next;
  }
}

static void fix_T_bc(int tau){
  int index;
  int ix,iy,iz;
  suNf *u;
  if (--tau<0) tau+= GLB_T;
  lprintf("meson_measurements",15,"Setting Dirichlet boundary conidtion at global time slice %d (local %d)\n",tau,tau-zerocoord[0]);
  if((zerocoord[0]-1<=tau && zerocoord[0]+T>tau) || (zerocoord[0]==0 && tau==GLB_T-1)) { 
    for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
	  if (tau==zerocoord[0]-1 || (zerocoord[0]==0 && tau==GLB_T-1)){
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


/********************************
*	Point Sources		*
*********************************/

#define corr_ind(px,py,pz,n_mom,tc,nm,cm) ((px)*(n_mom)*(n_mom)*(24)*(nm)+(py)*(n_mom)*(24)*(nm)+(pz)*(24)*(nm)+ ((cm)*(24)) +(tc))
void measure_spectrum_ff_pt(int tau, int nm, double* m, int n_mom,int nhits,int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm*NF,&glattice);
  init_propagator_ff_eo(nm, m, precision);
  int k;
  lprintf("MAIN",0,"Point Source at (%d,0,0,0) \n",tau);
  for (k=0;k<NF;++k){
    create_point_source(source,tau,k);
    calc_propagator_ff_eo(prop + 4*k,source,4);//4 for spin components
    if (n_mom>1){
      measure_point_mesons_momenta(meson_correlators,prop+4*k, source, nm, tau, n_mom);
    }
    else{
      measure_mesons(meson_correlators,prop+4*k, source, nm, tau);
    }
  }

  flip_scalar_field(ff_pi);
  
  for (k=0;k<NF;++k){
    create_point_source(source,tau,k);
    calc_propagator_ff_eo(prop + 4*k,source,4);//4 for spin components
    if (n_mom>1){
      measure_point_mesons_momenta(meson_correlators,prop+4*k, source, nm, tau, n_mom);
    }
    else{
      measure_mesons(meson_correlators,prop+4*k, source, nm, tau);
    }
  }
  print_mesons(meson_correlators,1.,conf_num,nm,m,GLB_T,n_mom,"DEFAULT_POINT");
  flip_scalar_field(ff_pi);

 
  free_propagator_ff_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
}


/********************************
*	SEMWall Sources		*
*********************************/

void measure_spectrum_ff_semwall(int nm, double* m, int nhits,int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  int tau,k;
  init_propagator_ff_eo(nm, m, precision);
  for (k=0;k<nhits;++k){
    tau=create_diluted_source_equal_eo(source);
    calc_propagator_ff_eo(prop,source,4);//4 for spin dilution
    measure_mesons(meson_correlators,prop,source,nm,tau);

    //Change the sign of the pi field
    flip_scalar_field(ff_pi);
    //create_diluted_source_equal_atau_eo(source,tau);
    calc_propagator_ff_eo(prop,source,4);//4 for spin dilution
    measure_mesons(meson_correlators,prop,source,nm,tau);
    flip_scalar_field(ff_pi); //back to normal
  }
  print_mesons(meson_correlators,nhits*GLB_VOL3/2.,conf_num,nm,m,GLB_T,1,"DEFAULT_SEMWALL");


  free_propagator_ff_eo();
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
}


void measure_spectrum_semwall_ff_ext(int nm, double* m, int nhits,int conf_num, double precision){
  int k,l,tau;
  spinor_field* source = alloc_spinor_field_f(4,&glat_even);
  spinor_field* prop_p = alloc_spinor_field_f(8*nm,&glattice);
  spinor_field* prop_a = prop_p+4*nm;

  int dilution = 4; //4 for spin dilution
  init_propagator_ff_eo(nm, m, precision);
  for (k=0;k<nhits;++k){
    tau=create_diluted_source_equal_eo(source);
    lprintf("MEASURE_SPECTRUM_SEMWALL_EXT",10,"SEM wall source (noise) at time slice %d.\n",tau);
    calc_propagator_ff_eo(prop_p,source,dilution);
    measure_mesons_ext(meson_correlators,prop_p,source,nm,tau,0);
    flip_T_bc(tau);
    calc_propagator_ff_eo(prop_a,source,dilution);
    flip_T_bc(tau);
    for (l=0;l<dilution*nm;++l) {
      spinor_field_add_assign_f(&prop_p[l],&prop_a[l]);
      spinor_field_mul_f(&prop_p[l],0.5,&prop_p[l]);
    }
    measure_mesons_ext(meson_correlators,prop_p,source,nm,tau,1);
    for (l=0;l<dilution*nm;++l) spinor_field_sub_assign_f(&prop_p[l],&prop_a[l]);
    measure_mesons_ext(meson_correlators,prop_p,source,nm,tau,2);

    //Change the sign of the pi field
    flip_scalar_field(ff_pi);
    tau=create_diluted_source_equal_eo(source);
    lprintf("MEASURE_SPECTRUM_SEMWALL_EXT",10,"SEM wall source (noise) at time slice %d.\n",tau);
    calc_propagator_ff_eo(prop_p,source,dilution);
    measure_mesons_ext(meson_correlators,prop_p,source,nm,tau,0);
    flip_T_bc(tau);
    calc_propagator_ff_eo(prop_a,source,dilution);
    flip_T_bc(tau);
    for (l=0;l<dilution*nm;++l) {
      spinor_field_add_assign_f(&prop_p[l],&prop_a[l]);
      spinor_field_mul_f(&prop_p[l],0.5,&prop_p[l]);
    }
    measure_mesons_ext(meson_correlators,prop_p,source,nm,tau,1);
    for (l=0;l<dilution*nm;++l) spinor_field_sub_assign_f(&prop_p[l],&prop_a[l]);
    measure_mesons_ext(meson_correlators,prop_p,source,nm,tau,2);
    flip_scalar_field(ff_pi); //back to normal
  }
  print_mesons(meson_correlators,1.*nhits*GLB_VOL3,conf_num,nm,m,3*GLB_T,1,"EXTENDED_SEMWALL");
  free_propagator_ff_eo();
  free_spinor_field_f(source);
  free_spinor_field_f(prop_p);
}


/****************************************
*	Disconnected Measurements	*
*****************************************/

void measure_spectrum_discon_ff_semwall(int nm, double* m, int nhits, int degree_hopping, int nhits_hopping,int conf_num, double precision){
  spinor_field* source = alloc_spinor_field_f(4,&glattice);
  spinor_field* prop =  alloc_spinor_field_f(4*nm,&glattice);
  int k,beta;
  char label[100];
  init_propagator_ff_eo(nm, m, precision);

//hopping parameter expansion
// 1/D = sum(i=0 to k) (1/A*D)^i*1/A + (1/A*H)^k 1/D
//the first part, cheap
  for (k=0;k<nhits;++k){
    for(int kk=0;kk<nhits_hopping;++kk) {
      flip_corrs(triplet_discon_correlators); //For dd-uu
      create_noise_source_equal_eo(source);
      calc_propagator_ff_hopping_series(prop,source,degree_hopping,4); //4 for spin dilution
      measure_mesons(discon_correlators,prop,source,nm,0);
      measure_mesons(triplet_discon_correlators,prop,source,nm,0);
     
      //Change the sign of the pi field
      flip_scalar_field(ff_pi);
      flip_corrs(triplet_discon_correlators); //For dd-uu
      //create_noise_source_equal_eo(source);
      calc_propagator_ff_hopping_series(prop,source,degree_hopping,4); //4 for spin dilution
      measure_mesons(discon_correlators,prop,source,nm,0);
      measure_mesons(triplet_discon_correlators,prop,source,nm,0);
      sprintf(label,"src %d DISCON_SEMWALL_HOP",k);
      flip_scalar_field(ff_pi); //back to normal

      flip_corrs(triplet_discon_correlators); //For dd-uu
      create_noise_source_equal_oe(source);
      calc_propagator_ff_hopping_series(prop,source,degree_hopping,4); //4 for spin dilution
      measure_mesons(discon_correlators,prop,source,nm,0);
      measure_mesons(triplet_discon_correlators,prop,source,nm,0);
     
      //Change the sign of the pi field
      flip_scalar_field(ff_pi);
      flip_corrs(triplet_discon_correlators); //For dd-uu
      //create_noise_source_equal_eosource);
      calc_propagator_ff_hopping_series(prop,source,degree_hopping,4); //4 for spin dilution
      measure_mesons(discon_correlators,prop,source,nm,0);
      measure_mesons(triplet_discon_correlators,prop,source,nm,0);
      sprintf(label,"src %d DISCON_SEMWALL_HOP",k);
      flip_scalar_field(ff_pi); //back to normal
 
    }
    //Now printing the average of 10 sources for each hit
    if(nhits_hopping>0) print_mesons(discon_correlators,-(2.0*(double)nhits_hopping)/2.,conf_num,nm,m,GLB_T,1,label);
    if(nhits_hopping>0) print_mesons(triplet_discon_correlators,-(2.0*(double)nhits_hopping)/2.,conf_num,nm,m,GLB_T,1,label);
    
  }

//The second term in the hopping expansion,
  for (k=0;k<nhits;++k){
      create_noise_source_equal_eo(source);
      calc_propagator_ff_hopping_oe(prop,source,degree_hopping,4);//5 for order of hopping parameter expansion,  4 for spin dilution
      measure_mesons(triplet_discon_correlators,prop,source,nm,0);
      measure_mesons(discon_correlators,prop,source,nm,0);
     
      //Change the sign of the pi field
      flip_scalar_field(ff_pi);
      flip_corrs(triplet_discon_correlators); //For dd-uu
      //create_noise_source_equal_eo(source);
      calc_propagator_ff_hopping_oe(prop,source,degree_hopping,4);//5 for order of hopping parameter expansion,  4 for spin dilution
      measure_mesons(triplet_discon_correlators,prop,source,nm,0);
      measure_mesons(discon_correlators,prop,source,nm,0);
      sprintf(label,"src %d DISCON_SEMWALL",k);
      flip_scalar_field(ff_pi); //back to normal

      //print_mesons(triplet_discon_correlators,1./2.,conf_num,nm,m,GLB_T,1,label);

      flip_corrs(triplet_discon_correlators);
      create_noise_source_equal_oe(source);
      calc_propagator_ff_hopping_oe(prop,source,degree_hopping,4);//5 for order of hopping parameter expansion,  4 for spin dilution
      measure_mesons(triplet_discon_correlators,prop,source,nm,0);
      measure_mesons(discon_correlators,prop,source,nm,0);
     
      //Change the sign of the pi field
      flip_scalar_field(ff_pi);
      flip_corrs(triplet_discon_correlators); //For dd-uu
      calc_propagator_ff_hopping_oe(prop,source,degree_hopping,4);//5 for order of hopping parameter expansion,  4 for spin dilution
      measure_mesons(triplet_discon_correlators,prop,source,nm,0);
      measure_mesons(discon_correlators,prop,source,nm,0);
      sprintf(label,"src %d DISCON_SEMWALL",k);
      flip_scalar_field(ff_pi); //back to normal
      
      print_mesons(triplet_discon_correlators,2./2.,conf_num,nm,m,GLB_T,1,label);

  }
  //print_mesons(discon_correlators,nhits*GLB_VOL3/2.,conf_num,nm,m,GLB_T,1,"DISCON_SEMWALL");
  free_propagator_ff_eo(); 
  free_spinor_field_f(source);
  free_spinor_field_f(prop);
}





