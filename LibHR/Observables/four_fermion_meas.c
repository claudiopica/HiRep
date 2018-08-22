/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File four_fermion_meas.c
*
* Averages of the auxiliary fields
*
*******************************************************************************/

#include "global.h"
#include "suN.h"
#include "communications.h"
#include "logger.h"
#include "observables.h"
#include "geometry.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "error.h"
#include "linear_algebra.h"
#include "spinor_field.h"
#include "communications.h"
#include "memory.h"
#include "random.h"

#include "update.h"


void spinor_scalarfield_mult_add_assign(spinor_field *out,scalar_field *sigma,double rho, spinor_field *in);


//Basic observables involving auxiliary fields
void ff_observables()
{
  double s_ave=0, ss_ave=0, p_ave=0, pp_ave=0;
  static int cnfg=1;
   
  _MASTER_FOR_SUM(ff_sigma->type,ix,s_ave,ss_ave,p_ave,pp_ave) {
      double s = *_FIELD_AT(ff_sigma,ix);
      double p = *_FIELD_AT(ff_pi,ix);
      s_ave += s;
      p_ave += p;
      ss_ave += s*s;
      pp_ave += p*p;
    }

    global_sum(&s_ave,1);
    global_sum(&p_ave,1);
    global_sum(&ss_ave,1);
    global_sum(&pp_ave,1);


    //Normalization as in the action, making <sigma> roughly equal to <psi^dagger psi>
    s_ave /= GLB_VOLUME;
    p_ave /= GLB_VOLUME;
    ss_ave /= GLB_VOLUME;
    pp_ave /= GLB_VOLUME;

  
  lprintf("AUXFIELDS",0,"sigma %1.8e sigma^2 %1.8e pi %1.8e pi^2 %1.8e\n",s_ave, ss_ave, p_ave, pp_ave);


#if 0 //rotate the pi field to zero, order a error
  double norm = 1./(s_ave*s_ave+p_ave*p_ave);
  if(G_4fermi > 0) { //Correlators 
    int ix,t,x,y,z,tc;
    for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
      ix=ipt(t,x,y,z);
      double s = *_FIELD_AT(ff_sigma,ix);
      double p = *_FIELD_AT(ff_pi,ix);
      
      double s_rot = (p*p_ave+s*s_ave)*norm;
      double p_rot = (p*s_ave-s*p_ave)*norm;

      *_FIELD_AT(ff_sigma,ix) = s_rot;
      *_FIELD_AT(ff_pi,ix) = p_rot;
    }
  }

if(G_4fermi > 0) { 
    _MASTER_FOR_SUM(ff_sigma->type,ix,s_ave,ss_ave,p_ave,pp_ave) {
      double s = *_FIELD_AT(ff_sigma,ix);
      double p = *_FIELD_AT(ff_pi,ix);
      s_ave += s;
      p_ave += p;
      ss_ave += s*s;
      pp_ave += p*p;
    }

    global_sum(&s_ave,1);
    global_sum(&p_ave,1);
    global_sum(&ss_ave,1);
    global_sum(&pp_ave,1);


    //Normalization as in the action, making <sigma> roughly equal to <psi^dagger psi>
    s_ave /= GLB_VOLUME*G_4fermi*G_4fermi;
    p_ave /= GLB_VOLUME*G_4fermi*G_4fermi;
    ss_ave /= GLB_VOLUME*G_4fermi*G_4fermi;
    pp_ave /= GLB_VOLUME*G_4fermi*G_4fermi;

  }
  lprintf("AUXFIELDS_ROT",0,"sigma %1.8e sigma^2 %1.8e pi %1.8e pi^2 %1.8e\n",s_ave, ss_ave, p_ave, pp_ave);

#endif

    int ix,t,x,y,z,tc;
    double sc[GLB_T],pc[GLB_T],spc[GLB_T];
    for (t=0; t<GLB_T; t++) sc[t]=pc[t]=spc[t]=0;
    for (t=0; t<T; t++) {
      tc = zerocoord[0]+t;
      for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
	ix=ipt(t,x,y,z);
	double s,p;
	s = *_FIELD_AT(ff_sigma,ix);
	p = *_FIELD_AT(ff_pi,ix);
	
	sc[tc] += s;
	pc[tc] += p;
	spc[tc]+= s*p;
	
      }
    } 
    
    // Just print out the sum over 3d volume
    lprintf("AUX_SIGMA_T",0," #%d =",cnfg);
    for (t=0; t<GLB_T; t++) {
     global_sum(&sc[t],1);
     lprintf("AUX_SIGMA_T",0," %1.8e ", sc[t] );
    }
    lprintf("AUX_PI_T",0," #%d =",cnfg);
    for (t=0; t<GLB_T; t++) {
     global_sum(&pc[t],1);
     lprintf("AUX_PI_T",0," %1.8e", pc[t]);
    }
    lprintf("AUX_SP_T",0," #%d =",cnfg);
    for (t=0; t<GLB_T; t++) {
     global_sum(&spc[t],1);
     lprintf("AUX_SP_T",0," %1.8e", spc[t]);
    }

    
  
#ifndef WITH_MPI 
  //Point to Point with no mpi
    for (t=0; t<GLB_T; t++) sc[t]=pc[t]=0;
    for (t=0; t<T; t++) {
      for (int t0=0; t0<GLB_T; t0++) {
        tc = (t+t0)%GLB_T;
        for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
          ix=ipt(t,x,y,z);
	  int ix0=ipt(tc,x,y,z);
	  double s,p,s0,p0;
	  s = *_FIELD_AT(ff_sigma,ix);
	  p = *_FIELD_AT(ff_pi,ix);
	  s0 = *_FIELD_AT(ff_sigma,ix0);
	  p0 = *_FIELD_AT(ff_pi,ix0);
	  
	  sc[t0] += s*s0;
	  pc[t0] += p*p0;
	
        }
      }
    }
    
    // Just print out the sum over 3d volume
    lprintf("AUX_SIGMA_PP",0," #%d =",cnfg);
    for (t=0; t<GLB_T; t++) {
     global_sum(&sc[t],1);
     lprintf("AUX_SIGMA_PP",0," %1.8e ", sc[t] );
    }
    lprintf("AUX_PI_PP",0," #%d =",cnfg);
    for (t=0; t<GLB_T; t++) {
     global_sum(&pc[t],1);
     lprintf("AUX_PI_PP",0," %1.8e", pc[t]);
    }

    lprintf("AUX_PI_PP",0,"\n");
    
  
#endif
  
  cnfg++;

}




