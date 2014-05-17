/******************************************************************************
*
* Testing the spin structure & printing the gamma matrices
*
******************************************************************************/

#define MAIN_PROGRAM

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



#ifdef WITH_MPI
#error: check_diracoperator_4 only works only on serial jobs
#endif


#if defined(BC_T_ANTIPERIODIC) || defined(BC_X_ANTIPERIODIC) || defined(BC_Y_ANTIPERIODIC) || defined(BC_Z_ANTIPERIODIC)
#error This test does not work with antiperiodic boundary conditions
#endif

#if defined(BASIC_SF) || defined(ROTATED_SF)
#error This test does not work with SF
#endif

void D(spinor_field *out, spinor_field *in){
   Dphi_(out,in);
}

static int init_gamma = 0;
complex my_gamma[4][4][4];

void compute_gamma(int g[4], int ic) {
  int p[4], c[4], shift[4];
  double dbl;
  complex locmy_gamma[4][4][4];
  
  spinor_field *in, *out;
  in=alloc_spinor_field_f(1,&glattice);
  out=alloc_spinor_field_f(1,&glattice);

  for(int mu=0; mu<4; mu++)
  for(int beta=0; beta<4; beta++)
  for(int alpha=0; alpha<4; alpha++) {
    locmy_gamma[mu][alpha][beta].re = 0.;
    locmy_gamma[mu][alpha][beta].im = 0.;
  }
  
  p[0] = g[0]/T; p[1] = g[1]/X; p[2] = g[2]/Y; p[3] = g[3]/Z;
  
  for(int mu=0; mu<4; mu++) {
    shift[0]=shift[1]=shift[2]=shift[3]=0;
    shift[mu]=1;
    
    _MASTER_FOR(&glattice,ix) {
      for(int nu=0; nu<4; nu++) {
        _suNg_zero(*_4FIELD_AT(u_gauge,ix,nu));
      }
    }
    if(COORD[0]==p[0] && COORD[1]==p[1] && COORD[2]==p[2] && COORD[3]==p[3]) {
      c[0] = g[0]%T; c[1] = g[1]%X; c[2] = g[2]%Y; c[3] = g[3]%Z;
      int ix=ipt(c[0],c[1],c[2],c[3]);
      _suNg_unit(*_4FIELD_AT(u_gauge,ix,mu));
    }
    start_gf_sendrecv(u_gauge);
    represent_gauge_field();

    dbl=0.;
    _MASTER_FOR_SUM(&glattice,ix,dbl) {
      for(int nu=0; nu<4; nu++) {
        double tmp;
        _suNf_sqnorm(tmp,*_4FIELD_AT(u_gauge_f,ix,nu));
        dbl+=tmp;
      }
    }
    global_sum(&dbl,1);
    if(dbl!=NF)
      lprintf("ERROR",0,"u_gauge_f sqnorm=%f\n",dbl);

    for(int beta=0; beta<4; beta++) {
      spinor_field_zero_f(in);
      _MASTER_FOR(&glattice,ix) {
        _FIELD_AT(in,ix)->c[beta].c[ic].re=1.;
      }

      dbl=spinor_field_sqnorm_f(in);
      if(dbl!=GLB_T*GLB_X*GLB_Y*GLB_Z)
        lprintf("ERROR",0,"source sqnorm=%f\n",dbl);
      
      D(out,in);

      dbl=spinor_field_sqnorm_f(out);
      if(dbl==0.)
        lprintf("ERROR",0,"vanishing out sqnorm\n");

      for(c[0]=0; c[0]<T; c[0]++)
      for(c[1]=0; c[1]<X; c[1]++)
      for(c[2]=0; c[2]<Y; c[2]++)
      for(c[3]=0; c[3]<Z; c[3]++) {
        int ix=ipt(c[0],c[1],c[2],c[3]);
        if(c[0]+zerocoord[0]==g[0] &&
           c[1]+zerocoord[1]==g[1] &&
           c[2]+zerocoord[2]==g[2] &&
           c[3]+zerocoord[3]==g[3])
        {
          for(int alpha=0; alpha<4; alpha++) {
            locmy_gamma[mu][alpha][beta].re += _FIELD_AT(out,ix)->c[alpha].c[ic].re;
            locmy_gamma[mu][alpha][beta].im += _FIELD_AT(out,ix)->c[alpha].c[ic].im;
          }
        } else if(c[0]+zerocoord[0]==(g[0]+shift[0])%GLB_T &&
                  c[1]+zerocoord[1]==(g[1]+shift[1])%GLB_X &&
                  c[2]+zerocoord[2]==(g[2]+shift[2])%GLB_Y &&
                  c[3]+zerocoord[3]==(g[3]+shift[3])%GLB_Z)
        {
          for(int alpha=0; alpha<4; alpha++) {
            locmy_gamma[mu][alpha][beta].re -= _FIELD_AT(out,ix)->c[alpha].c[ic].re;
            locmy_gamma[mu][alpha][beta].im -= _FIELD_AT(out,ix)->c[alpha].c[ic].im;
          }
        } else {
          _spinor_prod_re_f(dbl,*_FIELD_AT(out,ix),*_FIELD_AT(out,ix));
          if(dbl!=0.)
            lprintf("ERROR",0,"Non-zero spinor! g=(%d,%d,%d,%d) ic=%d mu=%d beta=%d ix=%d sqnorm=%e\n",g[0],g[1],g[2],g[3],ic,mu,beta,ix,dbl);
        }
      }
    }
  }

  global_sum((double*)locmy_gamma,2*4*4*4);
/*
  counter=0;
  for(mu=0; mu<4; mu++)
  for(beta=0; beta<4; beta++)
  for(alpha=0; alpha<4; alpha++) {
    global_sum(&locmy_gamma[mu][alpha][beta].re,counter);
    global_sum(&locmy_gamma[mu][alpha][beta].im,counter);
    counter+=2;
  }
*/
  if(init_gamma==0) {
    for(int mu=0; mu<4; mu++) {
      lprintf("GAMMA",0,"gamma[%d] = \n",mu);
      for(int alpha=0; alpha<4; alpha++) {
        lprintf("GAMMA",0,"[ ");
        for(int beta=0; beta<4; beta++) {
          my_gamma[mu][alpha][beta]=locmy_gamma[mu][alpha][beta];
          lprintf("GAMMA",0,"(%.2f,%.2f) ",my_gamma[mu][alpha][beta].re,my_gamma[mu][alpha][beta].im);
        }
        lprintf("GAMMA",0,"]\n");
      }
    }
    init_gamma=1;
  } else {
    for(int mu=0; mu<4; mu++) {
      dbl=0.;
      for(int beta=0; beta<4; beta++)
      for(int alpha=0; alpha<4; alpha++) {
        dbl+=(locmy_gamma[mu][alpha][beta].re-my_gamma[mu][alpha][beta].re)*(locmy_gamma[mu][alpha][beta].re-my_gamma[mu][alpha][beta].re);
        dbl+=(locmy_gamma[mu][alpha][beta].im-my_gamma[mu][alpha][beta].im)*(locmy_gamma[mu][alpha][beta].im-my_gamma[mu][alpha][beta].im);
      }
      if(dbl!=0.) {
        lprintf("ERROR",0,"Wrong gamma matrix! g=(%d,%d,%d,%d) ic=%d mu=%d err2=%e\n",g[0],g[1],g[2],g[3],ic,mu,dbl);
        lprintf("ERROR",0,"gamma[%d] = \n",mu);
        for(int alpha=0; alpha<4; alpha++) {
          lprintf("ERROR",0,"[ ");
          for(int beta=0; beta<4; beta++) {
            lprintf("ERROR",0,"(%.2f,%.2f) ",locmy_gamma[mu][alpha][beta].re,locmy_gamma[mu][alpha][beta].im);
          }
          lprintf("ERROR",0,"]\n");
        }
      }
    }
  }
  
  free_spinor_field_f(in);
  free_spinor_field_f(out);
}

/* FUNZIONA SOLO CON PBC */

int main(int argc,char *argv[])
{
  char tmp[256];
  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */

  logger_setlevel(0,100); /* log all */
  if (PID!=0) { 
    logger_disable();}
  else{
    sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }
   
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  read_input(glb_var.read,"test_input");

  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

  /* setup lattice geometry */
  geometry_mpi_eo();
/*  test_geometry_mpi_eo(); */

  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  int g[4], ic;

  for(ic=0; ic<NF; ic++)
    for(g[0]=0; g[0]<GLB_T; g[0]++)
      for(g[1]=0; g[1]<GLB_X; g[1]++)
	for(g[2]=0; g[2]<GLB_Y; g[2]++)
	  for(g[3]=0; g[3]<GLB_Z; g[3]++) {
	    lprintf("MAIN",0,".");
	    
	    compute_gamma(g,ic);
	  }


  finalize_process();

}
