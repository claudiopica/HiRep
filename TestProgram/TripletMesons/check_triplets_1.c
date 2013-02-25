
/******************************************************************************
 *
 * File check_triplets_1.c
 *
 * Check of the triplet mesons (free case)
 *
 * Author: Agostino Patella
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


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif

#error "Old version of Mesons, it should be updated"

enum MesonT        {A=0,Pi, Rho,B, Pi2,Rho2,Xt,Yt};
static int gid[8]= {1, -1, -1, -1,  1,  1,  1, -1}; /* gid = tr \bar{Gamma} Gamma / 4 */
static int g0[8] = {1,  1,  1, -1, -1, -1,  1, -1}; /* g0 = tr Gamma Gamma / 4 */
static int g1[8] = {1,  1, -1, -1,  1, -1, -1, +1}; /* g1 = tr gamma_1 \bar{Gamma} gamma_1 Gamma / 4 */
static int g2[8] = {1,  1,  1,  1,  1,  1, -1, -1}; /* g2 = tr gamma_2 \bar{Gamma} gamma_2 Gamma / 4 */
static int g3[8] = {1,  1,  1,  1,  1,  1, -1, -1}; /* g3 = tr gamma_3 \bar{Gamma} gamma_3 Gamma / 4 */
static int gb[8] = {1,  1,  1, -1, -1, -1,  1, -1};
static char nameT[8][256] = {"a", "pi", "rho", "b", "pi2", "rho2", "forbidden triplet 0+-", "forbidden triplet 1++" };

/* Mesons parameters */
typedef struct _input_mesons {
  char mstring[256];

  /* for the reading function */
  input_record_t read[2];
  
} input_mesons;

#define init_input_mesons(varname)					\
  {									\
    .read={								\
      {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring}, \
      {NULL, NULL, INT_T, NULL}						\
    }									\
  }


static double mass;
void free_correlators(double **triplets);


input_glb glb_ip = init_input_glb(glb_ip);
input_mesons mes_ip = init_input_mesons(mes_ip);


int main(int argc,char *argv[])
{
  int i, t;
  int g[4];
  double *ex_triplets[8];
  double *pta_triplets[8];
  char pame[256];
  spinor_field **pta_qprop=0;


  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  /* if (PID!=0) { logger_disable(); } */
  if (PID!=0) {
    logger_disable();}
  else{
    sprintf(pame,">out_%d",PID); logger_stdout(pame);
    sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
  }
  
  logger_map("DEBUG","debug");

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  read_input(glb_ip.read,"input_file");
  read_input(mes_ip.read,"input_file");

  strcpy(pame,mes_ip.mstring);
  mass=atof(strtok(pame, ";"));       


  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  lprintf("MAIN",0,"Zerocoord{%d,%d,%d,%d}\n",zerocoord[0],zerocoord[1],zerocoord[2],zerocoord[3]);

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

  error(!(GLB_X==GLB_Y && GLB_X==GLB_Z),1,"main", "This test works only for GLB_X=GLB_Y=GLB_Z");

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);

  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_ip.rlxd_level,glb_ip.rlxd_seed);
  rlxd_init(glb_ip.rlxd_level,glb_ip.rlxd_seed+PID);
 
  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif
  represent_gauge_field();

  lprintf("MAIN",0,"Mass = %f\n",mass);

  ex_triplets[0]=(double*)malloc(8*GLB_T*sizeof(double));
  pta_triplets[0]=(double*)malloc(8*GLB_T*sizeof(double));
  for(i=0; i<8; i++) {
    ex_triplets[i]=ex_triplets[0]+i*GLB_T;
    pta_triplets[i]=pta_triplets[0]+i*GLB_T;
  }
  pta_qprop=(spinor_field**)malloc(sizeof(spinor_field*));
  pta_qprop[0]=alloc_spinor_field_f(4*NF,&glattice);


  /* CALCOLO ESPLICITO */
  free_correlators(ex_triplets);


  /* MESONI CON PROPAGATORE POINT-TO-ALL */

  g[0]=g[1]=g[2]=g[3]=0;
  pta_qprop_QMR(g,pta_qprop, 1, &mass, 1e-9);
	
  id_correlator(pta_triplets[A], g[0], pta_qprop[0]);
  g0_correlator(pta_triplets[Xt], g[0], pta_qprop[0]);
  g5_correlator(pta_triplets[Pi], g[0], pta_qprop[0]);
  g0g5_correlator(pta_triplets[Pi2], g[0], pta_qprop[0]);
  g1_correlator(pta_triplets[Rho], g[0], pta_qprop[0]);
  g0g1_correlator(pta_triplets[Rho2], g[0], pta_qprop[0]);
  g5g1_correlator(pta_triplets[Yt], g[0], pta_qprop[0]);
  g0g5g1_correlator(pta_triplets[B], g[0], pta_qprop[0]);


  /* STAMPA */
   
  /* 	error = 0.; */

  lprintf("TEST",0,"\nANALITICO\tPOINT-TO-ALL\tERROR (must be less than 1e-9)\n");
  for(i=0; i<8; i++) {
    lprintf("TEST",0,"TRIPLET CORRELATOR %s\n", nameT[i]);
    for(t=0; t<GLB_T; t++)
      lprintf("TEST",0,"%e\t%e\t%e\n",
	      ex_triplets[i][t],
	      pta_triplets[i][t],
	      fabs(ex_triplets[i][t]-pta_triplets[i][t]));
  }


  /*   error = sqrt(error); */
  /*   printf("Error = %f\n", error); */

  finalize_process();

 
  exit(0);
}


double square(double x) {
  return x*x;
}


double re_ev(double k[4]) {
  return mass+4.0+cos((2.0*M_PI*k[0])/GLB_T)+cos((2.0*M_PI*k[1])/GLB_X)+cos((2.0*M_PI*k[2])/GLB_Y)+cos((2.0*M_PI*k[3])/GLB_Z);
}


double im_ev(double k[4]) {
  return sqrt(
	      square(sin((2.0*M_PI*k[0])/GLB_T))+
	      square(sin((2.0*M_PI*k[1])/GLB_X))+
	      square(sin((2.0*M_PI*k[2])/GLB_Y))+
	      square(sin((2.0*M_PI*k[3])/GLB_Z)));
}



void free_correlators(double **triplets) {
  double A2[GLB_T], B2[4][GLB_T];
  complex A[GLB_T], B[4][GLB_T];
  double norm2, re, im, ct, st, z;
  int i, t;
  double k[4];
  double sigma[4] = {0.,0.,0.,0.};
#ifdef ANTIPERIODIC_BC_T
  sigma[0] = .5;
#endif
#ifdef ANTIPERIODIC_BC_X
  sigma[1] = .5;
#endif
#ifdef ANTIPERIODIC_BC_Y
  sigma[2] = .5;
#endif
#ifdef ANTIPERIODIC_BC_Z
  sigma[3] = .5;
#endif
  
  z = - 4.*NF / square(GLB_T*GLB_X*GLB_Y*GLB_Z);
   
  lprintf("FREE",0,"sigma = (%f,%f,%f,%f)\n",sigma[0],sigma[1],sigma[2],sigma[3]);
  for(t = 0; t < GLB_T; t++) {
    A2[t] = 0.;
    B2[0][t] = B2[1][t] = B2[2][t] = B2[3][t] = 0.;
  }
   
  for(k[1] = sigma[1]; k[1] < GLB_X+sigma[1]-.5; k[1] += 1.)
    for(k[2] = sigma[2]; k[2] < GLB_Y+sigma[2]-.5; k[2] += 1.)
      for(k[3] = sigma[3]; k[3] < GLB_Z+sigma[3]-.5; k[3] += 1.) {
      
      
	for(t = 0; t < GLB_T; t++) {
	  A[t].re = A[t].im = 0.;
	  B[0][t].re = B[0][t].im = 0.;
	  B[1][t].re = B[1][t].im = 0.;
	  B[2][t].re = B[2][t].im = 0.;
	  B[3][t].re = B[3][t].im = 0.;
	}
	for(k[0] = sigma[0]; k[0] < GLB_T+sigma[0]-.5; k[0] += 1.) {
	  re = re_ev(k);
	  im = im_ev(k);
	  norm2 = re*re+im*im;

	  for(t = 0; t < GLB_T; t++) {
            ct = cos((2.0*M_PI*t*k[0])/GLB_T);
            st = sin((2.0*M_PI*t*k[0])/GLB_T);
            
            A[t].re += re * ct / norm2;
            A[t].im += re * st / norm2;
            B[0][t].re += sin((2.0*M_PI*k[0])/GLB_T) * ct / norm2;
            B[0][t].im += sin((2.0*M_PI*k[0])/GLB_T) * st / norm2;
            B[1][t].re += sin((2.0*M_PI*k[1])/GLB_X) * ct / norm2;
            B[1][t].im += sin((2.0*M_PI*k[1])/GLB_X) * st / norm2;
            B[2][t].re += sin((2.0*M_PI*k[2])/GLB_Y) * ct / norm2;
            B[2][t].im += sin((2.0*M_PI*k[2])/GLB_Y) * st / norm2;
            B[3][t].re += sin((2.0*M_PI*k[3])/GLB_Z) * ct / norm2;
            B[3][t].im += sin((2.0*M_PI*k[3])/GLB_Z) * st / norm2;
	  }
	}
	for(t = 0; t < GLB_T; t++) {
#define cmplx_norm2(z) ( (z).re * (z).re + (z).im * (z).im )
	  A2[t] += cmplx_norm2(A[t]) * z;
	  B2[0][t] += cmplx_norm2(B[0][t]) * z;
	  B2[1][t] += cmplx_norm2(B[1][t]) * z;
	  B2[2][t] += cmplx_norm2(B[2][t]) * z;
	  B2[3][t] += cmplx_norm2(B[3][t]) * z;
#undef cmplx_norm2
	}
      
      }

  lprintf("FREE",10,"A2\tB2[0]\tB2[1]\tB2[2]\tB[3]\n",A2[0]);
  for(t = 0; t < GLB_T; t++)
    lprintf("FREE",10,"%e\t%e\t%e\t%e\t%e\n",A2[t],B2[0][t],B2[1][t],B2[2][t],B2[3][t]);
   
  for(i = 0; i < 8; i++) {
    for(t = 0; t < GLB_T; t++) {
      triplets[i][t] = gb[i]*(gid[i]*A2[t] - g0[i]*B2[0][t] - g1[i]*B2[1][t] - g2[i]*B2[2][t] - g3[i]*B2[3][t]);
    }
  }

  lprintf("FREE",0,"Exact free correlators computed.\n");
   
}

