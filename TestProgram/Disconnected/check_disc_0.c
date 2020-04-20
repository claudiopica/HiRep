
/******************************************************************************
*
* File check_disc_0.c
*
* Check of the  disc loops (free case): discon volume source (type = 0)
*
* Author: Vincent Drach
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
#include "setup.h"
#include "spectrum.h"
#include "clover_tools.h"
#include "disconnected.h"
#include "communications.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif



static double complex gid[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. }; /* gid = tr Gamma  */
static double complex g0[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. }; /*  g0 = tr gamma_0 Gamma */
static double complex g1[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. }; /*  g1 = tr gamma_1 Gamma */
static double complex g2[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. }; /*  g2 = tr gamma_2 Gamma */
static double complex g3[16] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. }; /*  g3 = tr gamma_3 Gamma */

char* mes_channel_names[16]={"g5","g1","g2","g3","-ig0g5","-ig0g1","-ig0g2","-ig0g3","id","-ig5g1","-ig5g2","-ig5g3","g0","-ig5g0g1","-ig5g0g2","-ig5g0g3"};

#define mult_mat(r,A,B) \
{ \
  int _i, _j, _k; \
  double complex wm[4][4]; \
  for(_i=0; _i<4; _i++) \
  for(_j=0; _j<4; _j++) { \
    _complex_0(wm[_i][_j]); \
    for(_k=0; _k<4; _k++) { \
      wm[_i][_j] += A[_i][_k]*B[_k][_j];\
    } \
  } \
  for(_i=0; _i<4; _i++) \
  for(_j=0; _j<4; _j++) { \
    r[_i][_j] = wm[_i][_j];\
  } \
}

#define mult_mat_minusI(r,A) \
{ \
  int _i, _j; \
  for(_i=0; _i<4; _i++) \
  for(_j=0; _j<4; _j++) { \
    r[_i][_j] = -I*A[_i][_j];\
  } \
}



#define set_zero_mat(A) \
{ \
  int _i,_j; \
  for(_i=0; _i<4; _i++) \
  for(_j=0; _j<4; _j++) { \
    A[_i][_j]=0.; \
  } \
}
#define trace_mat(r,A) \
{ \
  int _i; \
  _complex_0(r); \
  for(_i=0; _i<4; _i++) { \
    r+=A[_i][_i]; \
  } \
}

/* Mesons parameters */
typedef struct _input_mesons {
  char mstring[256];

  /* for the reading function */
  input_record_t read[7];
  double precision;
  int nhits;
  int source_type;
  int n_mom;
  double csw;


} input_mesons;

#define init_input_mesons(varname)					\
{									\
  .read={								\
    {"Fermion mass", "disc:mass = %s", STRING_T, (varname).mstring}, \
    {"inverter precision", "disc:precision = %lf", DOUBLE_T, &(varname).precision}, \
    {"number of inversions per cnfg", "disc:nhits0 = %d", INT_T, &(varname).nhits}, \
    {"maximum component of momentum", "disc:n_mom = %d", INT_T, &(varname).n_mom}, \
    {"csw coefficient", "disc:csw = %lg",DOUBLE_T, &(varname).csw},	\
    {NULL, NULL,INT_T, NULL}				\
  }							\
}


static double mass;
void free_loops(double complex **loops,int n_mom);

input_mesons mes_ip = init_input_mesons(mes_ip);

char char_t[100];
FILE *fp;
char path[1035];


/*Gamma / 4 */
double complex get_gid(double complex Gamma[4][4]){
  double complex r;
  trace_mat(r,Gamma );
  return r;
}
/* gmu = tr Gamma gamma_mu   */
double complex get_gmu(double complex Gamma[4][4],int mu){
  int  sign;
  double complex tmp[4][4];
  double complex r;
  double complex gmu[4][4];


  if (mu==1)  g1_debug(gmu, &sign);
  if (mu==2)  g2_debug(gmu, &sign);
  if (mu==3)  g3_debug(gmu, &sign);
  if (mu==0)  g0_debug(gmu, &sign);

  mult_mat(tmp,Gamma,gmu);

  trace_mat(r,tmp);
  return r;
}


int compare_corr_abs(double complex  * corr_ex, double complex * corr_num, int ip, char* name[16], double tol ){
    int retval = 0;
    int nGamma = 16;
    for(int n=0; n<nGamma; n++){  
           if(cabs(corr_ex[n] - corr_num[n]) > tol) 
            {
                lprintf("TEST",0,"Mismatch %s, ip=%d, absolute diff: %e, numeric = %e + I*(%e), analytic = %e + I*(%e) \n",name[n],ip,cabs(corr_ex[n] - corr_num[n]),creal(corr_num[n]),cimag(corr_num[n]),creal(corr_ex[n]),cimag(corr_ex[n]));
                retval += 1;
            }
            else
            {
                 lprintf("TEST",0,"Match %s, ip=%d, numeric = %e + I*(%e), analytic = %e + I*(%e) \n",name[n],ip,creal(corr_num[n]),cimag(corr_num[n]),creal(corr_ex[n]),cimag(corr_ex[n]));
            }
    }  
    return retval;
}

int main(int argc,char *argv[])
{
  int i,sign;
  char pame[256];
  int source_type=0;
  int return_value=0;
  int n_Gamma=16;
  int n_mom_tot = 0;
  int nhits;
  double complex*** out_corr;
  double complex **ex_loops;
  double complex **mean_loops;
  double abs_tol=1e-1;
  double rel_tol_scalar_loop=1e-3;
  struct timeval start, end, etime;


  double complex g[16][4][4];
  double complex tmp[4][4];
  g5_debug(g[0], &sign);
  g1_debug(g[1], &sign);
  g2_debug(g[2], &sign);
  g3_debug(g[3], &sign);
  g0g5_debug(tmp, &sign);
  mult_mat_minusI(g[4],tmp);
  g0g1_debug(tmp, &sign);
  mult_mat_minusI(g[5],tmp);
  g0g2_debug(tmp, &sign);
  mult_mat_minusI(g[6],tmp);
  g0g3_debug(tmp, &sign);
  mult_mat_minusI(g[7],tmp);
  id_debug(g[8],&sign);
  g5g1_debug(tmp, &sign);
  mult_mat_minusI(g[9],tmp);
  g5g2_debug(tmp, &sign);
  mult_mat_minusI(g[10],tmp);
  g5g3_debug(tmp, &sign);
  mult_mat_minusI(g[11],tmp);
  g0_debug(g[12], &sign);
  mult_mat(g[13],g[0],g[5]);//  -i g5g0g1
  mult_mat(g[14],g[0],g[6]);//  -i g5g0g2
  mult_mat(g[15],g[0],g[7]);//  -i g5g0g3

  for (i=0;i<16;i++){
    gid[i] =  get_gid(g[i]);
    g0[i] =  get_gmu(g[i],0);
    g1[i] =  get_gmu(g[i],1);
    g2[i] =  get_gmu(g[i],2);
    g3[i] =  get_gmu(g[i],3);
  }

  /* setup process id and communications */
  setup_process(&argc,&argv);

  setup_gauge_fields();

  read_input(mes_ip.read,get_input_filename());
 
  strcpy(pame,mes_ip.mstring);
  mass=atof(strtok(pame, ";"));
  nhits  = mes_ip.nhits;
  n_mom_tot = mes_ip.n_mom*mes_ip.n_mom*mes_ip.n_mom;
  lprintf("MAIN",0,"disc:mass = %f\n",mass);
  lprintf("MAIN",0,"nhits = %i\n",nhits);
  lprintf("MAIN",0,"Inverter precision = %e\n",mes_ip.precision);
  lprintf("MAIN",0,"Number of momenta = %d\n",mes_ip.n_mom);
  #if defined(WITH_CLOVER) ||  defined(WITH_EXPCLOVER)
  set_csw(&mes_ip.csw);
  #endif
  
  gettimeofday(&start,0);

  unit_u(u_gauge);
  represent_gauge_field();
  #ifdef REPR_FUNDAMENTAL 
  apply_BCs_on_represented_gauge_field(); //This is a trick: the BCs are not applied in the case the REPR is fundamental because represent_gauge field assumes that the right BCs are already applied on the fundamental field!
  #endif
  
  

  lprintf("MAIN",0,"source type is fixed to 0: volume source   \n");
  lprintf("MAIN",0,"Measuring D(t) =  sum_x psibar(x) Gamma psi(x)\n");
  lprintf("MAIN",0,"Zerocoord{%d,%d,%d,%d}\n",zerocoord[0],zerocoord[1],zerocoord[2],zerocoord[3]);

  error(!(GLB_X==GLB_Y && GLB_X==GLB_Z),1,"main", "This test works only for GLB_X=GLB_Y=GLB_Z");

  lprintf("CORR",0,"Number of noise vector : nhits = %i \n", nhits);
  
  
  out_corr=(double complex***) malloc(sizeof(double complex**)*n_mom_tot);
  
  for(int i=0; i<n_mom_tot; i++)    out_corr[i]=(double complex**) malloc(sizeof(double complex*)*n_Gamma);
  for(int i=0; i<n_mom_tot; i++) for(int j=0; j<n_Gamma; j++) out_corr[i][j]=(double complex*) calloc(GLB_T,sizeof(double complex));
  

  ex_loops = (double complex **)malloc(sizeof(double complex*)*n_mom_tot) ;
  for(int i=0; i<n_mom_tot; i++) ex_loops[i] = (double complex *)calloc(n_Gamma,sizeof(double complex));
  
  mean_loops = (double complex **)malloc(sizeof(double complex*)*n_mom_tot) ;
  for(int i=0; i<n_mom_tot; i++) mean_loops[i] = (double complex *)calloc(n_Gamma,sizeof(double complex));
 
  measure_loops(1, &mass, nhits,0,  mes_ip.precision,source_type,mes_ip.n_mom,out_corr);

  //stochastic & time average 
  for(int i=0; i<n_mom_tot; i++) for(int j=0; j<n_Gamma; j++)   for(int t=0; t<GLB_T; t++) mean_loops[i][j] += out_corr[i][j][t]/(nhits*GLB_T);
  
  //  /* CALCOLO ESPLICITO */
  free_loops(ex_loops,mes_ip.n_mom);
  
  for(i=0; i<n_mom_tot; i++) return_value += compare_corr_abs(ex_loops[i],mean_loops[i],i,mes_channel_names,abs_tol);
  global_sum_int(&return_value,1);
 

  // The zero momentum scalar loops should agree to a much higher precision
  if (   cabs(ex_loops[0][8] - mean_loops[0][8])/cabs(ex_loops[0][8]) > rel_tol_scalar_loop  )  return_value +=1;

  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("MAIN",0,"Configuration : analysed in [%ld sec %ld usec]\n",etime.tv_sec,etime.tv_usec);
  
  global_sum_int(&return_value,1);

  lprintf("MAIN", 0, "return_value= %d\n ",  return_value);

  finalize_process();

  return return_value;
}


double square(double x) {
  return x*x;
}


double denom(double k[4]) {
  double res;

  res = square( mass+4.0 - cos((2.0*M_PI*k[0])/GLB_T)-cos((2.0*M_PI*k[1])/GLB_X)-cos((2.0*M_PI*k[2])/GLB_Y)-cos((2.0*M_PI*k[3])/GLB_Z) )
  +  square(sin((2.0*M_PI*k[0])/GLB_T)) + square(sin((2.0*M_PI*k[1])/GLB_X)) + square(sin((2.0*M_PI*k[2])/GLB_Y)) + square(sin((2.0*M_PI*k[3])/GLB_Z)) ;
  return res;
}
double complex compute_phase(int px,int py, int pz){

  double complex res,phase;
  double x,y,z,pdotx;

  _complex_0(res);
  for(x = 0.; x < GLB_X; x += 1.) {
    for(y = 0.; y < GLB_Y; y += 1.) {
      for(z = 0.; z < GLB_Z; z += 1.) {
        pdotx = 2.0*M_PI*( ((double) px)*x/GLB_X + ((double) py)*y/GLB_Y + ((double) pz)*z/GLB_Z);
        phase = cos(pdotx) +I*sin(pdotx);
        res += phase;

      }
    }
  }
  //printf("%e %e %e \n",creal(res*conj(res)),cimag(res*conj(res)), ((double) px)*x/GLB_X);
  return res;
}
void free_loops(double complex **loops,int n_mom) {
  double  A, B[4];
  double  tmp,norm;
  double complex phase;
  int i,ll;
  int px,py,pz;
  double k[4];
  double sigma[4] = {0.,0.,0.,0.};

  #ifdef BC_T_ANTIPERIODIC
  sigma[0] = .5;
  #endif
  #ifdef BC_X_ANTIPERIODIC
  sigma[1] = .5;
  #endif
  #ifdef BC_Y_ANTIPERIODIC
  sigma[2] = .5;
  #endif
  #ifdef BC_Z_ANTIPERIODIC
  sigma[3] = .5;
  #endif


  lprintf("FREE",0,"sigma = (%f,%f,%f,%f)\n",sigma[0],sigma[1],sigma[2],sigma[3]);
  norm = (double) NF/GLB_T;
  A=B[0]=B[1]=B[2]=B[3]=0.;

  for(k[1] = sigma[1]; k[1] < GLB_X+sigma[1]-.5; k[1] += 1.)
  for(k[2] = sigma[2]; k[2] < GLB_Y+sigma[2]-.5; k[2] += 1.)
  for(k[3] = sigma[3]; k[3] < GLB_Z+sigma[3]-.5; k[3] += 1.) {

    for(k[0] = sigma[0]; k[0] < GLB_T+sigma[0]-.5; k[0] += 1.) {

      tmp= denom(k);

      A += (mass+4.0 - cos((2.0*M_PI*k[0])/GLB_T)-cos((2.0*M_PI*k[1])/GLB_X)-cos((2.0*M_PI*k[2])/GLB_Y)-cos((2.0*M_PI*k[3])/GLB_Z))/tmp;
      B[0] += sin((2.0*M_PI*k[0])/GLB_T)/tmp;
      B[1] += sin((2.0*M_PI*k[1])/GLB_X)/tmp;
      B[2] += sin((2.0*M_PI*k[2])/GLB_Y)/tmp;
      B[3] += sin((2.0*M_PI*k[3])/GLB_Z)/tmp;

    }
  }


  lprintf("FREE",10,"A=%e\t B[0=]%e\t B[1]=%e\t B[2]=%e\t B[3]=%e\n",A,B[0],B[1],B[2],B[3]);
  ll=0;
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz) {

    phase = compute_phase(px,py,pz); // this is a delta (p,x)
    phase = phase/(GLB_X*GLB_Y*GLB_Z);

    for(i = 0; i < 16; i++)   loops[ll][i] = phase*(gid[i]*A-I*(g0[i]*B[0] + g1[i]*B[1] + g2[i]*B[2] + g3[i]*B[3]))*norm ; // = sum_vec{x} psi(x) Gamma psibar(x)]
    ll++;
  }

  lprintf("FREE",0,"Exact free correlators computed.\n");
}
