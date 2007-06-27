/*******************************************************************************
*
* Test of modules
*
*******************************************************************************/

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

int nhb,nor,nit,nth,nms,level,seed;
float beta;

void zero_corr(double *c){
  int i;
  for(i=0;i<T;++i)
    c[i]=0.;
}

void add_corr(double *c1, float *c2){
  int i;
  for(i=0;i<T;++i)
    c1[i]+=c2[i];
}

void print_usage(int argc,char *argv[]){
  printf("Usage: %s prop_file\n",argv[0]);
}

int main(int argc,char *argv[])
{
  int i,k,n;

  const int ncorrs = 1;
  
  suNf_spinor ***quark_prop;
  double **corr[ncorrs];
  FILE *propfile;

  float *m;
  int nm;

  if(argc!=2) {
    print_usage(argc,argv);
    return 0;
  }

  printf("Gauge group: SU(%d)\n",NG);
  printf("Fermion representation: dim = %d\n",NF);
  printf("The lattice size is %dx%d^3\n",T,L);

  printf("Computing Mesons corr functions\nfrom file: [%s]\n",argv[1]);

  error((propfile = fopen(argv[1], "rb"))==NULL,1,"Main",
  "Failed to open propagator file for reading");
  fread(&nm,(size_t) sizeof(int),1,propfile);
  m=(float*)malloc(sizeof(float)*nm);
  for(i=0;i<nm;++i)
  fread(m+i,(size_t) sizeof(float),1,propfile);

  printf("Found %d masses: ",nm);
  for(i=0;i<nm;++i)
  printf("%e ",m[i]);
  printf("\n\n");

  fflush(stdout);

  /*
  rlxs_init(level,seed);
  */

  geometry_eo_lexi();
  /*geometry_blocked();*/
  test_geometry();

  /* setup for quark propagators measures */
  quark_prop=(suNf_spinor***)malloc(sizeof(suNf_spinor**)*nm);
  for (k=0;k<nm;++k){
    quark_prop[k]=(suNf_spinor**)malloc(sizeof(suNf_spinor*)*4*NF);
    for (n=0;n<4*NF;++n)
      quark_prop[k][n]=(suNf_spinor*)malloc(sizeof(suNf_spinor)*VOLUME);
  }
  
  for (k=0;k<ncorrs;++k){
    corr[k]=(double**)malloc(sizeof(double*)*nm);
    for (n=0;n<nm;++n){
      corr[k][n]=(double*)malloc(sizeof(double)*T);
      zero_corr(corr[k][n]);
    }
  }


  i=0;
  /* Misure */
  do {
  
    float tmpcorr[T];
    
    for (n=0;n<4*NF;++n){
      for (k=0;k<nm;++k){
        error((fread(quark_prop[k][n],(size_t) sizeof(suNf_spinor),
            (size_t)(VOLUME),propfile)!=(VOLUME))&&(!feof(propfile)),
            1,"Main", "Failed to read form quark propagator!");
      }
    }
    
    for (k=0;k<nm;++k){
      if (!feof(propfile)){
        g5_correlator(tmpcorr, quark_prop[k]); /* misura pi */
        add_corr(corr[0][k],tmpcorr);
        
        printf("[%d] mass=%2.4f pi_corr= ",i,m[k]);
        for(n=0;n<T;++n) {
          printf("%e ",tmpcorr[n]);
        }
        printf("\n");
        fflush(stdout);
      }
    }
   
    ++i;
    
  } while(!feof(propfile));

  fclose(propfile);


  for (k=0;k<nm;++k){
    for (n=0;n<4*NF;++n)
      free(quark_prop[k][n]);
    free(quark_prop[k]);
  }
  free(quark_prop);
  
  for (k=0;k<ncorrs;++k){
    for (n=0;n<nm;++n)
      free(corr[k][n]);
    free(corr[k]);
  }

  free(m);

  return 0;
}
