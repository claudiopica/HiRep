
/******************************************************************************
*
* File check13.c
*
* Check of the mesons with the Dublin algorithm (free case)
*
* Author: Agostino Patella
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "linear_algebra.h"
#include "update.h"
#include "inverters.h"
#include "dirac.h"
#include "representation.h"
#include "global.h"
#include "logger.h"
#include "observables.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif


static double mass = 0.01;
static const double EPSILON = 1.e-12;

enum MesonT {A, Pi, Rho, B, Pi2, Rho2, X, Y};
static int g0[8] = {1, -1, -1, 1, -1, -1, 1, 1};
static int gs[8] = {3, -3, -1, -1, 3, 1, -3, 1};
static char nameT[8][256] = {"a", "pi", "rho", "b", "pi2", "rho2", "forbidden triplet 0+-", "forbidden triplet 1++" };
enum MesonS {F, Eta, Phi, H, Eta2, Phi2, Xs, Ys};
static int p[8] = {16, 0, 0, 0, 0, 0, 0, 0 };
static char nameS[8][256] = {"f", "eta", "phi", "h", "eta2", "phi2", "forbidden singlet 0+-", "forbidden singlet 1++" };


void free_correlators(double triplets[8][T], double singlets[8][T]);


int main(int argc,char *argv[])
{
	int i, n, t;
   FILE *log=NULL;   
	double m[1];
	double error;
	int nm=1;
	
	enum{  n_correlators = 17 };
	double **correlator[n_correlators];
	char corr_name[n_correlators][256];
	
	
	for(i=0; i<n_correlators; i++) {
		correlator[i] = (double**)malloc(sizeof(double*)*nm);
		correlator[i][0] = (double*)malloc(sizeof(double)*nm*T);
		for(n=1;n<nm;n++)
			correlator[i][n] = correlator[i][n-1] + T;
	}
	
	m[0] = mass;

   log=freopen("check13.log","w",stdout);
   printf("\n");
   printf("Dubliner mesons (free case)\n");
   printf("---------------------------------------\n\n");
   printf("The lattice size is %dx%d^3\n",T,L);
   printf("size of the gluon rep: %d, size of the fermion rep: %d\n",NG,NF);
   printf("mass of the fermion: %f\n",mass);
   
   rlxd_init(1,12345);

	logger_setlevel(0,10);

   geometry_eo_lexi();
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();

   set_spinor_len(VOLUME);

	dublin_meson_correlators(correlator, corr_name, n_correlators, nm, m);
   
   for(i = 0; i < n_correlators; i++) {
      printf("CORRELATOR %s\n", corr_name[i]);
      for(t = 0; t < T; t++)
         printf("%f\n", correlator[i][0][t]);
   }


   printf("\n");
   fclose(log);
   
   exit(0);
}


double square(double x) {
   return x*x;
}


double re_ev(int k[4]) {
   return mass+4.0+cos((2.0*M_PI*k[0])/T)+cos((2.0*M_PI*k[1])/L)+cos((2.0*M_PI*k[2])/L)+cos((2.0*M_PI*k[3])/L);
}


int smod(int a, int b) {
   b = (b>=0) ? b : -b;
   return (a>=0) ? a%b : b-1-((-a-1)%b);
}


double im_ev(int k[4]) {
   if( T%2 == 0 ) {
      if( smod(k[0],T/2)==0 && smod(k[1],L/2)==0 && smod(k[2],L/2)==0 && smod(k[3],L/2)==0 )
         return 0.0;
   } else {
      if( smod(k[0],T)==0 && smod(k[1],L)==0 && smod(k[2],L)==0 && smod(k[3],L)==0 )
         return 0.0;
   }   
   return sqrt(
      square(sin((2.0*M_PI*k[0])/T))+
      square(sin((2.0*M_PI*k[1])/L))+
      square(sin((2.0*M_PI*k[2])/L))+
      square(sin((2.0*M_PI*k[3])/L)));
}


void free_correlators(double triplets[8][T], double singlets[8][T]) {
   double A2[T], B2[2][T], A2bar;
   double A[T], B[2][T], Abar;
   double norm2, re, im, ct, st;
   int i, t, k[4];
   
   Abar = 0.;
   for(t = 0; t < T; t++) {
      A2[t] = 0.;
      B2[0][t] = B2[1][t] = 0.;
   }
   
   for(k[1] = 0; k[1] < L; k[1]++)
   for(k[2] = 0; k[2] < L; k[2]++)
   for(k[3] = 0; k[3] < L; k[3]++) {
      
      
      for(t = 0; t < T; t++) {
         A[t] = 0.;
         B[0][t] = B[1][t] = 0.;
      }
      for(k[0] = 0; k[0] < T; k[0]++) {
         re = re_ev(k);
         im = im_ev(k);
         norm2 = re*re+im*im;
         
         Abar += re / norm2;

         for(t = 0; t < T; t++) {
            ct = cos((2.0*M_PI*t*k[0])/T);
            st = sin((2.0*M_PI*t*k[0])/T);
            
            A[t] += re * ct / norm2;
            B[0][t] += sin((2.0*M_PI*k[0])/T) * st / norm2;
            B[1][t] += sin((2.0*M_PI*k[1])/L) * st / norm2;
         }
      }
      for(t = 0; t < T; t++) {
         A2[t] += A[t] * A[t];
         B2[0][t] += B[0][t] * B[0][t];
         B2[1][t] += B[1][t] * B[1][t];
      }
      
   }
   
   A2bar = Abar * Abar;
   
   for(t = 0; t < T; t++) for(i = 0; i < 8; i++) {
      triplets[i][t] = - 4*g0[i] * (A2[t] + g0[i]*B2[0][t] + gs[i]*B2[1][t]);
      singlets[i][t] = triplets[i][t] + p[i] * A2bar;
   }

}



