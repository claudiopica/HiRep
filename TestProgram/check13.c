
/******************************************************************************
*
* File check13.c
*
* Check of the mesons (both mesons.c and dublin_mesons.c) (free case)
*
* Author: Agostino Patella
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
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


static double mass = .5;
static const double EPSILON = 1.e-12;

enum MesonT {A=0, Pi, Rho, B, Pi2, Rho2, X, Y};
static int g0[8] = {1, -1, -1, 1, -1, -1, 1, 1};
static int gs[8] = {3, -3, -1, -1, 3, 1, -3, 1};
static int gb[8] = {1, -1, -1, -1, 1, 1, 1, -1};
static char nameT[8][256] = {"a", "pi", "rho", "b", "pi2", "rho2", "forbidden triplet 0+-", "forbidden triplet 1++" };
enum MesonS {F=8, Eta, Phi, H, Eta2, Phi2, Xs, Ys};
static int p[8] = {1, 0, 0, 0, 0, 0, 0, 0 };
static char nameS[8][256] = {"f", "eta", "phi", "h", "eta2", "phi2", "forbidden singlet 0+-", "forbidden singlet 1++" };


void free_correlators(double triplets[8][T], double singlets[8][T]);
void stat(int ndata, double *data, double *ave, double *error);


int main(int argc,char *argv[])
{
	int i, n, t, x, counter;
  FILE *log=NULL;   
	double m[1];
	int nm=1;
  double triplets[8][T], singlets[8][T];

	enum{  n_correlators = 17 };
	double ***dcorr_histo[n_correlators];
	double **dcorr_ave[n_correlators];
	double **dcorr_err[n_correlators];
	double **dcorr_tmp[n_correlators];
	char dcorr_name[n_correlators][256];

/* 	FILE *propfile; */
	suNf_spinor **quark_prop;
	double ptacorr[8][T];
	long int propsize;
	double norm;
	suNf_spinor *test;
	
	int nmeasures = 100;
	
	for(i=0; i<n_correlators; i++) {
		dcorr_tmp[i] = (double**)malloc(sizeof(double*)*nm);
		dcorr_ave[i] = (double**)malloc(sizeof(double*)*nm);
		dcorr_err[i] = (double**)malloc(sizeof(double*)*nm);
		dcorr_histo[i] = (double***)malloc(sizeof(double**)*nm);
		dcorr_tmp[i][0] = (double*)malloc(sizeof(double)*nm*T);
		dcorr_ave[i][0] = (double*)malloc(sizeof(double)*nm*T);
		dcorr_err[i][0] = (double*)malloc(sizeof(double)*nm*T);
		dcorr_histo[i][0] = (double**)malloc(sizeof(double*)*nm*T);
		dcorr_histo[i][0][0] = (double*)malloc(sizeof(double)*nm*T*nmeasures);
		for(n=0;n<nm;n++) {
			dcorr_tmp[i][n] = dcorr_tmp[i][0] + n*T;
			dcorr_ave[i][n] = dcorr_ave[i][0] + n*T;
			dcorr_err[i][n] = dcorr_err[i][0] + n*T;
			dcorr_histo[i][n] = dcorr_histo[i][0] + n*T;
			for(t = 0; t < T; t++) {
			   dcorr_histo[i][n][t] = dcorr_histo[i][0][0] + (n*T+t)*nmeasures;
			}
		}
	}
	
	m[0] = mass;

  log=freopen("check13.log","w",stdout);
  printf("\n");
  printf("Mesons (free case)\n");
  printf("---------------------------------------\n\n");
  printf("The lattice size is %dx%d^3\n",T,L);
  printf("size of the gluon rep: %d, size of the fermion rep: %d\n",NG,NF);
  printf("mass of the fermion: %f\n",mass);
  
  rlxd_init(1,12345);

	logger_setlevel(0,10015);

  geometry_eo_lexi();
  u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f();
#endif
  represent_gauge_field();

  set_spinor_len(VOLUME);

/* CALCOLO ESPLICITO */
  free_correlators(triplets, singlets);
  for(i=0; i<8; i++) {
  	for(t = 0; t < T; t++) {
  	   triplets[i][t] *= gb[i];
  	   singlets[i][t] *= gb[i];
  	}
  }

/* MESONI DUBLINESI */	
	for(counter = 0; counter < nmeasures; counter++) {
      printf("Counter = %d\n", counter);
   	dublin_meson_correlators(dcorr_tmp, dcorr_name, n_correlators, nm, m);
	   for(i=0; i<n_correlators; i++) {
		   for(n=0;n<nm;n++) {
		   	for(t = 0; t < T; t++)
		   	   dcorr_histo[i][n][t][counter] = dcorr_tmp[i][n][t];
		   }
	   }
	}
	
	
	
	for(i=0; i<n_correlators; i++) {
	  for(n=0;n<nm;n++) {
	  	for(t = 0; t < T; t++)
	  	   stat(nmeasures, dcorr_histo[i][n][t], dcorr_ave[i][n]+t, dcorr_err[i][n]+t);
	  }
	}

/* MESONI CON PROPAGATORE POINT-TO-ALL */

/*
	error((propfile = fopen("check13.prop", "wb"))==NULL,1,"Main",
			"Failed to open propagator file for writing");
	quark_propagator_QMR(propfile,0,nm,m,1.e-9);
	lprintf("MAIN",0,"MVM for last propagator: %ld\n",getMVM());

	fclose(propfile);

	error((propfile = fopen("check13.prop", "rb"))==NULL,1,"Main",
			"Failed to open propagator file for reading");
*/
	
	test = alloc_spinor_field_f();
	quark_prop = (suNf_spinor**)malloc(sizeof(suNf_spinor*)*4*NF);
	for(n = 0; n < 4*NF; n++)
		quark_prop[n] = alloc_spinor_field_f();

	propsize=sizeof(suNf_spinor)*VOLUME;

	for(n = 0; n < 4*NF; n++) {
		quark_propagator(n, nm, m, &(quark_prop[n]), 1e-9);

/*
		error((fread(quark_prop[n],(size_t) sizeof(suNf_spinor),
			(size_t)(VOLUME),propfile)!=(VOLUME))&&(!feof(propfile)),
			1,"Main", "Failed to read form quark propagator!");
*/		
		norm = spinor_field_sqnorm_f(quark_prop[n]);
		lprintf("MAIN",10,"Propagator norm [n=%d] = %e\n",n,norm);
		
		g5Dphi(m[0],test,quark_prop[n]);
		norm = 0.;
		for(x = 0; x < VOLUME; x++) {
			double *tf = (double*)(test+x);
			for(i = 0; i < sizeof(suNf_spinor)/sizeof(double); i++) {
				if(x == 0 && i == 2*n) norm = (tf[i]-1.)*(tf[i]-1.);
				else norm += tf[i]*tf[i];
			}
		}
		norm = sqrt(norm);
		lprintf("MAIN",10,"Testing inversion [n=%d] = %e\n",n,norm);

	}
	
	id_correlator(ptacorr[A], quark_prop);
	g0_correlator(ptacorr[X], quark_prop);
	g5_correlator(ptacorr[Pi], quark_prop);
	g0g5_correlator(ptacorr[Pi2], quark_prop);
	g1_correlator(ptacorr[Rho], quark_prop);
	g0g1_correlator(ptacorr[Rho2], quark_prop);
	g5g1_correlator(ptacorr[Y], quark_prop);
	g0g5g1_correlator(ptacorr[B], quark_prop);
/* 	CORR(g5_g0g5_correlator_im); */

/* 	fclose(propfile); */


/* STAMPA */
   
/* 	error = 0.; */

	printf("\nANALITICO\tPOINT-TO-ALL\tALL_TO_ALL\n");
	for(n = 0; n < 8; n++) {
		for(i = 0; i < n_correlators; i++)
			if(strcmp(dcorr_name[i],nameT[n])==0) break;

		printf("CORRELATOR TRIPLET %s\n", nameT[n]);
		if(i < n_correlators) {
			for(t = 0; t < T; t++)
				printf("%.8f\t%.8f\t%.8f (%.8f)\n",
							triplets[n][t],
							ptacorr[n][t],
							dcorr_ave[i][0][t], dcorr_err[i][0][t]);
		} else {
			for(t = 0; t < T; t++)
				printf("%.8f\t%.8f\n",
							triplets[n][t],
							ptacorr[n][t]);
		}
	}


	printf("\nANALITICO\tALL_TO_ALL\n");
	for(n = 0; n < 8; n++) {
		for(i = 0; i < n_correlators; i++)
			if(strcmp(dcorr_name[i],nameS[n])==0) break;

		printf("CORRELATOR SINGLET %s\n", nameS[n]);
		if(i < n_correlators) {
			for(t = 0; t < T; t++)
				printf("%.8f\t%.8f (%.8f)\n",
							singlets[n][t],
							dcorr_ave[i][0][t], dcorr_err[i][0][t]);
		} else {
				printf("%.8f\n",
							singlets[n][t]);
		}
	}

/*   error = sqrt(error); */
/*   printf("Error = %f\n", error); */


  printf("\n");
  fclose(log);
  
  exit(0);
}


double square(double x) {
   return x*x;
}


double re_ev(double k[4]) {
   return mass+4.0+cos((2.0*M_PI*k[0])/T)+cos((2.0*M_PI*k[1])/L)+cos((2.0*M_PI*k[2])/L)+cos((2.0*M_PI*k[3])/L);
}


double im_ev(double k[4]) {
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
   int i, t;
   double k[4];
   double sigma = 0.;
#ifdef ANTIPERIODIC_BC_T
   sigma = .5;
#endif
   
   lprintf("FREE_CORRELATORS",0,"sigma = %f\n",sigma);
   Abar = 0.;
   for(t = 0; t < T; t++) {
      A2[t] = 0.;
      B2[0][t] = B2[1][t] = 0.;
   }
   
   for(k[1] = 0.; k[1] < L-.5; k[1] += 1.)
   for(k[2] = 0.; k[2] < L-.5; k[2] += 1.)
   for(k[3] = 0.; k[3] < L-.5; k[3] += 1.) {
      
      
      for(t = 0; t < T; t++) {
         A[t] = 0.;
         B[0][t] = B[1][t] = 0.;
      }
      for(k[0] = sigma; k[0] < T+sigma-.5; k[0] += 1.) {
         re = re_ev(k);
         im = im_ev(k);
         norm2 = re*re+im*im;
         
         Abar += re / norm2;

         for(t = 0; t < T; t++) {
            ct = cos((2.0*M_PI*t*k[0])/T);
            st = sin((2.0*M_PI*t*k[0])/T);
            
            A[t] += re * ct / norm2;
            B[0][t] += sin((2.0*M_PI*k[0])/T) * st / norm2;
            B[1][t] += sin((2.0*M_PI*k[1])/L) * ct / norm2;
         }
      }
      for(t = 0; t < T; t++) {
         A2[t] += A[t] * A[t];
         B2[0][t] += B[0][t] * B[0][t];
         B2[1][t] += B[1][t] * B[1][t];
      }
      
   }
   
   lprintf("FREE_CORRELATORS",0,"Abar = %f\n",Abar);
   A2bar = Abar * Abar;
   lprintf("FREE_CORRELATORS",10,"A2[0] = %f\n",A2[0]);
   lprintf("FREE_CORRELATORS",10,"B2[0][0] = %f\n",B2[0][0]);
   lprintf("FREE_CORRELATORS",10,"B2[1][0] = %f\n",B2[1][0]);
   
   for(i = 0; i < 8; i++) {
      lprintf("FREE_CORRELATORS",10,"i=%d  g0=%d  gs=%d\n",i,g0[i],gs[i]);
      for(t = 0; t < T; t++) {
         triplets[i][t] = - 4*NF*g0[i] * (A2[t] - g0[i]*B2[0][t] - gs[i]*B2[1][t]) / (VOLUME*VOLUME);
         singlets[i][t] = triplets[i][t] + 16*p[i] * NF*NF* A2bar / (VOLUME*VOLUME);
      }
   }

   lprintf("FREE_CORRELATORS",0,"Exact free correlators computed.\n");
   
}



void stat(int ndata, double *data, double *ave, double *error) {
	int n, bin_size, nbins, i, k;
	double tmp;
	double* bin_ave = (double*)malloc(sizeof(double)*ndata);

	*ave = 0.0;
	for(n = 0; n < ndata; n++) *ave += data[n];
	*ave /= ndata;
	
	*error = 0.0;
	for(bin_size = 1; bin_size <= ndata; bin_size ++){
		nbins = ndata/bin_size;
		for(i = 0; i < nbins; i++){
			bin_ave[i] = 0.0;
			for(k = 0; k < bin_size; k++)
				bin_ave[i] += data[k+i*bin_size];
			bin_ave[i] /= bin_size;
		}
		tmp = 0.0;
		for(i = 0; i < nbins; i++)
			tmp += ( bin_ave[i] - *ave ) * ( bin_ave[i] - *ave );
		tmp = sqrt(tmp)/nbins;
		if(tmp > *error) *error = tmp;
	}
	
	free(bin_ave);
}


