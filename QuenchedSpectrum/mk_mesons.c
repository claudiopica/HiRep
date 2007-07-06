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

#define CORR(name) \
				name(tmpcorr, quark_prop);\
				printf("[%d] mass=%2.6f " #name "= ",i,m[k]);\
				for(n=0;n<T;++n) {\
					printf("%e ",tmpcorr[n]);\
				}\
				printf("\n");\
				fflush(stdout)

int main(int argc,char *argv[])
{
	int i,k,n;

	suNf_spinor **quark_prop;
	FILE *propfile;
	long int start;
	long int propsize;
	long int cur_pos;

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
	quark_prop=(suNf_spinor**)malloc(sizeof(suNf_spinor*)*4*NF);
	for (n=0;n<4*NF;++n)
		quark_prop[n]=(suNf_spinor*)malloc(sizeof(suNf_spinor)*VOLUME);

	propsize=sizeof(suNf_spinor)*VOLUME;
	start=ftell(propfile);

	i=0;
	/* Misure */
	do {

		float tmpcorr[T];

		for (k=0;k<nm;++k){
			/* read propagator for 1 mass from file */
			for (n=0;n<4*NF;++n){
				cur_pos=start+n*nm*propsize;
				error(fseek(propfile,cur_pos,SEEK_SET)!=0,1,"Main","Fseek failed!");
				error((fread(quark_prop[n],(size_t) sizeof(suNf_spinor),
								(size_t)(VOLUME),propfile)!=(VOLUME))&&(!feof(propfile)),
						1,"Main", "Failed to read form quark propagator!");
			}
			if (!feof(propfile)){
				CORR(id_correlator);
				CORR(g0_correlator);
				CORR(g5_correlator);
				CORR(g0g5_correlator);
				CORR(g1_correlator);
				CORR(g2_correlator);
				CORR(g3_correlator);
				CORR(g0g1_correlator);
				CORR(g0g2_correlator);
				CORR(g0g3_correlator);
				CORR(g5g1_correlator);
				CORR(g5g2_correlator);
				CORR(g5g3_correlator);
				CORR(g0g5g1_correlator);
				CORR(g0g5g2_correlator);
				CORR(g0g5g3_correlator);
			}
			start+=propsize;
		}

		start+=(nm)*(4*NF-1)*propsize;
		++i;

	} while(!feof(propfile));

	fclose(propfile);


	for (n=0;n<4*NF;++n)
		free(quark_prop[n]);
	free(quark_prop);

	free(m);

	return 0;
}
