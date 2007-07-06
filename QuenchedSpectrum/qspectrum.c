
/*******************************************************************************
*
* Computation of quenched quark propagator
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

float mass;

void read_cmdline(int argc, char*argv[])
{
	int i, ai=0, ao=0;
	FILE *in=NULL, *out=NULL;
	for (i=1;i<argc;i++){
		if (strcmp(argv[i],"-i")==0)
			ai=i+1;
		else if (strcmp(argv[i],"-o")==0)
			ao=i+1;
	}

	if (ao!=0)
		out=freopen(argv[ao],"w",stdout);
	else
		out=stdout;

	error(ai==0,1,"suN.c",
			"Syntax: suN -i <input file> [-o <output file>]");

	in=freopen(argv[ai],"r",stdin);
	error(in==NULL,1,"run1.c","Cannot open input file");

	scanf("beta %f nhb %d nor %d nit %d nth %d nms %d level %d seed %d",
			&beta,&nhb,&nor,&nit,&nth,&nms,&level,&seed);
	fclose(in);

}

int main(int argc,char *argv[])
{
	int i,n;

	char propname[256];
	FILE *propfile;
	float *m;

	float kappa[]={0.160,0.161,0.162,0.163,0.164};
	int nm=sizeof(kappa)/sizeof(float);
	m = (float *) malloc(sizeof(kappa));
	for (i=0; i<nm; ++i)
	{
		m[i] = 0.5/kappa[i] - 4.0;
	}

	read_cmdline(argc, argv); 

	printf("Gauge group: SU(%d)\n",NG);
	printf("Fermion representation: dim = %d\n",NF);
	printf("The lattice size is %dx%d^3\n",T,L);
	printf("beta = %2.4f\n",beta);
	printf("nth  = %d\tNumber of thermalization cycles\n",nth);
	printf("nms  = %d\tNumber of measure cycles\n",nms);
	printf("nit  = %d\tNumber of hb-or iterations per cycle\n",nit);
	printf("nhb  = %d\tNumber of heatbaths per iteration\n",nhb);   
	printf("nor  = %d\tNumber or overrelaxations per iteration\n",nor);
	printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
	printf("Computing quark prop for %d masses: ",nm);
	for(i=0;i<nm;++i)
		printf("%e ",m[i]);
	printf("corresponding to the following kappa:");
	for(i=0;i<nm;++i)
		printf("%e ",kappa[i]);
	printf("\n\n");

	fflush(stdout);

	rlxs_init(level,seed);

	geometry_eo_lexi();
	/*geometry_blocked();*/
	test_geometry();

	u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
	u_gauge_f=alloc_gfield_f();
#endif

	represent_gauge_field();

	sprintf(propname,"quark_prop_qmr_%3.5f_%d_%d_%d_%d",beta,T,L,L,L);
	error((propfile = fopen(propname, "wb"))==NULL,1,"Main",
			"Failed to open propagator file for writing");
	fwrite(&nm,(size_t) sizeof(int),1,propfile);
	for(i=0;i<nm;++i)
		fwrite(m+i,(size_t) sizeof(float),1,propfile);

	/* Termalizzazione */
	for (i=0;i<nth;++i) {
		update(beta,nhb,nor);
		if ((i%10)==0)
			printf("%d",i);
		else
			printf(".");
		fflush(stdout);
	}
	if(i) printf("%d\nThemalization done.\n",i);
	represent_gauge_field();

	/* Misure */
	for (i=0;i<nms;++i){ /* nms misure */
		printf("[%d] <p> = %1.6f\n",i,avr_plaquette());
		fflush(stdout);
		quark_propagator_QMR(propfile,0,nm,m);

		for (n=0;n<nit;n++) /* nit updates */
			update(beta,nhb,nor);
		represent_gauge_field();

	}


	fclose(propfile);


	free_field(u_gauge);

#ifndef REPR_FUNDAMENTAL
	free_field(u_gauge_f);
#endif

	return 0;
}
