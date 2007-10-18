/******************************************************************************
*
* Test of modules
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

int nhb,nor,nit,nth,nms,level,seed;
double beta;

static double hmass=0.1;


void D(suNf_spinor *out, suNf_spinor *in){
	Dphi(hmass,out,in);
}

void H(suNf_spinor *out, suNf_spinor *in){
	g5Dphi(hmass,out,in);
}

static suNf_spinor *tmp;
void M(suNf_spinor *out, suNf_spinor *in){
	g5Dphi(-hmass,tmp,in); 
	g5Dphi(-hmass,out,tmp);
}


int main(int argc,char *argv[])
{
	int i;
	double tau;
	suNf_spinor *s1, *s2;
	suNf_spinor **res;

	mshift_par par;
	MINRES_par MINRESpar2;

	int cgiters;

	printf("Gauge group: SU(%d)\n",NG);
	printf("Fermion representation: dim = %d\n",NF);
	geometry_eo_lexi();
	test_geometry();
	printf("The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
	printf("\n");

	level=1;
	seed=123;
	printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
	fflush(stdout);

	rlxd_init(level,seed);

	u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
	u_gauge_f=alloc_gfield_f();
#endif
	represent_gauge_field();

	printf("Generating a random gauge field... ");
	fflush(stdout);
	random_u();
	printf("done.\n");
	represent_gauge_field();

	set_spinor_len(VOLUME);

	par.n = 6;
	par.shift=(double*)malloc(sizeof(double)*(par.n));
	par.err2=1.e-28;
	par.max_iter=0;
	res=(suNf_spinor**)malloc(sizeof(suNf_spinor*)*(par.n));
	res[0]=alloc_spinor_field_f(par.n+3);
	for(i=1;i<par.n;++i)
		res[i]=res[i-1]+VOLUME;
	s1=res[par.n-1]+VOLUME;
	s2=s1+VOLUME;
	tmp=s2+VOLUME;

	par.shift[0]=+0.1;
	par.shift[1]=-0.21;
	par.shift[2]=+0.05;
	par.shift[3]=-0.01;
	par.shift[4]=-0.15;
	par.shift[5]=-0.05;

	gaussian_spinor_field(s1);

	/* TEST MINRES_M */

	printf("\n");
	printf("Testing MINRES multishift\n");
	printf("-------------------------\n");

	cgiters=MINRES_mshift(&par, &H, s1, res);
	printf("Converged in %d iterations\n",cgiters);

	for(i=0;i<par.n;++i){
		H(s2,res[i]);
		spinor_field_mul_add_assign_f(s2,-par.shift[i],res[i]);
		spinor_field_sub_assign_f(s2,s1);
		tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
		printf("test MINRES[%d] = %e (req. %e)\n",i,tau,par.err2);
	}

	/* TEST MINRES_M */

	printf("\n");
	printf("Testing MINRES \n");
	printf("-------------- \n");

	MINRESpar2.err2=par.err2;
	MINRESpar2.max_iter=0;

	cgiters=MINRES(&MINRESpar2, &H, s1, res[0],0);
	for(i=1;i<par.n;++i){
		hmass=0.1-par.shift[i-1];
		cgiters+=MINRES(&MINRESpar2, &H, s1, res[i],res[i-1]);
	}
	printf("Converged in %d iterations\n",cgiters);

	hmass=0.1;
	for(i=0;i<par.n;++i){
		if(i!=0)
			hmass=0.1-par.shift[i-1];
		H(s2,res[i]);
		spinor_field_sub_f(s2,s2,s1);
		tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
		printf("test MINRES[%d] = %e (req. %e)\n",i,tau,par.err2);
	}

	free_field(res[0]);
	free(res);
	free(par.shift);

	exit(0);
}
