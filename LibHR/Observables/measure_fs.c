/*
 * Measures the scalar-fermion meson two-point function
 */

#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"
#include "io.h"
#include "random.h"
#include "communications.h"
#include "ranlux.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "propagator.h"
#include <string.h>
#include "meson_observables.h"
#define PI 3.141592653589793238462643383279502884197

static void scalar_prop_scalar(complex result[4][4], suNg_vector *S_src, suNg_vector *S_snk, suNf_propagator* psi ){
	suNg_vector S_dag;
	vector_star(&S_dag, S_snk);
	for(int mu=0; mu<4;mu++)
		for(int nu=0; nu<4;nu++)
		{
			for(int i=0;i<NF;i++)
				for(int j=0;j<NF;j++)
				{
					complex c1;
					_complex_mul(c1,(*S_src).c[i], S_dag.c[j]);
					_complex_mul_assign(result[mu][nu], (*psi).c[j].c[nu].c[mu].c[i],c1) ;
				}
		}
}

//Measures meson correlators made of 1 scalar and 1 fermion
//void contract_fs(spinor_field* psi0, int tau){
//	complex corr_fs[GLB_T][4][4];
//	memset(corr_fs, 0, sizeof(corr_fs));
//	suNf_propagator Stmp;
//	suNg_vector *S_src = pu_scalar(ipt(tau,0,0,0)); //CHANGE LATER
//	for(int t = 0; t < T; t++)
//	{
//		int tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
//	//	for(int i=0;i<4;i++){
//	//		for(int j=0;j<4;j++){
//	//			_complex_0(corr_fs[tc][i][j]);
//	//		}
//	//	}
//
//		for(int x = 0; x < X; x++) for(int y = 0; y < Y; y++) for(int z = 0; z < Z; z++)
//		{
//			int ix = ipt(t, x, y, z);
//			suNg_vector *S_snk = pu_scalar(ix); //CHANGE LATER
//			for(int a = 0; a < NF; ++a)
//				for(int beta = 0; beta < 4; beta++)
//				{
//					int idx = beta + a*4;
//					_propagator_assign(Stmp, *_FIELD_AT(&psi0[idx],ix), a, beta);
//				//	lprintf("DEBUG SCALAR SPECTRUM",0, "%d %d %d %d %d %d %f %f %f %f %f %f %f %f \n", a, beta, t, x, y, z, Stmp.c[a].c[0].c[beta].c[a].re, Stmp.c[a].c[0].c[beta].c[a].im, Stmp.c[a].c[1].c[beta].c[a].re, Stmp.c[a].c[1].c[beta].c[a].im, Stmp.c[a].c[2].c[beta].c[a].re, Stmp.c[a].c[2].c[beta].c[a].im, Stmp.c[a].c[3].c[beta].c[a].re, Stmp.c[a].c[3].c[beta].c[a].im);
//				}
//			scalar_prop_scalar(corr_fs[tc], S_src, S_snk, &Stmp );
//		}
//	}
//	global_sum((double*)corr_fs, 2*4*4*GLB_T);
//	for(int t = 0; t < GLB_T; t++) for(int i = 0; i < 4; i++)
//	{
//
//		lprintf("CORR_FS", 0, "%d %d %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e \n",
//				t,i,
//				corr_fs[t][i][0].re,
//				corr_fs[t][i][0].im,
//				corr_fs[t][i][1].re,
//				corr_fs[t][i][1].im,
//				corr_fs[t][i][2].re,
//				corr_fs[t][i][2].im,
//				corr_fs[t][i][3].re,
//				corr_fs[t][i][3].im
//		       );
//	}
//}

void contract_fs(spinor_field* psi0, int tau){
	complex corr_fs[GLB_T][4][4];
	memset(corr_fs, 0, sizeof(corr_fs));
	suNf_propagator Stmp;
	suNg_vector *S_src;
	// this only works for tau=0! 
  	suNg_vector S_zero=*pu_scalar(ipt(0,0,0,0));
        tau=0;
        //broadcast from mpi process 0
        bcast((double *) &S_zero, sizeof(S_zero)/sizeof(double)); 
 	S_src=&S_zero;        

	for(int t = 0; t < T; t++)
	{
		int tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;

		for(int x = 0; x < X; x++) for(int y = 0; y < Y; y++) for(int z = 0; z < Z; z++)
		{
			int ix = ipt(t, x, y, z);
			suNg_vector *S_snk = pu_scalar(ix); 
			for(int a = 0; a < NF; ++a)
				for(int beta = 0; beta < 4; beta++)
				{
					int idx = beta + a*4;
					_propagator_assign(Stmp, *_FIELD_AT(&psi0[idx],ix), a, beta);
				}
			scalar_prop_scalar(corr_fs[tc], S_src, S_snk, &Stmp );
		}
	}
	global_sum((double*)corr_fs, 2*4*4*GLB_T);
	for(int t = 0; t < GLB_T; t++) for(int i = 0; i < 4; i++)
	{

		lprintf("CORR_FS", 0, "%d %d %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e \n",
				t,i,
				corr_fs[t][i][0].re,
				corr_fs[t][i][0].im,
				corr_fs[t][i][1].re,
				corr_fs[t][i][1].im,
				corr_fs[t][i][2].re,
				corr_fs[t][i][2].im,
				corr_fs[t][i][3].re,
				corr_fs[t][i][3].im
		       );
	}
}

void contract_s_left(spinor_field* psi0, int tau){
	complex corr_fs[GLB_T][4][4];
	memset(corr_fs, 0, sizeof(corr_fs));

	for(int t = 0; t < T; t++)
	{
		int tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;

		for(int x = 0; x < X; x++) for(int y = 0; y < Y; y++) for(int z = 0; z < Z; z++)
		{
			int ix = ipt(t, x, y, z);
			suNg_vector *S_snk = pu_scalar(ix); 
			for(int alpha = 0; alpha < 4; alpha++) //spin source
				for(int beta = 0; beta < 4; beta++) //spin sink
					for(int a = 0; a < NF; ++a) //color
					{
						_complex_prod_assign(corr_fs[tc][alpha][beta], (*S_snk).c[a],(*_FIELD_AT(&psi0[alpha],ix)).c[beta].c[a]);
					}
		}
	}
	global_sum((double*)corr_fs, sizeof(corr_fs)/sizeof(double));
	for(int t = 0; t < GLB_T; t++) for(int i = 0; i < 4; i++)
	{

		lprintf("CORR_FS_SSRC", 0, "%d %d %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e \n",
				t,i,
				corr_fs[t][i][0].re,
				corr_fs[t][i][0].im,
				corr_fs[t][i][1].re,
				corr_fs[t][i][1].im,
				corr_fs[t][i][2].re,
				corr_fs[t][i][2].im,
				corr_fs[t][i][3].re,
				corr_fs[t][i][3].im
		       );
	}
}

void measure_fs_pt(double* m, double precision){
	spinor_field* source = alloc_spinor_field_f(4*NF,&glattice); //This isn't glat_even so that the odd sites will be set to zero explicitly
	spinor_field* prop =  alloc_spinor_field_f(4*NF,&glattice);
	int nm=1;
	int tau=0;

	init_propagator_eo(nm, m, precision);

	// create point source
	create_full_point_source(source,tau);

	// calc invert
	calc_propagator(prop,source,4*NF);//4x3 for QCD 

	// perform contraction
	contract_fs(prop,tau);
  
	// free 
	free_propagator_eo(); 
	free_spinor_field_f(source);
	free_spinor_field_f(prop);
}

void measure_fs_scSrc(double* m, double precision){
	spinor_field* source = alloc_spinor_field_f(4,&glattice); //This isn't glat_even so that the odd sites will be set to zero explicitly
	spinor_field* prop =  alloc_spinor_field_f(4,&glattice);
	int nm=1;
	int tau=0;


	init_propagator_eo(nm, m, precision);

	// create point source
	create_scalar_source(source,tau);

	// calc invert
	calc_propagator(prop,source,4);//4x3 for QCD 

	// perform contraction
	contract_s_left(prop,tau);
  
	// free 
	free_propagator_eo(); 
	free_spinor_field_f(source);
	free_spinor_field_f(prop);
}

