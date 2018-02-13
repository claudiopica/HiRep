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




void measure_SUS(int tau){
	complex SUS[GLB_T][4];
	double SS[GLB_T];
	memset(SUS, 0, sizeof(SUS));
	memset(SS, 0, sizeof(SS));
	suNg_vector *Sx, *Sup, UtimesSup;
	complex SUSup;
	double SStmp;
	suNg *U;
	for(int t = 0; t < T; t++)
	{
		int tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
		for(int x = 0; x < X; x++) for(int y = 0; y < Y; y++) for(int z = 0; z < Z; z++)
		{
			int ix = ipt(t, x, y, z);

			Sx = _FIELD_AT(u_scalar,ix);
			_vector_prod_re_g(SStmp, *Sx, *Sx);
			SS[tc] += SStmp;

			for(int mu = 0; mu < 4; mu++)
			{
				Sup = _FIELD_AT(u_scalar, iup(ix,mu));
				U = _4FIELD_AT(u_gauge, ix, mu);
				_suNg_multiply(UtimesSup, *U, *Sup);
				_vector_prod_re_g(SUSup.re, *Sx, UtimesSup);
				_vector_prod_im_g(SUSup.im, *Sx, UtimesSup);
				_complex_add_assign(SUS[tc][mu], SUSup);
			}
		}
	}
	global_sum((double*)SUS,8*GLB_T);
	global_sum((double*)SS,GLB_T);
	for(int t = 0; t < GLB_T; t++)
	{
		int tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
		lprintf("Hi!",0,"SS[%d] = %f\n",tc, SS[tc]);
	}
	for(int t = 0; t < GLB_T; t++)
	{
		int tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;
		for(int mu = 0; mu < 4; mu++)
		{
			lprintf("Hi!",0,"SUS[%d][%d] = %f %f\n",tc,mu, SUS[tc][mu].re, SUS[tc][mu].im);
		}
	}
}


complex average_SUS(){

	suNg_vector *Sx, *Sup, UtimesSup;
	complex SUSup, av_SUS;
	suNg *U;

	_complex_0(av_SUS);

	_MASTER_FOR(&glattice,ix){
		Sx = _FIELD_AT(u_scalar,ix);
		for(int mu = 0; mu < 4; mu++)
		{
			Sup = _FIELD_AT(u_scalar, iup(ix,mu));
			U = _4FIELD_AT(u_gauge, ix, mu);
			_suNg_multiply(UtimesSup, *U, *Sup);
			_vector_prod_re_g(SUSup.re, *Sx, UtimesSup);
			_vector_prod_im_g(SUSup.im, *Sx, UtimesSup);
			_complex_add_assign(av_SUS, SUSup);
		}
	}
	global_sum((double*)&av_SUS,sizeof(av_SUS)/sizeof(double));

	_complex_mulr(av_SUS,1/(4.0*GLB_VOLUME),av_SUS);

	return av_SUS;
}


suNg_vector average_S(){

	suNg_vector S_av, Sc;
	_vector_zero_g(S_av);
	_MASTER_FOR(&glattice,ix){
		Sc = *pu_scalar(ix);
		_vector_add_assign_g(S_av,Sc);
	}
	global_sum((double*)&S_av,sizeof(S_av)/sizeof(double));

	_vector_mul_g(S_av,1/(double)(GLB_VOLUME),S_av);

	return S_av;
}


double average_SdagS(){

	suNg_vector Sc;
	double S_cond = 0;
	double S_sq = 0;
	_MASTER_FOR_SUM(&glattice,ix,S_cond){
		Sc = *pu_scalar(ix);
		_vector_prod_re_g(S_sq,Sc,Sc);
		S_cond += S_sq;
	}
	global_sum(&S_cond,1);

	return S_cond/(double)GLB_VOLUME;
}
