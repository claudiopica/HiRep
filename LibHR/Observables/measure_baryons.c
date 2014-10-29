/***************************************************************************\
* Copyright (c) 2014 Vincent Drach, Ari Hietanen, Martin Hansen             *
*                                                                           *
\***************************************************************************/

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

#define corr_left_mult_g5(C) \
do { \
for(int i = 2; i < 4; i++) \
for(int j = 0; j < 4; j++) \
{ \
C[i][j].re = -C[i][j].re; \
C[i][j].im = -C[i][j].im; \
} \
} while(0)

#define corr_right_mult_g5(C) \
do { \
for(int i = 0; i < 4; i++) \
for(int j = 2; j < 4; j++) \
{ \
C[i][j].re = -C[i][j].re; \
C[i][j].im = -C[i][j].im; \
} \
} while(0)


// The function computes: out = Cg_\mu * in * g0(Cg_\nu)^\dagger*g0
void propagator_mul_left_right(suNf_propagator *out, suNf_propagator *in, int mu, int nu)
{
	suNf_propagator tmp;

	switch(mu)
	{
		case 0:
			_g2_propagator(tmp, *in);
			_propagator_mul_assign(tmp, -1);
			break;
		case 1:
			_g5g3_propagator(tmp, *in);
			_propagator_mul_assign(tmp, -1);
			break;
		case 2:
			_g0_propagator(tmp, *in);
			break;
		case 3:
			_g5g1_propagator(tmp, *in);
			break;
		default:
			break;
	}

	switch(nu)
	{
		case 0:
			_propagator_g2(*out, tmp);
			break;
		case 1:
			_propagator_g5g3(*out, tmp);
			break;
		case 2:
			_propagator_g0(*out, tmp);
			break;
		case 3:
			_propagator_g5g1(*out, tmp);
			_propagator_mul_assign(*out, -1);
			break;
		default:
			break;
	}
}

// baryons 1 ; C[beta][delta]  = S^ii' StildeT^jj' S^ kk'  
void _propagator_baryon1_mul2(complex C[4][4],suNf_propagator S,suNf_propagator Stilde,int i,int ip,int j,int jp,int k,int kp)
{
	complex tmp[4][4];

	for(int alpha = 0; alpha < 4; ++alpha)
	{
		for(int beta = 0;beta < 4; ++beta)
		{
			_complex_0(tmp[alpha][beta]);
			for(int gamma = 0; gamma < 4; gamma++)
			{
				_complex_mul_assign(tmp[alpha][beta], S.c[i].c[alpha].c[gamma].c[ip], Stilde.c[j].c[beta].c[gamma].c[jp]);
			}
		}
	} 

	for(int alpha = 0; alpha < 4; ++alpha)
	{
		for(int beta = 0; beta < 4; ++beta)
		{
			_complex_0(C[alpha][beta]);
			for(int gamma = 0; gamma < 4; gamma++)
			{
				_complex_mul_assign(C[alpha][beta], tmp[alpha][gamma], S.c[k].c[gamma].c[beta].c[kp]);
			}
		}
	}
}

// baryons 2 ; C[beta][delta]  = S^ii'_*tr(_S^ jj' StilteT^kk'
void _propagator_baryon2_mul2(complex C[4][4],suNf_propagator S,suNf_propagator Stilde,int i,int ip,int j,int jp,int k,int kp)
{
	complex tmp;
	_complex_0(tmp);

	for(int alpha = 0; alpha < 4; ++alpha)
	{
		for(int gamma = 0; gamma < 4; gamma++)
		{
			_complex_mul_assign(tmp, S.c[j].c[alpha].c[gamma].c[jp], Stilde.c[k].c[alpha].c[gamma].c[kp]);
		}
	}

	for(int alpha = 0; alpha < 4; ++alpha)
	{
		for(int beta = 0; beta < 4; ++beta)
		{
			_complex_mul(C[alpha][beta], S.c[i].c[alpha].c[beta].c[ip], tmp);
		}
	}
}

void contract_baryons(spinor_field* psi0, int tau)
{
	int ix, tc, idx;

	suNf_propagator S, Stmp;
	suNf_propagator Snucleon[2][2]; // Stilde for the nucleon
	suNf_propagator Sdelta[4][4]; // Stilde for the delta

	complex C1[4][4]; // to initialize
	complex C2[4][4]; // to initialize

	complex corr_nucleon[GLB_T][2][2][4][4];
	complex corr_delta[GLB_T][4][4][4][4];

	double col_factor[NF][NF][NF][NF][NF][NF];
	double eps[NF][NF][NF];

	lprintf("contract_baryons", 50, "Performing baryon contraction ");
	error(NG != 3, 1, "contract_baryons [measure_baryons.c]", "The code does not work for that number of colors!\n");
	
#ifdef REPR_FUNDAMENTAL 

	for(int a = 0; a < NF; ++a)
	for(int b = 0; b < NF; ++b)
	for(int c = 0; c < NF; ++c)
	{
		eps[a][b][c] = 0;
	}

	eps[0][1][2] = 1.;
	eps[1][2][0] = 1.;
	eps[2][0][1] = 1.;
	eps[1][0][2] = -1.;
	eps[0][2][1] = -1.;
	eps[2][1][0] = -1.;

#elif REPR_SYMMETRIC

	// def of the color factor for the sextet 
	double eS[NF][NG][NG];

	for(int i = 0; i < NF; ++i)
	for(int j = 0; j < NG; ++j)
	for(int k = 0; k < NG; ++k)
	{
		eS[i][j][k] = 0;
	}

	eS[0][0][0] = 1.;
	eS[1][0][1] = 1./sqrt(2);
	eS[1][1][0] = 1./sqrt(2);
	eS[2][1][1] = 1.;
	eS[3][2][0] = 1./sqrt(2);
	eS[3][0][2] = 1./sqrt(2);
	eS[4][1][2] = 1./sqrt(2);
	eS[4][2][1] = 1./sqrt(2);
	eS[5][2][2] = 1.;


	double tr_eS[6];
	double tmp1,tmp2,tmp3,tmp4,tmp5;
	double tmpmat1[3][3];
	double tmpmat2[3][3];

	for(int k = 0; k < NF; ++k)
	{ 
		tr_eS[k] = eS[k][0][0] + eS[k][1][1] + eS[k][2][2];
	}

	for(int a = 0; a < NF; ++a)
	for(int b = 0; b < NF; ++b)
	for(int c = 0; c < NF; ++c)
	{
		tmp1 = tr_eS[a] * tr_eS[b] * tr_eS[c];

		if(b == c)
		{
			tmp2 = -tr_eS[a];
		}
		else
		{
			tmp2 = 0;
		}

		if(a == c)
		{
			tmp3 = -tr_eS[b];
		}
		else
		{
			tmp3 = 0;
		}

		if(a == b)
		{
			tmp4 = -tr_eS[c];
		}
		else
		{
			tmp4 = 0;
		}

		for(int k = 0; k < 3; ++k)
		for(int l = 0; l < 3; ++l)
		{
			tmpmat1[k][l] = 0.;
			for(int m = 0; m < 3; ++m)
			{
				tmpmat1[k][l] += eS[b][k][m]*eS[c][m][l] +  eS[c][k][m]*eS[b][m][l];
			}
		}

		for(int k = 0; k < 3; ++k)
		for(int l = 0; l < 3; ++l)
		{
			tmpmat2[k][l] = 0.;
			for(int m = 0; m < 3; ++m)
			{
				tmpmat2[k][l] += eS[a][k][m]*tmpmat1[m][l];
			}
		}

		tmp5 = tmpmat2[0][0] +  tmpmat2[1][1] + tmpmat2[2][2];
		eps[a][b][c] = tmp1 + tmp2 + tmp3 + tmp4 + tmp5;

		if(fabs(eps[a][b][c]) < 1e-15)
		{
			eps[a][b][c] = 0;
		}

		lprintf("COLOR_FACTOR", 0, " %i %i %i %3.10e \n",a,b,c,eps[a][b][c]);
	}

#else 
	error(1, 1, "contract_baryons [measure_baryons.c]", "The code does not work for that representation\n");
#endif

	for(int a = 0; a < NF; ++a)
	for(int ap = 0; ap < NF; ++ap)
	for(int b = 0; b < NF; ++b)
	for(int bp = 0; bp < NF; ++bp)
	for(int c = 0; c < NF; ++c)
	for(int cp = 0; cp < NF; ++cp)
	{
		col_factor[a][b][c][ap][bp][cp] = eps[a][b][c]*eps[ap][bp][cp];
	}

	// Zero variables
	memset(corr_nucleon, 0, sizeof(corr_nucleon));
	memset(corr_delta, 0, sizeof(corr_delta));

	// LOOP VOLUME 
	for(int t = 0; t < T; t++)
	{
		tc = (zerocoord[0] + t + GLB_T - tau) % GLB_T;

		for(int x = 0; x < X; x++)
		for(int y = 0; y < Y; y++)
		for(int z = 0; z < Z; z++)
		{
			ix = ipt(t, x, y, z);

			// get S^ab_beta_gamma
			for(int a = 0; a < NF; ++a)
			for(int beta = 0; beta < 4; beta++)
			{
				idx = beta + a*4;
				_propagator_assign(Stmp, *_FIELD_AT(&psi0[idx],ix), a, beta);
			}

			// Why do we need to tranpose here?
			_propagator_transpose(S, Stmp);

			// Stilde = Gamma S Gamma_tilde
			// C  =  g0 g2
			// C g5 =  g1 g3
			// Nucleon case: Stilde = -Cg5 S Cg5 = - g5g0g2 S g5g0g2
 			// we forgot about the minus sign that will contribute only to a global sign...

			// Nucleon propagator
			_g5g0g2_propagator(Stmp, S);
			_propagator_g5g0g2(Snucleon[0][0], Stmp);

			_g5g0g2_propagator(Stmp, S);
			_propagator_g0g2(Snucleon[0][1], Stmp);

			_g0g2_propagator(Stmp, S);
			_propagator_g5g0g2(Snucleon[1][0], Stmp);

			_g0g2_propagator(Stmp, S);
			_propagator_g0g2(Snucleon[1][1], Stmp);




			// Delta propagator
			for(int mu = 0; mu < 4; mu++)
			for(int nu = 0; nu < 4; nu++)
			{
				propagator_mul_left_right(&Sdelta[mu][nu], &S, mu, nu);
			}

			for(int a = 0; a < NF; ++a)
			for(int ap = 0; ap < NF; ++ap)
			for(int b = 0; b < NF; ++b)
			for(int bp = 0; bp < NF; ++bp)
			for(int c = 0; c < NF; ++c)
			for(int cp = 0; cp < NF; ++cp)
			{
				if(fabs(col_factor[a][b][c][ap][bp][cp]) < 1e-5)
				{
					continue;
				}

				// Nucleon
				for (int k = 0; k < 2; k++)
				for (int l = 0; l < 2; l++)
				{
				_propagator_baryon1_mul2(C1, S, Snucleon[k][l], c, cp, b, bp, a, ap); // nucleon term 1
				_propagator_baryon2_mul2(C2, S, Snucleon[k][l], c, ap, a, cp, b, bp); // nucleon term 2
				
				if (k==1) corr_left_mult_g5(C1);
				if (k==1) corr_left_mult_g5(C2);

				if (l==1) corr_right_mult_g5(C1);
				if (l==1) corr_right_mult_g5(C2);

				// Multiply by color factor and accumulate in the correlator
					for(int i = 0; i < 4; i++)
					for(int j = 0; j < 4; j++)
					{
					corr_nucleon[tc][k][l][i][j].re -= col_factor[a][b][c][ap][bp][cp] * (C1[i][j].re - C2[i][j].re);
					corr_nucleon[tc][k][l][i][j].im -= col_factor[a][b][c][ap][bp][cp] * (C1[i][j].im - C2[i][j].im);
					}
			  }


				// Delta
				for(int mu = 0; mu < 4; mu++)
				for(int nu = 0; nu < 4; nu++)
				{
					_propagator_baryon1_mul2(C1, S, Sdelta[mu][nu], c, cp, b, bp, a, ap); // delta term 1
					_propagator_baryon2_mul2(C2, S, Sdelta[mu][nu], c, ap, a, cp, b, bp); // delta term 2

					// Multiply by color factor and accumulate in the correlator
					for(int i = 0; i < 4; i++)
					for(int j = 0; j < 4; j++)
					{
						corr_delta[tc][mu][nu][i][j].re -= 2*col_factor[a][b][c][ap][bp][cp] * (C2[i][j].re - 2*C1[i][j].re);
						corr_delta[tc][mu][nu][i][j].im -= 2*col_factor[a][b][c][ap][bp][cp] * (C2[i][j].im - 2*C1[i][j].im);
					}
				}
			} // end cp
		} // END SPATIAL LOOP
	} // END TEMPORAL LOOP

 global_sum((double*)corr_nucleon, 2*2*2*4*4*GLB_T);
 global_sum((double*)corr_delta, 2*4*4*4*4*GLB_T);

	// Print nucleon correlator
	for(int t = 0; t < GLB_T; t++)
	for(int k = 0; k < 2; k++)
	for(int l = 0; l < 2; l++)
	for(int i = 0; i < 4; i++)
	{

		lprintf("CORR_NUCLEON", 0, "%d %d %d %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e \n",
				t,k,l,
				corr_nucleon[t][k][l][i][0].re,
				corr_nucleon[t][k][l][i][0].im,
				corr_nucleon[t][k][l][i][1].re,
				corr_nucleon[t][k][l][i][1].im,
				corr_nucleon[t][k][l][i][2].re,
				corr_nucleon[t][k][l][i][2].im,
				corr_nucleon[t][k][l][i][3].re,
				corr_nucleon[t][k][l][i][3].im
		);
	}

	// Print delta correlator
	for(int t = 0; t < GLB_T; t++)
	for(int mu = 0; mu < 4; mu++)
	for(int nu = 0; nu < 4; nu++)
	for(int i = 0; i < 4; i++)
	{
		lprintf("CORR_DELTA", 0, "%d %d %d %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e \n",
				t,
				mu,
				nu,
				corr_delta[t][mu][nu][i][0].re,
				corr_delta[t][mu][nu][i][0].im,
				corr_delta[t][mu][nu][i][1].re,
				corr_delta[t][mu][nu][i][1].im,
				corr_delta[t][mu][nu][i][2].re,
				corr_delta[t][mu][nu][i][2].im,
				corr_delta[t][mu][nu][i][3].re,
				corr_delta[t][mu][nu][i][3].im
		);
	}

	lprintf("contract_baryons", 50, "Measuring DONE!\n");
}
