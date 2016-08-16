/***************************************************************************
* Copyright (c) 2014, Martin Hansen                                        *
* All rights reserved.                                                     *
*                                                                          *
* Chronological inverter using the MRE algorithm                           *
* arXiv: hep-lat/9509012                                                   *
***************************************************************************/

#include "update.h"
#include "linear_algebra.h"
#include "memory.h"
#include "global.h"
#include "dirac.h"
#include "logger.h"
#include <math.h>

#define cabs(a) sqrt(a.re*a.re + a.im*a.im)
#define MAX 15

// Global variables
static spinor_field *v[MAX];
static spinor_field *Dv;
static int num_init = 0;

// Variables used in in the LU solver
static complex A[MAX][MAX];
static complex b[MAX];
static complex x[MAX];
static complex y[MAX];
static int mutate[MAX];

void gram_schmidt(mre_par *par, int p, int max)
{
	complex rij;
	double rii;

	for(int i = 0; i < max; i++)
	{
		spinor_field_copy_f(v[i], &par->s[p][i]);
	}

	for(int i = 0; i < max; i++)
	{
		rii = spinor_field_sqnorm_f(v[i]);
		rii = 1.0/sqrt(rii);
		spinor_field_mul_f(v[i], rii, v[i]);

		for(int j = i+1; j < max; j++)
		{
			rij = spinor_field_prod_f(v[i], v[j]);
			_complex_minus(rij, rij);
			spinor_field_mulc_add_assign_f(v[j], rij, v[i]);
		}
	}
}

void lu_solve(int max)
{
	double big;
	int row;
	complex ctmp;
	int itmp;

	// Setup mutate
	for(int i = 0; i < max; i++)
	{
		mutate[i] = i;
	}

	// LU factorization
	for(int i = 0; i < max; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			for(int k = 0; k < j; k++)
			{
				_complex_mul_sub_assign(A[j][i], A[j][k], A[k][i]);
			}
		}

		big = cabs(A[i][i]);
		row = i;

		for(int j = i+1; j < max; j++)
		{
			for(int k = 0; k < i; k++)
			{
				_complex_mul_sub_assign(A[j][i], A[j][k], A[k][i]);
			}

			if(cabs(A[j][i]) > big)
			{
				big = cabs(A[j][i]);
				row = j;
			}
		}

		if(big < 1.0e-14)
		{
			lprintf("MRE", 10, "LU decomposition failed: matrix is singular\n");
			return;
		}

		if(row != i)
		{
			for(int k = 0; k < max; k++)
			{
				ctmp = A[row][k];
				A[row][k] = A[i][k];
				A[i][k] = ctmp;
			}

			itmp = mutate[row];
			mutate[row] = mutate[i];
			mutate[i] = itmp;
		}

		for(int k = i+1; k < max; k++)
		{
			_complex_div(ctmp, A[k][i], A[i][i]);
			A[k][i] = ctmp;
		}
	}

	// Forward substitution
	for(int i = 0; i < max; i++)
	{
		y[i] = b[mutate[i]];
		y[i].im = -y[i].im;
		for(int k = 0; k < i; k++)
		{
			_complex_mul_sub_assign(y[i], A[i][k], y[k]);
		}
	}

	// Backward substitution
	for(int i = max-1; i >= 0; i--)
	{
		x[i] = y[i];
		for(int k = i+1; k < max; k++)
		{
			_complex_mul_sub_assign(x[i], A[i][k], x[k]);
		}
		_complex_div(ctmp, x[i], A[i][i]);
		x[i] = ctmp;
	}
}

int factorial(int k)
{
	int f = 1;
	for(; k > 1; k--) f *= k;
	return f;
}

int coefficient(int k, int n)
{
	int c = factorial(n)/(factorial(k)*factorial(n-k));
	if((k - 1) % 2) c = -c;
	return c;
}

void mre_init(mre_par *par, int max, double prec)
{
	if(max == 0)
	{
		par->max = 0;
		par->init = 0;
		return;
	}
	else
	{
		par->max = (max < MAX) ? max : MAX;
		par->init = 1;
	}

	par->s[0] = alloc_spinor_field_f(par->max, &glat_default);
	par->s[1] = alloc_spinor_field_f(par->max, &glat_default);
	par->num[0] = 0;
	par->num[1] = 0;

	if(num_init == 0)
	{
		Dv = alloc_spinor_field_f(1, &glat_default);
	}

	for(int i = num_init; i < par->max; i++)
	{
		v[i] = alloc_spinor_field_f(1, &glat_default);
		num_init++;
	}

	if(prec > 1e-14)
	{
		lprintf("MRE", 10,  "WARNING: Inverter precision should be at least 1e-14 to ensure reversibility!\n");
	}

	lprintf("MRE", 10, "Enabled chronological inverter with %d past solutions\n", par->max);
}

void mre_store(mre_par *par, int p, spinor_field *in)
{
	spinor_field *tmp;

	if(num_init == 0 || par->init == 0 || par->max <= 0)
	{
		return;
	}

	if(p > 1)
	{
		lprintf("MRE", 10, "Cannot store solution: p = %d is not valid\n", p);
		return;
	}

	// Shift all vectors by one and store the new vector at position zero
	tmp = &par->s[p][par->max - 1];

	for(int i = (par->max - 1); i > 0; i--)
	{
		par->s[p][i] = par->s[p][i-1];
	}

	par->s[p][0] = *tmp;
	spinor_field_copy_f(&par->s[p][0], in);
	par->num[p]++;
}

void mre_guess(mre_par *par, int p, spinor_field *out, spinor_operator D, spinor_field *pf)
{
	int max;

	if(num_init == 0 || par->init == 0 || par->max <= 0)
	{
		return;
	}

	if(p > 1)
	{
		lprintf("MRE", 10, "Cannot guess solution: p = %d is not valid\n", p);
		return;
	}

	spinor_field_zero_f(out);
	max = (par->num[p] > par->max) ? par->max : par->num[p];
	gram_schmidt(par, p, max);

	for(int i = 0; i < max; i++)
	{
		D(Dv, v[i]);

		for(int j = 0; j < max; j++)
		{
			A[j][i] = spinor_field_prod_f(Dv, v[j]);
		}

		b[i] = spinor_field_prod_f(v[i], pf);
	}

	lu_solve(max);

	for(int i = 0; i < max; i++)
	{
//		spinor_field_mul_add_assign_f(out, coefficient(i+1, max), s[p][i]);
		spinor_field_mulc_add_assign_f(out, x[i], v[i]);
	}
}
