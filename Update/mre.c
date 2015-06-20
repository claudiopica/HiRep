/***************************************************************************\
* Copyright (c) 2014, Martin Hansen                                         *
* All rights reserved.                                                      *
\***************************************************************************/

#include "linear_algebra.h"
#include "memory.h"
#include "global.h"
#include "dirac.h"
#include <math.h>

#define cabs(a) \
	sqrt(a.re*a.re + a.im*a.im)

#define complex_mul_sub_assign(c,a,b) \
	c.re -= a.re*b.re - a.im*b.im; \
	c.im -= a.re*b.im + a.im*b.re; \

#define complex_div_assign(a,b) \
if(1) \
{ \
	double re = a.re; \
	double sq = b.re*b.re + b.im*b.im; \
	a.re = (re*b.re + a.im*b.im) / sq; \
	a.im = (a.im*b.re - re*b.im) / sq; \
}

#define MAX  5
#define REPS 2
static spinor_field *s[REPS][MAX];
static spinor_field *v[MAX];
static spinor_field *Dv;
static int initialized = 0;
static int num[REPS];

// Variables used in in the LU solver
static complex A[MAX][MAX];
static complex b[MAX];
static complex x[MAX];
static complex y[MAX];
static int mutate[MAX];

// Preconditioning
#ifdef UPDATE_EO
#define TYPE &glat_even
#else
#define TYPE &glattice
#endif

void gram_schmidt(int p, int max)
{
	complex rij;
	double rii;

	for(int i = 0; i < max; i++)
	{
		spinor_field_copy_f(v[i], s[p][i]);
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
				complex_mul_sub_assign(A[j][i], A[j][k], A[k][i]);
			}
		}

		big = cabs(A[i][i]);
		row = i;

		for(int j = i+1; j < max; j++)
		{
			for(int k = 0; k < i; k++)
			{
				complex_mul_sub_assign(A[j][i], A[j][k], A[k][i]);
			}

			if(cabs(A[j][i]) > big)
			{
				big = cabs(A[j][i]);
				row = j;
			}
		}

		if(big < 1.0e-16)
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
			complex_div_assign(A[k][i], A[i][i]);
		}
	}

	// Forward substitution
	for(int i = 0; i < max; i++)
	{
		y[i] = b[mutate[i]];
		y[i].im = -y[i].im;
		for(int k = 0; k < i; k++)
		{
			complex_mul_sub_assign(y[i], A[i][k], y[k]);
		}
	}

	// Backward substitution
	for(int i = max-1; i >= 0; i--)
	{
		x[i] = y[i];
		for(int k = i+1; k < max; k++)
		{
			complex_mul_sub_assign(x[i], A[i][k], x[k]);
		}
		complex_div_assign(x[i], A[i][i]);
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
	int c;
	c = factorial(n)/(factorial(k)*factorial(n-k));
	if((k - 1) % 2) c = -c;
	return c;
}

void init()
{
	for(int i = 0; i < MAX; i++)
	{
		for(int j = 0; j < REPS; j++)
		{
			s[j][i] = alloc_spinor_field_f(1, TYPE);
			spinor_field_zero_f(s[j][i]);
			num[j] = 0;
		}

		v[i] = alloc_spinor_field_f(1, TYPE);
		spinor_field_zero_f(v[i]);
	}

	Dv = alloc_spinor_field_f(1, TYPE);
	initialized = 1;
}

void mre_store(int p, spinor_field *in)
{
	spinor_field *tmp;

	if(!initialized)
	{
		init();
	}

	if(p >= REPS)
	{
		lprintf("MRE", 10, "Cannot store solution: p = %d is not valid\n", p);
		return;
	}

	// Shift all vectors by one and store the new vector at position zero
	tmp = s[p][MAX-1];

	for(int i = MAX-1; i > 0; i--)
	{
		s[p][i] = s[p][i-1];
	}

	s[p][0] = tmp;
	spinor_field_copy_f(s[p][0], in);
	num[p]++;
}

void mre_guess(int p, spinor_field *out, spinor_operator D, spinor_field *pf)
{
	double c;
	int max;

	if(!initialized)
	{
		init();
	}

	if(p >= REPS)
	{
		lprintf("MRE", 10, "Cannot guess solution: p = %d is not valid\n", p);
		return;
	}

	spinor_field_zero_f(out);
	max = (num[p] > MAX) ? MAX : num[p];
	gram_schmidt(p, max);

	for(int i = 0; i < max; i++)
	{
		D.dbl(Dv, v[i]);

		for(int j = 0; j < max; j++)
		{
			A[j][i] = spinor_field_prod_f(Dv, v[j]);
		}

		b[i] = spinor_field_prod_f(v[i], pf);
	}

	lu_solve(max);

	for(int i = 0; i < max; i++)
	{
		c = coefficient(i+1, max);
//		spinor_field_mul_add_assign_f(out, c, s[p][i]);
		spinor_field_mulc_add_assign_f(out, x[i], v[i]);
	}
}
