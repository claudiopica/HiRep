/***************************************************************************
* Copyright (c) 2016, Martin Hansen                                        *
* All rights reserved.                                                     *
***************************************************************************/

#include "global.h"
#include "memory.h"
#include "suN.h"
#include "logger.h"
#include "communications.h"
#include "clover_tools.h"
#include "clover_exp.h"
#include "utils.h"
#include <math.h>
#include <string.h>

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
static double sigma;
static double csw_value;

#define clover(id, mu, ndx) \
	_4FIELD_AT(cl_term, id, mu)->c[ndx]

#define clover_re(id, mu, ndx) \
	creal(_4FIELD_AT(cl_term, id, mu)->c[ndx])

#define clover_im(id, mu, ndx) \
	cimag(_4FIELD_AT(cl_term, id, mu)->c[ndx])

#define clover_force(id, mu, ndx) \
	_6FIELD_AT(cl_force, id, mu)->c[ndx]

#define clover_force_re(id, mu, ndx) \
	_6FIELD_AT(cl_force, id, mu)->c[ndx].re

#define clover_force_im(id, mu, ndx) \
	_6FIELD_AT(cl_force, id, mu)->c[ndx].im

#define cplx(x, i, j) \
	(x[(i) * ((i) + 1) / 2 + (j)])

#define re(x, i, j) \
	creal(x[(i) * ((i) + 1) / 2 + (j)])

#define im(x, i, j) \
	cimag(x[(i) * ((i) + 1) / 2 + (j)])

double get_csw()
{
	return csw_value;
}

void set_csw(double * csw)
{
	csw_value = * csw;
	lprintf("CLOVER", 10, "Coefficient: reset to csw = %1.6f\n", csw_value);
}

static void clover_loop(int id, int mu, int nu, suNf *u)
{
	int o1, o2, o3;
	suNf s1, s2, s3;

	// Leaf 1
	o1 = iup(id, mu);
	o2 = iup(id, nu);
	_suNf_times_suNf(s1, *pu_gauge_f(id, mu), *pu_gauge_f(o1, nu));
	_suNf_times_suNf(s2, *pu_gauge_f(id, nu), *pu_gauge_f(o2, mu));
	_suNf_times_suNf_dagger(*u, s1, s2);

	// Leaf 2
	o1 = idn(id, mu);
	o2 = iup(o1, nu);
	_suNf_times_suNf(s1, *pu_gauge_f(o1, nu), *pu_gauge_f(o2, mu));
	_suNf_times_suNf_dagger(s2, *pu_gauge_f(id, nu), s1);
	_suNf_times_suNf(s3, s2, *pu_gauge_f(o1, mu));
	_suNf_add_assign(*u, s3);

	// Leaf 3
	o1 = idn(id, mu);
	o2 = idn(id, nu);
	o3 = idn(o1, nu);
	_suNf_times_suNf(s1, *pu_gauge_f(o3, nu), *pu_gauge_f(o1, mu));
	_suNf_dagger_times_suNf(s2, s1, *pu_gauge_f(o3, mu));
	_suNf_times_suNf(s3, s2, *pu_gauge_f(o2, nu));
	_suNf_add_assign(*u, s3);

	// Leaf 4
	o1 = idn(id, nu);
	o2 = iup(o1, mu);
	_suNf_dagger_times_suNf(s1, *pu_gauge_f(o1, nu), *pu_gauge_f(o1, mu));
	_suNf_times_suNf_dagger(s2, *pu_gauge_f(o2, nu), *pu_gauge_f(id, mu));
	_suNf_times_suNf(s3, s1, s2);
	_suNf_add_assign(*u, s3);
}

static void ldl(int N, double complex *A)
{
	for (int i = 0; i < N; i++)
	{
		for (int k = 0; k < i; k++)
		{
			cplx(A, i, i) -= (_complex_prod(cplx(A, i, k), cplx(A, i, k))) * re(A, k, k);
		}
		for (int j = i + 1; j < N; j++)
		{
			for (int k = 0; k < i; k++)
			{
				cplx(A, j, i) -= (_complex_prod(cplx(A, i, k), cplx(A, j, k))) * re(A, k, k);
			}
			cplx(A, j, i) /= re(A, i, i);
		}
	}
}

static void _compute_ldl_decomp(int id)
{
	int m, n, ij, ji;
	double complex *A, *B;

	// Setup pointers
	A = _FIELD_AT(cl_ldl, id)->up;
	B = _FIELD_AT(cl_ldl, id)->dn;

	// Construct matrices
	for (int i = 0; i < NF; i++)
	{
		m = i + NF;
		for (int j = 0; j < NF; j++)
		{
			n = j + NF;
			ij = i * NF + j;
			ji = j * NF + i;

			cplx(A, m, j) = conj(clover(id, 1, ji));
			cplx(B, m, j) = conj(clover(id, 3, ji));

			if (i >= j)
			{
				cplx(A, i, j) = clover(id, 0, ij);
				cplx(A, m, n) = -clover(id, 0, ij);
				cplx(B, i, j) = clover(id, 2, ij);
				cplx(B, m, n) = -clover(id, 2, ij);

				if (i == j)
				{
					cplx(A, i, j) += sigma;
					cplx(A, m, n) += sigma;
					cplx(B, i, j) += sigma;
					cplx(B, m, n) += sigma;
				}
			}
		}
	}

	// LDL factorization
	ldl(2 * NF, A);
	ldl(2 * NF, B);
}

static void _compute_clover_term(int id)
{
	suNf tmp[6];
	double csw;

	double complex atmp;
	double complex btmp;
	double complex ctmp;
	double complex dtmp;

	csw = csw_value;
	csw = -csw / 16.0;

	clover_loop(id, 0, 1, &tmp[0]);
	clover_loop(id, 0, 2, &tmp[1]);
	clover_loop(id, 0, 3, &tmp[2]);
	clover_loop(id, 1, 2, &tmp[3]);
	clover_loop(id, 1, 3, &tmp[4]);
	clover_loop(id, 2, 3, &tmp[5]);

	for (int i = 0; i < NF; i++)
	{
		for (int j = 0; j < NF; j++)
		{
			int ij = i * NF + j;
			int ji = j * NF + i;

#ifdef REPR_IS_REAL
			atmp = I * (tmp[2].c[ji] - tmp[2].c[ij]);
			btmp = I * (tmp[3].c[ij] - tmp[3].c[ji]);
			ctmp = tmp[1].c[ji] - tmp[1].c[ij] + I * (tmp[0].c[ji] - tmp[0].c[ij]);
			dtmp = tmp[4].c[ij] - tmp[4].c[ji] + I * (tmp[5].c[ji] - tmp[5].c[ij]);

#else
			atmp = I * (conj(tmp[2].c[ji]) - tmp[2].c[ij]);
			btmp = I * (conj(tmp[3].c[ij]) - tmp[3].c[ji]);
			ctmp = I * (conj(tmp[0].c[ji]) - tmp[0].c[ij]) - tmp[1].c[ij] + conj(tmp[1].c[ji]);
			dtmp = tmp[4].c[ij] - conj(tmp[4].c[ji]) + I * (conj(tmp[5].c[ji]) - tmp[5].c[ij]);
#endif
			clover(id, 0, ij) = csw * (atmp - conj(btmp));
			clover(id, 1, ij) = csw * (ctmp - dtmp);
			clover(id, 2, ij) = -csw * (atmp + conj(btmp));
			clover(id, 3, ij) = -csw * (ctmp + dtmp);
		}
	}
}

static void _compute_clover_force(int id, double coeff)
{
	double complex A[2 * NF][2 * NF];
	double complex B[2 * NF][2 * NF];
	memset(A, 0, sizeof(A));
	memset(B, 0, sizeof(B));

	double complex a11;
	double complex a12;
	double complex a21;
	double complex a22;
	double complex a33;
	double complex a34;
	double complex a43;
	double complex a44;

	double complex *U = _FIELD_AT(cl_ldl, id)->up;
	double complex *L = _FIELD_AT(cl_ldl, id)->dn;

	// Calculate inverse from LDL
	for (int n = 0; n < 2 * NF; n++)
	{
		A[n][n] = coeff;
		B[n][n] = coeff;

		for (int i = n; i < 2 * NF; i++)
		{
			for (int k = 0; k < i; k++)
			{
				A[i][n] -= cplx(U, i, k) * A[k][n];
				B[i][n] -= cplx(L, i, k) * B[k][n];
			}
		}
		for (int i = 2 * NF - 1; i >= n; i--)
		{
			A[i][n] /= re(U, i, i);
			B[i][n] /= re(L, i, i);
			for (int k = i + 1; k < 2 * NF; k++)
			{
				A[i][n] -= conj(cplx(U, k, i)) * A[k][n];
				B[i][n] -= conj(cplx(L, k, i)) * B[k][n];
			}
		}
	}

	// Construct force matrices
	for (int i = 0; i < NF; i++)
	{
		for (int j = 0; j < NF; j++)
		{
			int ij = i * NF + j;

			a21 = A[i + NF][j];
			a12 = conj(A[j + NF][i]);
			a43 = B[i + NF][j];
			a34 = conj(B[j + NF][i]);

			if (i < j)
			{
				a11 = conj(A[j][i]);
				a22 = conj(A[j + NF][i + NF]);
				a33 = conj(B[j][i]);
				a44 = conj(B[j + NF][i + NF]);
			}
			else
			{
				a11 = A[i][j];
				a22 = A[i + NF][j + NF];
				a33 = B[i][j];
				a44 = B[i + NF][j + NF];
			}

#ifdef REPR_IS_REAL
			clover_force(id, 0, ij) += cimag(a12) + cimag(a21) - cimag(a34) - cimag(a43); // X_01
			clover_force(id, 1, ij) += creal(a12) - creal(a21) + creal(a43) - creal(a34); // X_02
			clover_force(id, 2, ij) += cimag(a22) - cimag(a11) + cimag(a44) - cimag(a33); // X_12
			clover_force(id, 3, ij) += cimag(a11) - cimag(a22) + cimag(a44) - cimag(a33); // X_03
			clover_force(id, 4, ij) += creal(a12) - creal(a21) + creal(a34) - creal(a43); // X_13
			clover_force(id, 5, ij) -= cimag(a12) + cimag(a21) + cimag(a34) + cimag(a43); // X_23
#else
			clover_force(id, 0, ij) -= I * (a12 + a21 - a34 - a43); // X_01
			clover_force(id, 1, ij) += a12 - a21 + a43 - a34;		// X_02
			clover_force(id, 2, ij) -= I * (a22 - a11 + a44 - a33); // X_12
			clover_force(id, 3, ij) -= I * (a11 - a22 + a44 - a33); // X_03
			clover_force(id, 4, ij) += a12 - a21 + a34 - a43;		// X_13
			clover_force(id, 5, ij) += I * (a12 + a21 + a34 + a43); // X_23

#endif
		}
	}
}

void compute_force_logdet(double mass, double coeff)
{
	// Update LDL decomposition
	compute_ldl_decomp(4.0 + mass);

	// Loop odd sites
	_MASTER_FOR(&glat_odd, id)
	{
		_compute_clover_force(id, coeff);
	}
}

void clover_la_logdet(double nf, double mass, scalar_field *la)
{
	// Update LDL decomposition
	compute_ldl_decomp(4.0 + mass);

	// Add local action
	_MASTER_FOR(&glat_odd, id)
	{
		double prod = 1;
		double complex *A = _FIELD_AT(cl_ldl, id)->up;
		double complex *B = _FIELD_AT(cl_ldl, id)->dn;

		for (int n = 0; n < 2 * NF; n++)
		{
			prod *= re(A, n, n) * re(B, n, n);
		}

		*_FIELD_AT(la, id) -= nf * log(prod);
	}
}

void compute_clover_term()
{
	sigma = 0xF00F;
	_MASTER_FOR(&glattice, id)
	{
		_compute_clover_term(id);
	}
	apply_BCs_on_clover_term(cl_term);
}

void compute_ldl_decomp(double sigma0)
{
	if (sigma == sigma0)
	{
		return;
	}
	else
	{
		sigma = sigma0;
	}

	_MASTER_FOR(&glattice, id)
	{
		_compute_ldl_decomp(id);
	}
}

void clover_init(double csw)
{
	cl_term = alloc_clover_term(&glattice);
	cl_ldl = alloc_clover_ldl(&glattice);
	cl_force = alloc_clover_force(&glattice);

	sigma = 0xF00F;
	csw_value = csw;
	lprintf("CLOVER", 10, "Initial Coefficient: csw = %1.6f\n", csw_value);

#if defined(WITH_EXPCLOVER)
	init_clover_exp();
#endif
}

#endif //#ifdef WITH_CLOVER
