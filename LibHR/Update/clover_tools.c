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
#include "utils.h"
#include <math.h>
#include <string.h>

static double sigma;
static double csw_value;

#ifdef WITH_CLOVER

#define clover_re(id,mu,ndx) \
	_4FIELD_AT(cl_term,id,mu)->c[ndx].re

#define clover_im(id,mu,ndx) \
	_4FIELD_AT(cl_term,id,mu)->c[ndx].im

#define clover_force(id,mu,ndx) \
	_6FIELD_AT(cl_force,id,mu)->c[ndx]

#define clover_force_re(id,mu,ndx) \
	_6FIELD_AT(cl_force,id,mu)->c[ndx].re

#define clover_force_im(id,mu,ndx) \
	_6FIELD_AT(cl_force,id,mu)->c[ndx].im

#define re(x,i,j) \
	x[(i)*((i)+1)/2+(j)].re

#define im(x,i,j) \
	x[(i)*((i)+1)/2+(j)].im

double get_csw()
{
	return csw_value;
}

static void clover_loop(int id, int mu, int nu, suNf *u)
{
	int o1, o2, o3;
	suNf s1, s2, s3;

	// Leaf 1
	o1 = iup(id,mu);
	o2 = iup(id,nu);
	_suNf_times_suNf(s1, *pu_gauge_f(id,mu), *pu_gauge_f(o1,nu));
	_suNf_times_suNf(s2, *pu_gauge_f(id,nu), *pu_gauge_f(o2,mu));
	_suNf_times_suNf_dagger(*u, s1, s2);

	// Leaf 2
	o1 = idn(id,mu);
	o2 = iup(o1,nu);
	_suNf_times_suNf(s1, *pu_gauge_f(o1,nu), *pu_gauge_f(o2,mu));
	_suNf_times_suNf_dagger(s2, *pu_gauge_f(id,nu), s1);
	_suNf_times_suNf(s3, s2, *pu_gauge_f(o1,mu));
	_suNf_add_assign(*u, s3);

	// Leaf 3
	o1 = idn(id,mu);
	o2 = idn(id,nu);
	o3 = idn(o1,nu);
	_suNf_times_suNf(s1, *pu_gauge_f(o3,nu), *pu_gauge_f(o1,mu));
	_suNf_dagger_times_suNf(s2, s1, *pu_gauge_f(o3,mu));
	_suNf_times_suNf(s3, s2, *pu_gauge_f(o2,nu));
	_suNf_add_assign(*u, s3);

	// Leaf 4
	o1 = idn(id,nu);
	o2 = iup(o1,mu);
	_suNf_dagger_times_suNf(s1, *pu_gauge_f(o1,nu), *pu_gauge_f(o1,mu));
	_suNf_times_suNf_dagger(s2, *pu_gauge_f(o2,nu), *pu_gauge_f(id,mu));
	_suNf_times_suNf(s3, s1, s2);
	_suNf_add_assign(*u, s3);
}

static void ldl(int N, complex *A)
{
	for(int i = 0; i < N; i++)
	{
		for(int k = 0; k < i; k++)
		{
			re(A,i,i) -= (re(A,i,k)*re(A,i,k) + im(A,i,k)*im(A,i,k)) * re(A,k,k);
		}
		for(int j = i+1; j < N; j++)
		{
			for(int k = 0; k < i; k++)
			{
				re(A,j,i) -= (re(A,j,k)*re(A,i,k) + im(A,j,k)*im(A,i,k)) * re(A,k,k);
				im(A,j,i) -= (im(A,j,k)*re(A,i,k) - re(A,j,k)*im(A,i,k)) * re(A,k,k);
			}
			re(A,j,i) /= re(A,i,i);
			im(A,j,i) /= re(A,i,i);
		}
	}
}

static void _compute_ldl_decomp(int id)
{
	int m, n, ij, ji;
	complex *A, *B;

	// Setup pointers
	A = _FIELD_AT(cl_ldl,id)->up;
	B = _FIELD_AT(cl_ldl,id)->dn;

	// Construct matrices
	for(int i = 0; i < NF; i++)
	{
		m = i+NF;
		for(int j = 0; j < NF; j++)
		{
			n = j+NF;
			ij = i*NF+j;
			ji = j*NF+i;

			re(A,m,j) =  clover_re(id,1,ji);
			im(A,m,j) = -clover_im(id,1,ji);
			re(B,m,j) =  clover_re(id,3,ji);
			im(B,m,j) = -clover_im(id,3,ji);

			if(i >= j)
			{
				re(A,i,j) =  clover_re(id,0,ij);
				im(A,i,j) =  clover_im(id,0,ij);
				re(A,m,n) = -clover_re(id,0,ij);
				im(A,m,n) = -clover_im(id,0,ij);
				re(B,i,j) =  clover_re(id,2,ij);
				im(B,i,j) =  clover_im(id,2,ij);
				re(B,m,n) = -clover_re(id,2,ij);
				im(B,m,n) = -clover_im(id,2,ij);

				if(i == j)
				{
					re(A,i,j) += sigma;
					re(A,m,n) += sigma;
					re(B,i,j) += sigma;
					re(B,m,n) += sigma;
				}
			}
		}
	}

	// LDL factorization
	ldl(2*NF, A);
	ldl(2*NF, B);
}

static void _compute_clover_term(int id)
{
	suNf tmp[6];
	double csw;
	double atmp_re, atmp_im;
	double btmp_re, btmp_im;
	double ctmp_re, ctmp_im;
	double dtmp_re, dtmp_im;

	csw = csw_value;
	csw = -csw/16.0;

	clover_loop(id, 0, 1, &tmp[0]);
	clover_loop(id, 0, 2, &tmp[1]);
	clover_loop(id, 0, 3, &tmp[2]);
	clover_loop(id, 1, 2, &tmp[3]);
	clover_loop(id, 1, 3, &tmp[4]);
	clover_loop(id, 2, 3, &tmp[5]);

	for(int i = 0; i < NF; i++)
	{
		for(int j = 0; j < NF; j++)
		{
			int ij = i*NF+j;
			int ji = j*NF+i;

#ifdef REPR_IS_REAL
			atmp_re = 0;
			atmp_im = tmp[2].c[ji] - tmp[2].c[ij];
			btmp_re = 0;
			btmp_im = tmp[3].c[ij] - tmp[3].c[ji];
			ctmp_re = tmp[1].c[ji] - tmp[1].c[ij];
			ctmp_im = tmp[0].c[ji] - tmp[0].c[ij];
			dtmp_re = tmp[4].c[ij] - tmp[4].c[ji];
			dtmp_im = tmp[5].c[ji] - tmp[5].c[ij];
#else
			atmp_re = tmp[2].c[ij].im + tmp[2].c[ji].im;
			atmp_im = tmp[2].c[ji].re - tmp[2].c[ij].re;
			btmp_re = tmp[3].c[ij].im + tmp[3].c[ji].im;
			btmp_im = tmp[3].c[ij].re - tmp[3].c[ji].re;
			ctmp_re = tmp[0].c[ij].im + tmp[0].c[ji].im - tmp[1].c[ij].re + tmp[1].c[ji].re;
			ctmp_im = tmp[0].c[ji].re - tmp[0].c[ij].re - tmp[1].c[ij].im - tmp[1].c[ji].im;
			dtmp_re = tmp[4].c[ij].re - tmp[4].c[ji].re + tmp[5].c[ij].im + tmp[5].c[ji].im;
			dtmp_im = tmp[4].c[ij].im + tmp[4].c[ji].im - tmp[5].c[ij].re + tmp[5].c[ji].re;
#endif

			clover_re(id,0,ij) =  csw * (atmp_re - btmp_re);
			clover_im(id,0,ij) =  csw * (atmp_im + btmp_im);
			clover_re(id,1,ij) =  csw * (ctmp_re - dtmp_re);
			clover_im(id,1,ij) =  csw * (ctmp_im - dtmp_im);
			clover_re(id,2,ij) = -csw * (atmp_re + btmp_re);
			clover_im(id,2,ij) =  csw * (btmp_im - atmp_im);
			clover_re(id,3,ij) = -csw * (ctmp_re + dtmp_re);
			clover_im(id,3,ij) = -csw * (ctmp_im + dtmp_im);
		}
	}
}

static void _compute_clover_force(int id, double coeff)
{
	complex A[2*NF][2*NF];
	complex B[2*NF][2*NF];
	memset(A, 0, sizeof(A));
	memset(B, 0, sizeof(B));

	double a11_re, a11_im;
	double a12_re, a12_im;
	double a21_re, a21_im;
	double a22_re, a22_im;
	double a33_re, a33_im;
	double a34_re, a34_im;
	double a43_re, a43_im;
	double a44_re, a44_im;

	complex *U = _FIELD_AT(cl_ldl,id)->up;
	complex *L = _FIELD_AT(cl_ldl,id)->dn;

	// Calculate inverse from LDL
	for(int n = 0; n < 2*NF; n++)
	{
		A[n][n].re = coeff;
		B[n][n].re = coeff;

		for(int i = n; i < 2*NF; i++)
		{
			for(int k = 0; k < i; k++)
			{
				A[i][n].re -= re(U,i,k)*A[k][n].re - im(U,i,k)*A[k][n].im;
				A[i][n].im -= re(U,i,k)*A[k][n].im + im(U,i,k)*A[k][n].re;
				B[i][n].re -= re(L,i,k)*B[k][n].re - im(L,i,k)*B[k][n].im;
				B[i][n].im -= re(L,i,k)*B[k][n].im + im(L,i,k)*B[k][n].re;
			}
		}

		for(int i = 2*NF-1; i >= n; i--)
		{
			A[i][n].re /= re(U,i,i);
			A[i][n].im /= re(U,i,i);
			B[i][n].re /= re(L,i,i);
			B[i][n].im /= re(L,i,i);

			for(int k = i+1; k < 2*NF; k++)
			{
				A[i][n].re -= re(U,k,i)*A[k][n].re + im(U,k,i)*A[k][n].im;
				A[i][n].im -= re(U,k,i)*A[k][n].im - im(U,k,i)*A[k][n].re;
				B[i][n].re -= re(L,k,i)*B[k][n].re + im(L,k,i)*B[k][n].im;
				B[i][n].im -= re(L,k,i)*B[k][n].im - im(L,k,i)*B[k][n].re;
			}
		}
	}

	// Construct force matrices
	for(int i = 0; i < NF; i++)
	{
		for(int j = 0; j < NF; j++)
		{
			int ij = i*NF+j;
			a21_re =  A[i+NF][j].re;
			a21_im =  A[i+NF][j].im;
			a12_re =  A[j+NF][i].re;
			a12_im = -A[j+NF][i].im;
			a43_re =  B[i+NF][j].re;
			a43_im =  B[i+NF][j].im;
			a34_re =  B[j+NF][i].re;
			a34_im = -B[j+NF][i].im;

			if(i < j)
			{
				a11_re =  A[j][i].re;
				a11_im = -A[j][i].im;
				a22_re =  A[j+NF][i+NF].re;
				a22_im = -A[j+NF][i+NF].im;
				a33_re =  B[j][i].re;
				a33_im = -B[j][i].im;
				a44_re =  B[j+NF][i+NF].re;
				a44_im = -B[j+NF][i+NF].im;
			}
			else
			{
				a11_re = A[i][j].re;
				a11_im = A[i][j].im;
				a22_re = A[i+NF][j+NF].re;
				a22_im = A[i+NF][j+NF].im;
				a33_re = B[i][j].re;
				a33_im = B[i][j].im;
				a44_re = B[i+NF][j+NF].re;
				a44_im = B[i+NF][j+NF].im;
			}

#ifdef REPR_IS_REAL
			clover_force(id,0,ij) +=  a12_im + a21_im - a34_im - a43_im; // X_01
			clover_force(id,1,ij) +=  a12_re - a21_re + a43_re - a34_re; // X_02
			clover_force(id,2,ij) +=  a22_im - a11_im + a44_im - a33_im; // X_12
			clover_force(id,3,ij) +=  a11_im - a22_im + a44_im - a33_im; // X_03
			clover_force(id,4,ij) +=  a12_re - a21_re + a34_re - a43_re; // X_13
			clover_force(id,5,ij) += -a12_im - a21_im - a34_im - a43_im; // X_23
#else
			clover_force_re(id,0,ij) +=  a12_im + a21_im - a34_im - a43_im; // X_01
			clover_force_im(id,0,ij) += -a12_re - a21_re + a34_re + a43_re;
			clover_force_re(id,1,ij) +=  a12_re - a21_re + a43_re - a34_re; // X_02
			clover_force_im(id,1,ij) +=  a12_im - a21_im + a43_im - a34_im;
			clover_force_re(id,2,ij) +=  a22_im - a11_im + a44_im - a33_im; // X_12
			clover_force_im(id,2,ij) +=  a11_re - a22_re + a33_re - a44_re;
			clover_force_re(id,3,ij) +=  a11_im - a22_im + a44_im - a33_im; // X_03
			clover_force_im(id,3,ij) +=  a22_re - a11_re + a33_re - a44_re;
			clover_force_re(id,4,ij) +=  a12_re - a21_re + a34_re - a43_re; // X_13
			clover_force_im(id,4,ij) +=  a12_im - a21_im + a34_im - a43_im;
			clover_force_re(id,5,ij) += -a12_im - a21_im - a34_im - a43_im; // X_23
			clover_force_im(id,5,ij) +=  a12_re + a21_re + a34_re + a43_re;
#endif
		}
	}
}

void compute_force_logdet(double mass, double coeff)
{
	// Update LDL decomposition
	compute_ldl_decomp(4.0+mass);

	// Loop odd sites
	_MASTER_FOR(&glat_odd,id)
	{
		_compute_clover_force(id, coeff);
	}
}

void clover_la_logdet(double nf, double mass, scalar_field *la)
{
	// Update LDL decomposition
	compute_ldl_decomp(4.0+mass);

	// Add local action
	_MASTER_FOR(&glat_odd,id)
	{
		double prod = 1;
		complex *A = _FIELD_AT(cl_ldl,id)->up;
		complex *B = _FIELD_AT(cl_ldl,id)->dn;

		for(int n = 0; n < 2*NF; n++)
		{
			prod *= re(A,n,n);
			prod *= re(B,n,n);
		}

		*_FIELD_AT(la,id) -= nf*log(prod);
	}
}

void compute_clover_term()
{
	sigma = 0xF00F;
	_MASTER_FOR(&glattice,id)
	{
		_compute_clover_term(id);
	}
	apply_BCs_on_clover_term(cl_term);
}

void compute_ldl_decomp(double sigma0)
{
	if(sigma == sigma0)
	{
		return;
	}
	else
	{
		sigma = sigma0;
	}

	_MASTER_FOR(&glattice,id)
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
	lprintf("CLOVER", 10, "Coefficient: csw = %1.6f\n", csw_value);
}

#endif //#ifdef WITH_CLOVER