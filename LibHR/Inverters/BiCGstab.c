/***************************************************************************\
* Copyright (c) 2014, Claudio Pica, Martin Hansen                           *
* All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "linear_algebra.h"
#include "complex.h"
#include "update.h"
#include "memory.h"
#include "logger.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
 * Performs BiCGstab inversions: out = M^-1*in
 * Returns the number of iterations spend during the inversion
 */
int BiCGstab(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out)
{
	spinor_field *s;
	spinor_field *r, *r1, *o, *Ms, *Mo, *o0;
	spinor_field *sptmp;

	complex delta, phi;
	complex alpha, beta, chi;
	complex ctmp1, ctmp2;
	double rtmp1;
	double innorm2;

	int cgiter = 0;
	int notconverged = 1;

	s  = alloc_spinor_field_f(7, in->type);
	r  = s+1;
	r1 = s+2;
	o  = s+3;
	Ms = s+4;
	Mo = s+5;
	o0 = s+6;

	// Init recursion
	_complex_1(beta);
	_complex_1(chi);
	_complex_0(alpha);

	// Initial residue
	M(Ms, out);
	cgiter++;
	spinor_field_sub_f(r, in, Ms);

	spinor_field_copy_f(s, r);
	innorm2 = spinor_field_sqnorm_f(in);

	// Choose omega so that delta and phi are not zero
	gaussian_spinor_field(o0);
	spinor_field_copy_f(o, o0);
	delta = spinor_field_prod_f(o, r);

	M(Ms, s);
	ctmp1 = spinor_field_prod_f(o0, Ms);
	_complex_div(phi, ctmp1, delta);

	// BiCGstab recursion
	while((par->max_iter==0 || cgiter<par->max_iter) && notconverged)
	{
		cgiter += 2;

		// Compute beta
		_complex_inv(beta, phi); /* b=1/phi */
		_complex_minus(beta, beta); /* b=-b */

		// Compute omega and chi
		spinor_field_mulc_f(o,beta,Ms);
		spinor_field_add_assign_f(o,r);

		M(Mo,o);

		ctmp1 = spinor_field_prod_f(Mo,o);
		rtmp1 = 1./spinor_field_sqnorm_f(Mo);
		_complex_mulr(chi, rtmp1, ctmp1);

		// Compute r1
		spinor_field_mulc_f(r1, chi, Mo);
		spinor_field_sub_f(r1, o, r1);

		// Update delta and alpha
		_complex_mul(ctmp1, delta, chi);
		delta = spinor_field_prod_f(o0, r1);
		_complex_minus(ctmp1, ctmp1);
		_complex_mul(ctmp2, beta, delta);
		_complex_div(alpha, ctmp2, ctmp1); // alpha = -beta*delta/(delta*chi)

		// Update solution
		_complex_minus(ctmp1, beta);
		spinor_field_clc_add_assign_f(out, ctmp1, s, chi, o);

		// Compute s
		_complex_mul(ctmp1, alpha, chi);
		_complex_minus(ctmp1, ctmp1);
		spinor_field_clc_f(Mo, alpha, s, ctmp1, Ms);
		spinor_field_add_f(s, r1, Mo);

		// Exchange r <-> r1
		sptmp = r;
		r = r1;
		r1 = sptmp;

		// Update phi
		M(Ms, s);
		ctmp1 = spinor_field_prod_f(o0, Ms);
		_complex_div(phi, ctmp1, delta);

		rtmp1 = spinor_field_sqnorm_f(r);
		if(rtmp1 < par->err2*innorm2)
		{
			notconverged = 0;
		}
	}

	// Test result
	M(Ms, out);
	spinor_field_sub_assign_f(Ms, in);
	innorm2 = spinor_field_sqnorm_f(Ms) / innorm2;
	if(fabs(innorm2) > par->err2)
	{
		lprintf("INVERTER", 30, "BiCGstab failed: err2 = %1.8e > %1.8e\n",innorm2, par->err2);
	}
	else
	{
		lprintf("INVERTER", 20, "BiCGstab inversion: err2 = %1.8e < %1.8e\n", innorm2, par->err2);
	}
	
	// Free memory
	free_spinor_field_f(s);

	// Print log info
	lprintf("INVERTER", 10, "BiCGstab: cgiter = %d\n", cgiter);

	// Return number of iterations
	return cgiter;
}
