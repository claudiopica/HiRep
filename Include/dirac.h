/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef DIRAC_H
#define DIRAC_H

#include "suN_types.h"
#include "utils.h"


void Dphi_(spinor_field *out, spinor_field *in);
void Dphi(double m0, spinor_field *out, spinor_field *in);
void g5Dphi(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_sq(double m0, spinor_field *out, spinor_field *in);

void Dphi_flt_(spinor_field_flt *out, spinor_field_flt *in);
void Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

unsigned long int getMVM();
unsigned long int getMVM_flt();

// Dirac operators with clover term
void Cphi(double, spinor_field*, spinor_field*);
void g5Cphi(double, spinor_field*, spinor_field*);
void g5Cphi_sq(double, spinor_field*, spinor_field*);
void Cphi_eopre(double, spinor_field*, spinor_field*);
void g5Cphi_eopre(double, spinor_field*, spinor_field*);
void g5Cphi_eopre_sq(double, spinor_field*, spinor_field*);
void Cphi_diag(double, spinor_field*, spinor_field*);
void Cphi_diag_inv(double, spinor_field*, spinor_field*);

/* Even/Odd preconditioned matrix */
void Dphi_eopre(double m0, spinor_field *out, spinor_field *in);
void Dphi_oepre(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_eopre(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_eopre_sq(double m0, spinor_field *out, spinor_field *in);

void Dphi_eopre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphi_oepre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_eopre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_eopre_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

//Twisted mass operators
void Qhat_eopre(double m0, double mu, spinor_field* out, spinor_field *in);
void Qhat_eopre_sq(double m0, double mu, spinor_field *out, spinor_field *in);


/* Dirac operators used in the Update */

void set_dirac_mass(double mass); //this is the mass used in the following operators
double get_dirac_mass();
void set_twisted_mass(double mass); //this is the twisted mass used in the twisted mass operators (Q)
void H(spinor_field *out, spinor_field *in);
void H_flt(spinor_field_flt *out, spinor_field_flt *in);
void H2(spinor_field *out, spinor_field *in);
void D(spinor_field *out, spinor_field *in);
void D_flt(spinor_field_flt *out, spinor_field_flt *in);

void Qtm_p(spinor_field *out, spinor_field *in);
void Qtm_m(spinor_field *out, spinor_field *in);
void QpQm_tm(spinor_field *out, spinor_field *in);

void Qtm_p_alt(spinor_field *out, spinor_field *in);
void Qtm_m_alt(spinor_field *out, spinor_field *in);
void QpQm_tm_alt(spinor_field *out, spinor_field *in);

//Wrapper to invert the twisted mass operator Qtm_p
void tm_invert(spinor_field* out, spinor_field *in, mshift_par* mpar);
void tm_invert_alt(spinor_field* out, spinor_field *in, mshift_par* mpar);




/* Dirac operators with a four fermion interaction */
void Dphi_eopre_4f(double m0, spinor_field *out, spinor_field *in, double shift);
void Dphi_eopre_4f_dagger(double m0, spinor_field *out, spinor_field *in, double shift);
void Dphieopre_4f_sq(double m0, spinor_field *out, spinor_field *in, double shift);

void Dphi_eopre_4f_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphi_eopre_4f_dagger_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphieopre_4f_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

/* Dirac operators used in update */
void set_ff_dirac_mass(double mass); //this is the mass used in the following operators
void set_ff_dirac_shift(double mass); //The shift added to four fermion Hasenbush-Dirac operators (Dff, Dff_dagger and Df_sq)
void Dff(spinor_field *out, spinor_field *in);
void Dff_dagger(spinor_field *out, spinor_field *in);
void Dff_sq(spinor_field *out, spinor_field *in);


#endif
	
