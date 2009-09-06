/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef DIRAC_H
#define DIRAC_H

#include "suN_types.h"

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

/* Even/Odd preconditioned matrix */
void Dphi_eopre(double m0, spinor_field *out, spinor_field *in);
void Dphi_oepre(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_eopre(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_eopre_sq(double m0, spinor_field *out, spinor_field *in);

void Dphi_eopre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphi_oepre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_eopre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_eopre_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

/* p = out ; q = in */

void Dphi_old(double m0, suNf_spinor *p, suNf_spinor *q);
void g5Dphi_old(double m0, suNf_spinor *p, suNf_spinor *q);
void Dphi_flt_old(double m0, suNf_spinor_flt *p, suNf_spinor_flt *q);
void g5Dphi_flt_old(double m0, suNf_spinor_flt *p, suNf_spinor_flt *q);


#endif
	
