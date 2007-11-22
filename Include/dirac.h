#ifndef DIRAC_H
#define DIRAC_H

#include "suN_types.h"

typedef enum {
   EO=0,
   OE=1
} block_selector;

void Dphi_(block_selector B, spinor_field *out, spinor_field *in);
void Dphi(double m0, spinor_field *out, spinor_field *in);
void g5Dphi(double m0, spinor_field *out, spinor_field *in);

void Dphi_flt_(block_selector B, spinor_field_flt *out, spinor_field_flt *in);
void Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

unsigned long int getMVM();
unsigned long int getMVM_flt();

/* p = out ; q = in */

void Dphi_old(double m0, suNf_spinor *p, suNf_spinor *q);
void g5Dphi_old(double m0, suNf_spinor *p, suNf_spinor *q);
void Dphi_flt_old(double m0, suNf_spinor_flt *p, suNf_spinor_flt *q);
void g5Dphi_flt_old(double m0, suNf_spinor_flt *p, suNf_spinor_flt *q);


#endif
	
