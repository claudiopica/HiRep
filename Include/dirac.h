#ifndef DIRAC_H
#define DIRAC_H

#include "suN_types.h"

typedef enum {
   EO=0,
   OE=1
} block_selector;

void Dphi_(block_selector B, suNf_spinor *out, suNf_spinor *in);
void Dphi(double m0, suNf_spinor *out, suNf_spinor *in);
void g5Dphi(double m0, suNf_spinor *out, suNf_spinor *in);
void Dphi_eopre(double m0, suNf_spinor *out, suNf_spinor *in);
void g5Dphi_eopre(double m0, suNf_spinor *out, suNf_spinor *in);

void Dphi_flt_(block_selector B, suNf_spinor_flt *out, suNf_spinor_flt *in);
void Dphi_flt(double m0, suNf_spinor_flt *out, suNf_spinor_flt *in);
void g5Dphi_flt(double m0, suNf_spinor_flt *out, suNf_spinor_flt *in);
void Dphi_eopre_flt(double m0, suNf_spinor_flt *out, suNf_spinor_flt *in);
void g5Dphi_eopre_flt(double m0, suNf_spinor_flt *out, suNf_spinor_flt *in);

unsigned long int getMVM();
unsigned long int getMVM_flt();

/* p = out ; q = in */

void Dphi_old(double m0, suNf_spinor *p, suNf_spinor *q);
void g5Dphi_old(double m0, suNf_spinor *p, suNf_spinor *q);
void Dphi_flt_old(double m0, suNf_spinor_flt *p, suNf_spinor_flt *q);
void g5Dphi_flt_old(double m0, suNf_spinor_flt *p, suNf_spinor_flt *q);


#endif
	
