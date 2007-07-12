#ifndef DIRAC_H
#define DIRAC_H

#include "suN_types.h"

typedef enum {
   EO=0,
   OE=1
} block_selector;

void Dphi_(block_selector B, suNf_spinor *out, suNf_spinor *in);
void Dphi(float m0, suNf_spinor *out, suNf_spinor *in);
void g5Dphi(float m0, suNf_spinor *out, suNf_spinor *in);

unsigned long int getMVM();

/* p = out ; q = in */
void Dphi_old(float m0, suNf_spinor *p, suNf_spinor *q);
void g5Dphi_old(float m0, suNf_spinor *p, suNf_spinor *q);
void Dphi_dble_old(double m0, suNf_spinor_dble *p, suNf_spinor_dble *q);
void g5Dphi_dble_old(double m0, suNf_spinor_dble *p, suNf_spinor_dble *q);

#endif
	
