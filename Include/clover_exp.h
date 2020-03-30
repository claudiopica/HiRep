#ifndef CLOVER_EXP_H
#define CLOVER_EXP_H
#include "suN_types.h"
#include "hr_complex.h"
#include "spinor_field.h"
#include "suN.h"

void _su2Nfc_times_su2Nfc_herm(suNfc *C, suNfc *B, suNfc *A);
void _su2Nfc_times_su2Nf(suNfc *C, suNfc *B, suNfc *A);
void _su2Nfc_times_su2Nfc_assign(suNfc *C, suNfc *B, suNfc *A);
void _su2Nfc_times_su2Nfc_assign_herm(suNfc *C, suNfc *B, suNfc *A);
void _su2Nfc_times_su2Nfc_trace(double complex *trace, suNfc *B, suNfc *A);
void _su2Nfc_times_su2Nfc_trace_herm_sq(double complex *trace, suNfc *B);
void _su2Nfc_unit(suNfc *A);
void _su2Nfc_trace(double complex *p, suNfc *A);
void clover_exp(suNfc *Aplus, suNfc *expAplus);
void evaluate_sw_order(double * mass);

void clover_exp_taylor(suNfc *Aplus, suNfc *expAplus);
void doublehorner(double *C, suNfc *A);
void init_clover_exp();
int get_NNexp();
void factorialCoef(double *C);
#endif
