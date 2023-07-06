#ifndef CLOVER_EXP_H
#define CLOVER_EXP_H

#ifdef WITH_EXPCLOVER

#include "libhr_core.h"

#ifdef __cplusplus
extern "C" {
#endif

visible void _su2Nfc_times_su2Nfc_herm(suNfc *C, suNfc *B, suNfc *A);
visible void _su2Nfc_times_su2Nf(suNfc *C, suNfc *B, suNfc *A);
visible void _su2Nfc_times_su2Nfc_assign(suNfc *C, suNfc *B, suNfc *A);
visible void _su2Nfc_times_su2Nfc_assign_herm(suNfc *C, suNfc *B, suNfc *A);
visible void _su2Nfc_times_su2Nfc_trace(hr_complex *trace, suNfc *B, suNfc *A);
visible void _su2Nfc_times_su2Nfc_trace_herm_sq(hr_complex *trace, suNfc *B);
visible void _su2Nfc_unit(suNfc *A);
visible void _su2Nfc_trace(hr_complex *p, suNfc *A);
void factorialCoef(double *C, int NN, int NNexp);

void clover_exp(suNfc *Aplus, suNfc *expAplus);
void evaluate_sw_order(double *mass);
void doublehorner(double *C, suNfc *A);

void init_clover_exp();
int get_NNexp();
int get_NN();

#ifdef WITH_GPU
deviceonly void clover_exp_gpu(suNfc *Aplus, suNfc *expAplus, int NN, int NNexp);
#endif

#ifdef __cplusplus
}
#endif
#endif
#endif
