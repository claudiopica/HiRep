#ifndef MEASURE_RENORMALIZATION_H
#define MEASURE_RENORMALIZATION_H

#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

void measure_renormalization(spinor_field *psi_in, spinor_field *psi_out, int nm, int pt_in, int px_in, int py_in, int pz_in, int pt_out, int px_out, int py_out, int pz_out);
void print_renormalization(int conf, int nm, double *mass, char *label, int pt_in, int px_in, int py_in, int pz_in, int pt_out, int px_out, int py_out, int pz_out);

#ifdef __cplusplus
	}
#endif
#endif //MEASURE_RENORMALIZATION_H
