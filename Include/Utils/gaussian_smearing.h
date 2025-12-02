#ifndef GAUSS_SMEARING
#define GAUSS_SMEARING

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

void gaussian_smearing(spinor_field *restrict out, spinor_field *restrict in, suNf_field *gauge_f, double alpha);


#ifdef __cplusplus
}
#endif
#endif