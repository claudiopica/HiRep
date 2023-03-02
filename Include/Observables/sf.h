// header for files:
// - sf_action.c
// - sf_pcac.c

#ifndef SF_H
#define SF_H

#include "spinor_field.h"
#include "Utils/data_storage.h"

#ifdef __cplusplus
extern "C" {
#endif

// sf_action.c
double SF_action(double beta);

// sf_pcac.c
int SF_quark_propagator(spinor_field *in, double mass, spinor_field *out, double acc);
data_storage_array *SF_PCAC_wall_corr(double mass, double acc, storage_switch swh);

#ifdef __cplusplus
}
#endif
#endif //SF_H
