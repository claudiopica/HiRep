// Header file for:
// - measure_meson_ff.c
// - four_fermion_meas.c

#ifndef MEASURE_MESONS_FF_H
#define MEASURE_MESONS_FF_H

#include "spinor_field.h"
#include "meson_observables.h"

#ifdef __cplusplus
extern "C" {
#endif

// measure_meson_ff.c
/* Disconnected part of triplet correlators */
/* appears with four fermion interactions */
extern meson_observable *triplet_discon_correlators;

void init_triplet_discon_correlators(void);
void free_triplet_discon_observables(void);
void measure_mesons_disconnected(meson_observable *mo, spinor_field *psi0, spinor_field *eta);

// four_fermion_meas.c
void ff_observables(void);

#ifdef __cplusplus
}
#endif
#endif //MEASURE_MESONS_FF_H
