#ifndef CALC_PROP_FF_H
#define CALC_PROP_FF_H

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Functions that include four fermion interactions */
void init_propagator_ff_eo(int nm, double *m, double acc);
void free_propagator_ff_eo();

void calc_propagator_ff(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator_ff_eo(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator_ff_oe(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator_ff_hopping_eo(spinor_field *psi, spinor_field *eta, int hopping, int ndilute);
void calc_propagator_ff_hopping_oe(spinor_field *psi, spinor_field *eta, int hopping, int ndilute);
void calc_propagator_ff_hopping_series(spinor_field *psi, spinor_field *eta, int hopping, int ndilute);

#ifdef __cplusplus
}
#endif
#endif //CALC_PROP_FF_H
