#ifndef STOUT_SMEARING_H
#define STOUT_SMEARING_H

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

/* stout smearing */
void init_smearing(double, double);
double avr_smeared_plaquette();
void smear_gauge_field();
void smeared_gauge_force(suNg_av_field *, suNg_av_field *);

#ifdef __cplusplus
}
#endif
#endif //STOUT_SMEARING_H
