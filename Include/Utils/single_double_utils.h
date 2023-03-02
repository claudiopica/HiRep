#ifndef SINGLE_DOUBLE_UTILS_H
#define SINGLE_DOUBLE_UTILS_H

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

//single_double_utils.c
//void assign_u2ud(void); //TODO: why are these commented out?
//void assign_ud2u(void); //TODO: why are these commented out?
void assign_ud2u_f(void);
void assign_s2sd(spinor_field *out, spinor_field_flt *in);
void assign_sd2s(spinor_field_flt *out, spinor_field *in);

#ifdef __cplusplus
}
#endif
#endif //SINGLE_DOUBLE_UTILS_H
