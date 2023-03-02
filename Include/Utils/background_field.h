#ifndef BACKGROUND_FIELD_H
#define BACKGROUND_FIELD_H
#ifdef __cplusplus
extern "C" {
#endif

#include "spinor_field.h"

void apply_background_field_zdir(suNg_field *V, double Q, int n);

#ifdef __cplusplus
}
#endif
#endif //BACKGROUND_FIELD_H
