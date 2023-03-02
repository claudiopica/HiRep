#ifndef COPY_GFIELD_H
#define COPY_GFIELD_H

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

void suNg_field_copy(suNg_field *g1, suNg_field *g2);
void suNf_field_copy(suNf_field *g1, suNf_field *g2);
void suNg_scalar_field_copy(suNg_scalar_field *g1, suNg_scalar_field *g2);

#ifdef __cplusplus
}
#endif
#endif //COPY_GFIELD_H
