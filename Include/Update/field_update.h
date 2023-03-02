#ifndef FIELD_UPDATE_H
#define FIELD_UPDATE_H

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    suNg_field **field;
    suNg_av_field **momenta;
} field_gauge_par;

typedef struct {
    suNg_scalar_field **field;
    suNg_scalar_field **momenta;
} field_scalar_par;

//update_field.c
void update_gauge_field(double, void *);
void update_scalar_field(double, void *);

#ifdef __cplusplus
}
#endif
#endif //FIELD_UPDATE_H
