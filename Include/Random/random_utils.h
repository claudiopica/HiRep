/// Headerfile for:
/// - gauss.c
/// - random_suNg.c
/// - random_su2.c
/// - random_fields.c
/// - ran_utils.c
/// - random_momenta.c
/// - random_spinor_field.c

#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include "suN_types.h"
#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

// gauss.c
void gauss(double r[], int n);
void gauss_flt(float r[], int n);

// random_suNg.c
void random_suNg(suNg *u);
void random_suNf(suNf *u);
void gaussian_suNg_vector(suNg_vector *v);
//void random_suNg_unit_vector(suNg_vector *v); //TODO: not defined in libhr

// random_su2.c
void random_su2(double rho, double s[]);

// random_fields.c
void random_u(suNg_field *gf);
void random_u_f(suNf_field *);
void unit_u(suNg_field *gf);
void random_s(suNg_scalar_field *);
void zero_s(suNg_scalar_field *sf);

// ran_utils.c
void ranz2(double r[], int n);
void generate_random_point(int *pr);

//random_momenta.c
void gaussian_momenta(suNg_av_field *momenta);
void gaussian_scalar_momenta(suNg_scalar_field *momenta);

//random_spinor_field.c
void gaussian_spinor_field(spinor_field *s);
void gaussian_spinor_field_flt(spinor_field_flt *s);
void z2_spinor_field(spinor_field *s);

#ifdef __cplusplus
}
#endif
#endif //RAN_UTILS_H
