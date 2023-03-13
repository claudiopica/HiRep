#ifndef SINGLE_DOUBLE_UTILS_H
#define SINGLE_DOUBLE_UTILS_H

#include "libhr_core.h"

#ifdef __cplusplus
extern "C" {
#endif

void assign_ud2u_cpu(void);
void assign_u2ud_cpu(void);
void assign_ud2u_f_cpu(void);
void assign_u2ud_f_cpu(void);
void assign_s2sd_cpu(spinor_field *out, spinor_field_flt *in);
void assign_sd2s_cpu(spinor_field_flt *out, spinor_field *in);
void add_assign_s2sd_cpu(spinor_field *out, spinor_field_flt *in);
void add_assign_sd2s_cpu(spinor_field_flt *out, spinor_field *in);

#ifdef WITH_GPU
void assign_ud2u_gpu(void);
void assign_u2ud_gpu(void);
void assign_ud2u_f_gpu(void);
void assign_u2ud_f_gpu(void);
void assign_s2sd_gpu(spinor_field *out, spinor_field_flt *in);
void assign_sd2s_gpu(spinor_field_flt *out, spinor_field *in);
void add_assign_s2sd_gpu(spinor_field *out, spinor_field_flt *in);
void add_assign_sd2s_gpu(spinor_field_flt *out, spinor_field *in);
#endif

extern void (*assign_ud2u)(void);
extern void (*assign_u2ud)(void);
extern void (*assign_ud2u_f)(void);
extern void (*assign_u2ud_f)(void);
extern void (*assign_s2sd)(spinor_field *out, spinor_field_flt *in);
extern void (*assign_sd2s)(spinor_field_flt *out, spinor_field *in);
extern void (*add_assign_s2sd)(spinor_field *out, spinor_field_flt *in);
extern void (*add_assign_sd2s)(spinor_field_flt *out, spinor_field *in);

#ifdef __cplusplus
}
#endif
#endif
