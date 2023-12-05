#ifndef RANDOM_GENERIC_GPU_HPP
#define RANDOM_GENERIC_GPU_HPP

#ifdef __cplusplus

#include "inverters.h"
#include "libhr_core.h"

void random_field(spinor_field *s1);
void random_field(spinor_field_flt *s1);
void random_field(scalar_field *s1);
void random_field(suNg_field *s1);
void random_field(suNf_field *s1);
void random_field(suNg_field_flt *s1);
void random_field(suNf_field_flt *s1);
void random_field(suNg_scalar_field *s1);
void random_field(suNg_av_field *s1);
void random_field(gtransf *s1);
void random_field(ldl_field *s1);
void random_field(clover_term *s1);
void random_field(clover_force *s1);
void random_field(staple_field *s1);

#endif
#endif