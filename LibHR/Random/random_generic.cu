#include "libhr_core.h"
#include "inverters.h"
#include "random.h"
#include "utils.h"

void random_field(spinor_field *s1) {
    gaussian_spinor_field(s1);
}

void random_field(spinor_field_flt *s1) {
    gaussian_spinor_field_flt(s1);
}

void random_field(scalar_field *s1) {
    gaussian_scalar_field(s1);
}

void random_field(suNg_field *s1) {
    random_u(s1);
}

void random_field(suNf_field *s1) {
    random_u_f(s1);
}

void random_field(suNg_field_flt *s1) {
    random_suNg_field_flt_cpu(s1);
}

void random_field(suNf_field_flt *s1) {
    random_suNf_field_flt_cpu(s1);
}

void random_field(suNg_scalar_field *s1) {
    random_suNg_scalar_field_cpu(s1);
}

void random_field(suNg_av_field *s1) {
    random_suNg_av_field_cpu(s1);
}

void random_field(gtransf *s1) {
    random_gtransf_cpu(s1);
}

void random_field(ldl_field *s1) {
    random_ldl_field_cpu(s1);
}

void random_field(clover_term *s1) {
    random_clover_term_cpu(s1);
}

void random_field(clover_force *s1) {
    random_clover_force_cpu(s1);
}

void random_field(staple_field *s1) {
    random_staple_field_cpu(s1);
}