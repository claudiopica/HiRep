#ifndef BASIS_LINEAR_ALGEBRA_H
#define BASIS_LINEAR_ALGEBRA_H

// TODO: Merge everything into testing utilities.

// Testing utils
int check_diff_norm(double, double);
int check_diff_norm_zero(double);


// COPY 
void copy_gfield_cpu(suNg_field*, suNg_field*);
void copy_scalar_field_cpu(suNg_scalar_field*, suNg_scalar_field*);
void copy_gfield_flt_cpu(suNg_field_flt*, suNf_field_flt*);
void copy_gfield_f_cpu(suNf_field*, suNf_field*);
void copy_gfield_f_flt_cpu(suNf_field_flt*, suNf_field_flt*);
void copy_avfield_cpu(suNg_av_field*, suNg_av_field*);
void copy_sfield_cpu(scalar_field*, scalar_field*);
void copy_clover_ldl_cpu(ldl_field*, ldl_field*);
void copy_clover_term_cpu(suNfc_field*, suNfc_field*);
void copy_clover_force_cpu(suNf_field*, suNf_field*);

// SUB ASSIGN
void sub_assign_gfield_cpu(suNg_field*, suNg_field*);
void sub_assign_gfield_flt_cpu(suNg_field_flt*, suNg_field_flt*);
void sub_assign_gfield_f_cpu(suNf_field*, suNf_field*);
void sub_assign_gfield_f_flt_cpu(suNf_field_flt*, suNf_field_flt*);
void sub_assign_scalar_field_cpu(suNg_scalar_field*, suNg_scalar_field*);
void sub_assign_avfield_cpu(suNg_av_field*, suNg_av_field*);
void sub_assign_gtransf(suNg_field*, suNg_field*);
void sub_assign_clover_ldl(ldl_field*, ldl_field*);
void sub_assign_clover_term(suNfc_field*, suNfc_field*);
void sub_assign_clover_force(suNf_field*, suNf_field*);

// SQNORM
double sqnorm_gfield_cpu(suNg_field*);
double sqnorm_gfield_f_cpu(suNf_field*);
double sqnorm_gfield_flt_cpu(suNg_field_flt*);
double sqnorm_gfield_f_flt_cpu(suNf_field_flt*);
double sqnorm_scalar_field_cpu(suNg_scalar_field*);
double sqnorm_avfield_cpu(suNg_av_field*);
double sqnorm_gtransf_cpu(suNg_field*);
double sqnorm_clover_ldl_cpu(ldl_field*);
double sqnorm_clover_term_cpu(suNfc_field*);
double sqnorm_clover_force_cpu(suNf_field*);

// SET ZERO
double zero_gfield_cpu(suNg_field*);


#endif