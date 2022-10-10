#ifndef BASIS_LINEAR_ALGEBRA_H
#define BASIS_LINEAR_ALGEBRA_H

// COPY 
void copy_suNg_field_cpu(suNg_field*, suNg_field*);
void copy_suNg_scalar_field_cpu(suNg_scalar_field*, suNg_scalar_field*);
void copy_suNg_field_flt_cpu(suNg_field_flt*, suNf_field_flt*);
void copy_suNf_field_f_cpu(suNf_field*, suNf_field*);
void copy_suNf_field_f_flt_cpu(suNf_field_flt*, suNf_field_flt*);
void copy_suNg_av_field_cpu(suNg_av_field*, suNg_av_field*);
void copy_scalar_field_cpu(scalar_field*, scalar_field*);
void copy_ldl_field_cpu(ldl_field*, ldl_field*);
void copy_suNfc_field_cpu(suNfc_field*, suNfc_field*);

// SUB ASSIGN
void sub_assign_suNg_field_cpu(suNg_field*, suNg_field*);
void sub_assign_suNg_scalar_field_cpu(suNg_scalar_field*, suNg_scalar_field*);
void sub_assign_suNf_field_cpu(suNf_field*, suNf_field*);


// SQNORM
double sqnorm_suNg_field_cpu(suNg_field*);
double sqnorm_suNf_field_cpu(suNf_field*);


#endif