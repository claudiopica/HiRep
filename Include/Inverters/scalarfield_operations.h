/**
 * @file scalarfield_operations.h
 * @brief Linear algebra operations for scalar fields
 */

#ifndef SCALARFIELD_OPERATIONS_H
#define SCALARFIELD_OPERATIONS_H

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

void flip_scalar_field(scalar_field *s);
void set_scalar_field(scalar_field *s, double c);
void gaussian_scalar_field(scalar_field *s);
void spinor_scalarfield_mult_add_assign(spinor_field *out, scalar_field *sigma, double rho, spinor_field *in);
void spinor_scalarfield_ig5_mult_add_assign(spinor_field *out, scalar_field *pi, spinor_field *in);
void spinor_scalarfield_mig5_mult_add_assign(spinor_field *out, scalar_field *pi, spinor_field *in);
void spinor_sigma_pi_rho_minus_div_assign(spinor_field *out, scalar_field *sigma, scalar_field *pi, double rho,
                                          spinor_field *in);
void spinor_sigma_pi_rho_div_assign(spinor_field *out, scalar_field *sigma, scalar_field *pi, double rho, spinor_field *in);
void spinor_sigma_pi_dagger_rho_div_assign(spinor_field *out, scalar_field *sigma, scalar_field *pi, double rho,
                                           spinor_field *in);
void spinor_sigma_pi_dagger_rho_assign(spinor_field *out, scalar_field *sigma, scalar_field *pi, double rho, spinor_field *in);
void spinor_sigma_assign(spinor_field *out, scalar_field *sigma, spinor_field *in);
void spinor_scalarfield_mult_add_assign_flt(spinor_field_flt *out, scalar_field *sigma, double rho, spinor_field_flt *in);
void spinor_scalarfield_ig5_mult_add_assign_flt(spinor_field_flt *out, scalar_field *pi, spinor_field_flt *in);
void spinor_scalarfield_mig5_mult_add_assign_flt(spinor_field_flt *out, scalar_field *pi, spinor_field_flt *in);
void spinor_sigma_pi_rho_div_assign_flt(spinor_field_flt *out, scalar_field *sigma, scalar_field *pi, double rho,
                                        spinor_field_flt *in);
void spinor_sigma_pi_dagger_rho_div_assign_flt(spinor_field_flt *out, scalar_field *sigma, scalar_field *pi, double rho,
                                               spinor_field_flt *in);

#ifdef __cplusplus
}
#endif
#endif