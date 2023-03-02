/**
 * @file test_complex.h
 * @brief Tests for complex numbers
 */

#ifndef TEST_COMPLEX_H
#define TEST_COMPLEX_H

#ifdef __cplusplus
extern "C" {
#endif

int test_gpu_complex();
int test_overload_plus_rhs_double();
int test_overload_plus_rhs_double_gpu();
int test_overload_plus_lhs_double();
int test_overload_prod_rhs_double();
int test_overload_prod_lhs_double();
int test_overload_div_rhs_double();
int test_overload_div_lhs_double();
int test_overload_plus_hr_complex();
int test_overload_minus_hr_complex();
int test_overload_prod_hr_complex();
int test_overload_prod_hr_complex_gpu();
int test_overload_div_hr_complex();
int test_negate();
int test_cast_double();
int test_overload_plus_rhs_integer();
int test_overload_plus_lhs_integer();
int test_overload_div_rhs_integer();
int test_overload_div_lhs_integer();
int test_I_add();
int test_I_prod();

#ifdef __cplusplus
}
#endif
#endif
