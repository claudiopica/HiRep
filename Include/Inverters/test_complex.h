/**
 * @file test_complex.h
 * @brief Tests for complex numbers
 */

#ifndef TEST_COMPLEX_H
#define TEST_COMPLEX_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WITH_GPU
#ifndef WITH_MPI
// These tests only exist to check
// the integration of the definition
// of the custom c++ complex types
// with the c native complex numbers
// test without MPI.

// CPU consistency checks
int test_gpu_complex();
int test_overload_plus_rhs_double();
int test_overload_plus_lhs_double();
int test_overload_prod_rhs_double();
int test_overload_prod_lhs_double();
int test_overload_div_rhs_double();
int test_overload_div_lhs_double();
int test_overload_plus_hr_complex();
int test_overload_minus_hr_complex();
int test_overload_prod_hr_complex();
int test_overload_div_hr_complex();
int test_negate();
int test_cast_double();
int test_overload_plus_rhs_integer();
int test_overload_plus_lhs_integer();
int test_overload_div_rhs_integer();
int test_overload_div_lhs_integer();
int test_I_add();
int test_I_prod();

// Test CPU against GPU
int check_plus();
int check_minus();
int check_mul();
int check_div();

#endif
#endif
#ifdef __cplusplus
}
#endif
#endif
