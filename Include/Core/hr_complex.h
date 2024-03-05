/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file hr_complex.h
 * @brief Type definitions and macros for complex numbers
 */

#ifndef HR_COMPLEX_H
#define HR_COMPLEX_H

#ifdef __cplusplus
#include "gpu_complex.hpp"
#else
#include "cpu_complex.h"
#endif
#endif

#define _complex_re(a) creal(a)
#define _complex_im(a) cimag(a)
#define _complex_0(a) (a) = 0
#define _complex_1(a) (a) = 1
#define _complex_add_1(a) (a) += 1
#define _complex_i(a) (a) = I
#define _complex_star(a, b) (a) = conj(b)
#define _complex_star_minus(a, b) (a) = -conj(b)
#define _complex_star_assign(a) (a) = conj(a)
#define _complex_mul(a, b, c) (a) = (b) * (c)
#define _complex_mulr(a, r, b) (a) = (r) * (b)
#define _complex_add(a, b, c) (a) = (b) + (c)
#define _complex_sub(a, b, c) (a) = (b) - (c)
#define _complex_add_star(a, b, c) (a) = (b) + conj(c)
#define _complex_sub_star(a, b, c) (a) = (b)-conj(c)
#define _complex_div(a, b, c) (a) = (b) / (c)
#define _complex_inv(a, b) (a) = 1 / (b)
#define _complex_prod(a, b) conj(a) * b
#define _complex_prod_re(a, b) creal(conj(a) * (b))
#define _complex_prod_m1_re(a, b) creal(conj(1 - a) * (1 - b))
#define _complex_prod_im(a, b) cimag(conj(a) * b)
#define _complex_prod_assign(c, a, b) (c) += conj(a) * (b)
#define _complex_mul_star_star_assign(c, a, b) (c) += conj((a) * (b))
#define _complex_minus(a, b) (a) = -(b)
#define _complex_i_minus(a, b) (a) = -I * (b)
#define _complex_i_plus(a, b) (a) = I * (b)
#define _complex_i_add(a, b, c) (a) = (b) + I * (c)
#define _complex_i_sub(a, b, c) (a) = (b)-I * (c)
#define _complex_add_assign(a, b) (a) += (b)
#define _complex_sub_assign(a, b) (a) -= (b)
#define _complex_add_star_assign(a, b, c) (a) += (b) + conj(c)
#define _complex_i_add_assign(a, b) (a) += I * (b)
#define _complex_i_sub_assign(a, b) (a) -= I * (b)
#define _complex_mul_assign(a, b, c) (a) += (b) * (c)
#define _complex_mulcr_assign(a, r, b, c) (a) += (r) * (b) * (c)
#define _complex_mul_star(a, b, c) (a) = (b) * conj(c)
#define _complex_mul_star_assign(a, b, c) (a) += (b) * conj(c)
#define _complex_mul_star_assign_re(a, b, c) (a) += creal((b) * conj(c))
#define _complex_mul_sub_assign(a, b, c) (a) -= (b) * (c)
#define _complex_mulr_assign(a, r, b) (a) += (r) * (b)
#define _complex_rlc(a, r1, c1, r2, c2) (a) = (r1) * (c1) + (r2) * (c2)
#define _complex_rlc_assign(a, r1, c1, r2, c2) (a) += (r1) * (c1) + (r2) * (c2)
#define _complex_clc(a, z1, c1, z2, c2) (a) = (z1) * (c1) + (z2) * (c2)
#define _complex_clc_assign(a, z1, c1, z2, c2) (a) += (z1) * (c1) + (z2) * (c2)
