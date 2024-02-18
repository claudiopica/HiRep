/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file cpu_complex.h
 * @brief Type definitions and macros for complex numbers used in C
 */

#ifndef CPU_COMPLEX_H
#define CPU_COMPLEX_H

#include <tgmath.h>
// tgmath includes math.h and complex.h
// and defines type-generic macros for math functions
// e.g: float complex fc; creal(fc) invokes crealf(fc)

typedef double complex hr_complex;
typedef float complex hr_complex_flt;
typedef int complex hr_complex_int;

#endif
