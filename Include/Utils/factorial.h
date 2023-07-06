/***************************************************************************
 * Copyright (c) 2023, Sofie Martins                                        *
 * All rights reserved.                                                     *
 ***************************************************************************/

#ifndef FACTORIAL_H
#define FACTORIAL_H

#include "libhr_core.h"

#define MAX_FACTORIAL 100
// To match with clover exp
// But for the exp in the field update
// we need only 30

#ifdef __cplusplus
extern "C" {
#endif

void init_factorial(void);
void finalize_factorial(void);
double inverse_factorial(int);
deviceonly double inverse_factorial_gpu(int);

#ifdef __cplusplus
}
#endif
#endif