/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File single_double_utils_gpu.c
*
* Functions for conversion from single to double precision and viceversa wiht gpu
*
*******************************************************************************/

#include <stdlib.h>
#include "utils.h"
#include "suN.h"
#include "error.h"
#include "global.h"
#include "spinor_field.h"
#include "memory.h"


void assign_s2sd(spinor_field *out, spinor_field_flt *in) {
  spinor_field_copy_from_gpu_f_flt(in);
  assign_s2sd_cpu(out,in);
  spinor_field_copy_to_gpu_f(out);
}

void assign_sd2s(spinor_field_flt *out, spinor_field *in) {
  spinor_field_copy_from_gpu_f(in);
  assign_sd2s_cpu(out,in);
  spinor_field_copy_to_gpu_f_flt(out);
}
