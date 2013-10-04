/*******************************************************************************
*
* File gamma_spinor.h
*
* Macros for gamma x spinor
*
* 2013 Rudy Arthur, Ari Hietanen
*
*******************************************************************************/

#ifndef GAMMA_SPINOR_H
#define GAMMA_SPINOR_H

#include "suN_types.h"
#include "suN.h"

#define _spinor_plus_f(out,in)\
  (out).c[0] = (in).c[0];			\
       (out).c[1] = (in).c[1];			\
       (out).c[2] = (in).c[2];			\
       (out).c[3] = (in).c[3]

#define _spinor_i_plus_f(out,in) \
  _vector_i_plus_f((out).c[0],(in).c[0]); \
  _vector_i_plus_f((out).c[1],(in).c[1]); \
  _vector_i_plus_f((out).c[2],(in).c[2]); \
  _vector_i_plus_f((out).c[3],(in).c[3])

#define _spinor_i_minus_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[0]); \
  _vector_i_minus_f((out).c[1],(in).c[1]); \
  _vector_i_minus_f((out).c[2],(in).c[2]); \
  _vector_i_minus_f((out).c[3],(in).c[3])


#define _spinor_g0_f(out,in) \
  _vector_minus_f((out).c[0],(in).c[2]); \
  _vector_minus_f((out).c[1],(in).c[3]); \
  _vector_minus_f((out).c[2],(in).c[0]); \
  _vector_minus_f((out).c[3],(in).c[1])

#define _spinor_g1_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[3]); \
  _vector_i_minus_f((out).c[1],(in).c[2]); \
  _vector_i_plus_f((out).c[2],(in).c[1]); \
  _vector_i_plus_f((out).c[3],(in).c[0])

#define _spinor_g2_f(out,in) \
  _vector_minus_f((out).c[0],(in).c[3]); \
  (out).c[1] = (in).c[2]; \
  (out).c[2] = (in).c[1]; \
  _vector_minus_f((out).c[3],(in).c[0])

#define _spinor_g3_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[2]); \
  _vector_i_plus_f((out).c[1],(in).c[3]); \
  _vector_i_plus_f((out).c[2],(in).c[0]); \
  _vector_i_minus_f((out).c[3],(in).c[1])

#define _spinor_g0g5_f(out,in) \
  (out).c[0] = (in).c[2];			\
  (out).c[1] = (in).c[3];			\
  _vector_minus_f((out).c[2],(in).c[0]); \
  _vector_minus_f((out).c[3],(in).c[1])


#define _spinor_g5g0_f(out,in) \
  _vector_minus_f((out).c[0],(in).c[2]); \
  _vector_minus_f((out).c[1],(in).c[3]); \
  (out).c[2] = (in).c[0]; \
  (out).c[3] = (in).c[1]

#define _spinor_g5g1_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[3]); \
  _vector_i_minus_f((out).c[1],(in).c[2]); \
  _vector_i_minus_f((out).c[2],(in).c[1]); \
  _vector_i_minus_f((out).c[3],(in).c[0])

#define _spinor_g5g2_f(out,in) \
  _vector_minus_f((out).c[0],(in).c[3]); \
  (out).c[1] = (in).c[2]; \
  _vector_minus_f((out).c[2],(in).c[1]);\
  (out).c[3] = (in).c[0] 

#define _spinor_g5g3_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[2]); \
  _vector_i_plus_f((out).c[1],(in).c[3]); \
  _vector_i_minus_f((out).c[2],(in).c[0]); \
  _vector_i_plus_f((out).c[3],(in).c[1])


#define _spinor_g0g1_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[1]); \
  _vector_i_minus_f((out).c[1],(in).c[0]); \
  _vector_i_plus_f((out).c[2],(in).c[3]); \
  _vector_i_plus_f((out).c[3],(in).c[2])

#define _spinor_g0g2_f(out,in) \
  _vector_minus_f((out).c[0],(in).c[1]); \
  (out).c[1] = (in).c[0]; \
  (out).c[2] = (in).c[3]; \
  _vector_minus_f((out).c[3],(in).c[2])

#define _spinor_g0g3_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[0]); \
  _vector_i_plus_f((out).c[1],(in).c[1]); \
  _vector_i_plus_f((out).c[2],(in).c[2]); \
  _vector_i_minus_f((out).c[3],(in).c[3])


#define _spinor_g5g0g1_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[1]); \
  _vector_i_minus_f((out).c[1],(in).c[0]); \
  _vector_i_minus_f((out).c[2],(in).c[3]); \
  _vector_i_minus_f((out).c[3],(in).c[2])

#define _spinor_g5g0g2_f(out,in) \
  _vector_minus_f((out).c[0],(in).c[1]); \
  (out).c[1] = (in).c[0]; \
  _vector_minus_f((out).c[2],(in).c[3]); \
  (out).c[3] = (in).c[2]

#define _spinor_g5g0g3_f(out,in) \
  _vector_i_minus_f((out).c[0],(in).c[0]); \
  _vector_i_plus_f((out).c[1],(in).c[1]); \
  _vector_i_minus_f((out).c[2],(in).c[2]); \
  _vector_i_plus_f((out).c[3],(in).c[3])



#endif
