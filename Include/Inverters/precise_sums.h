/***************************************************************************\
* Copyright (c) 2022, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef PRECISE_SUMS_H
#define PRECISE_SUMS_H

#include <math.h>

///compute the compensating term with Neumaier trick
#define _comp_term(s,a,b)  (fabs(a)>fabs(b))?(((a)-(s))+(b)):(((b)-(s))+(a))

/// Exact addition of two single-length floating point numbers
/// The macro produces a double-length number (s,c) that satisfies
/// s+c = a+b exactly, with s = round(a + b)
#define _2Sum_acc(s, c, a) \
do {                           \
   double _a=a;                \
   double _t=s+_a;             \
   c += _comp_term(_t, s, _a); \
   s = _t;                     \
} while (0)                          

///second order 2Sum algo
#define _2Sum_acc2(s, c, cc, a) \
do {                                     \
   double _a=a;                          \
   double _t=s+_a;                       \
   double _ct = _comp_term(_t, s, _a);   \
   s = _t;                               \
   _t = c + _ct;                         \
   double _cct = _comp_term(_t, c, _ct); \
   c = _t;                               \
   cc += _cct;                           \
} while (0)

/// this is the versin in OpenQCD
#define _2Sum_acc_oqcd(s, c, a) \
do {                          \
   double _a=s+a;             \
   double _qp=_a-a;           \
   double _up=_a-_qp;         \
   double _b=(s-_qp)+(a-_up); \
   double _c=c+_b;            \
   double _d=_a+_c;           \
   s=_d;                      \
   c=_c-(_d-_a);              \
} while (0)

#define _2Sum_acc_naive(s, c, a) (s)+=(a) 


#endif