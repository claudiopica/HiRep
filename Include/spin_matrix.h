/*******************************************************************************
*
* File spin_matrix.h
*
* Type definitions and macros for spin_matrix
*
* 2013 Rudy Arthur, Ari Hietanen
*
*
*  |S_00 S_01| |1| = |S_00| = c[0]
*  |S_10 S_11| |0|   |S_10|
*
*  etc.
*  The notation I don't like _spinor_g5_f(r) = g5 r
*  becomes the more sensible _g5_spinmatrix(s) = g5 s and _spinmatrix_g5(s) = s g5
*
*******************************************************************************/

#ifndef SPIN_MATRIX_H
#define SPIN_MATRIX_H
#include "suN_types.h"
#include "suN.h"
#include "gamma_spinor.h"

typedef struct
{
   suNf_spinor c[4];
} suNf_spin_matrix;


/*  r=0  (r spinmatrix) */
#define _spinmatrix_zero(r) \
  _spinor_zero_f((r).c[0]); \
  _spinor_zero_f((r).c[1]); \
  _spinor_zero_f((r).c[2]); \
  _spinor_zero_f((r).c[3])

//r spinmatrix, s spinor, i index
#define _spinmatrix_assign_row(r, s, i) \
  (r).c[i] = (s)

//r spinmatrix, s spinor, i index
/*#define _spinmatrix_assign_col(r, s, i) \
  (r).c[0].c[i] = (s).c[0]; \
  (r).c[1].c[i] = (s).c[1]; \
  (r).c[2].c[i] = (s).c[2]; \
  (r).c[3].c[i] = (s).c[3]*/


//r = s + t 
#define _spinmatrix_add(r,s,t) \
  _spinor_add_f((r).c[0],(s).c[0],(t).c[0]);\
  _spinor_add_f((r).c[1],(s).c[1],(t).c[1]);\
  _spinor_add_f((r).c[2],(s).c[2],(t).c[2]);\
  _spinor_add_f((r).c[3],(s).c[3],(t).c[3])

//r = s - t 
#define _spinmatrix_sub(r,s,t) \
  _spinor_sub_f((r).c[0],(s).c[0],(t).c[0]);\
  _spinor_sub_f((r).c[1],(s).c[1],(t).c[1]);\
  _spinor_sub_f((r).c[2],(s).c[2],(t).c[2]);\
  _spinor_sub_f((r).c[3],(s).c[3],(t).c[3])
  

//s =  r  . gamma_5
#define _spinmatrix_g5(s, r) \
  _spinor_g5_f((s).c[0],(r).c[0]);		\
  _spinor_g5_f((s).c[1],(r).c[1]); \
  _spinor_g5_f((s).c[2],(r).c[2]); \
  _spinor_g5_f((s).c[3],(r).c[3])

//s = gamma5 . r 
#define _g5_spinmatrix(s, r) \
  (s).c[0] = (r).c[0]; \
  (s).c[1] = (r).c[1]; \
  _spinor_minus_f((s).c[2],(r).c[2]); \
  _spinor_minus_f((s).c[3],(r).c[3])

//s = r . gamma_0 
#define _spinmatrix_g0(s, r) \
  _spinor_g0_f((s).c[0],(r).c[0]);		\
  _spinor_g0_f((s).c[1],(r).c[1]); \
  _spinor_g0_f((s).c[2],(r).c[2]); \
  _spinor_g0_f((s).c[3],(r).c[3])

//s = gamma_0 . r
#define _g0_spinmatrix(s, r)	      \
  _spinor_minus_f((s).c[0],(r).c[2]); \
  _spinor_minus_f((s).c[1],(r).c[3]); \
  _spinor_minus_f((s).c[2],(r).c[0]); \
  _spinor_minus_f((s).c[3],(r).c[1])


//s = r . gamma_1 
#define _spinmatrix_g1(s, r) \
  _spinor_g1_f((s).c[0],(r).c[0]);		\
  _spinor_g1_f((s).c[1],(r).c[1]); \
  _spinor_g1_f((s).c[2],(r).c[2]); \
  _spinor_g1_f((s).c[3],(r).c[3])

//s = gamma_1 . r
#define _g1_spinmatrix(s, r) \
  _spinor_i_minus_f((s).c[0],(r).c[3]); \
  _spinor_i_minus_f((s).c[1],(r).c[2]); \
  _spinor_i_plus_f((s).c[2],(r).c[1]); \
  _spinor_i_plus_f((s).c[3],(r).c[0])


//s = r . gamma_2 
#define _spinmatrix_g2(s, r) \
  _spinor_g2_f((s).c[0],(r).c[0]);		\
  _spinor_g2_f((s).c[1],(r).c[1]); \
  _spinor_g2_f((s).c[2],(r).c[2]); \
  _spinor_g2_f((s).c[3],(r).c[3])

//s = gamma_2 . r
#define _g2_spinmatrix(s, r) \
  _spinor_minus_f((s).c[0],(r).c[3]); \
  _spinor_plus_f((s).c[1],(r).c[2]); \
  _spinor_plus_f((s).c[2],(r).c[1]); \
  _spinor_minus_f((s).c[3],(r).c[0])

//s = r . gamma_3
#define _spinmatrix_g3(s, r) \
  _spinor_g3_f((s).c[0],(r).c[0]);		\
  _spinor_g3_f((s).c[1],(r).c[1]); \
  _spinor_g3_f((s).c[2],(r).c[2]); \
  _spinor_g3_f((s).c[3],(r).c[3])

//s = gamma_3 . r
#define _g3_spinmatrix(s, r)	      \
  _spinor_i_minus_f((s).c[0],(r).c[2]); \
  _spinor_i_plus_f((s).c[1],(r).c[3]); \
  _spinor_i_plus_f((s).c[2],(r).c[0]); \
  _spinor_i_minus_f((s).c[3],(r).c[1])


//s = gamma_0 gamma_5 . r 
#define _spinmatrix_g0g5(s, r) \
  _spinor_g0g5_f((s).c[0],(r).c[0]); \
  _spinor_g0g5_f((s).c[1],(r).c[1]); \
  _spinor_g0g5_f((s).c[2],(r).c[2]); \
  _spinor_g0g5_f((s).c[3],(r).c[3])


//s = gamma_0 . gamma_5 . r
#define _g0g5_spinmatrix(s, r) \
  _spinor_plus_f((s).c[0],(r).c[2]); \
  _spinor_plus_f((s).c[1],(r).c[3]); \
  _spinor_minus_f((s).c[2],(r).c[0]);\
  _spinor_minus_f((s).c[3],(r).c[1]);

//s = gamma_5 gamma_0 . r 
#define _spinmatrix_g5g0(s, r) \
  _spinor_g5g0_f((s).c[0],(r).c[0]); \
  _spinor_g5g0_f((s).c[1],(r).c[1]); \
  _spinor_g5g0_f((s).c[2],(r).c[2]); \
  _spinor_g5g0_f((s).c[3],(r).c[3])

//s = gamma_5 . gamma_0 . r
#define _g5g0_spinmatrix(s, r) \
  _spinor_minus_f((s).c[0],(r).c[2]); \
  _spinor_minus_f((s).c[1],(r).c[3]); \
  (s).c[2] = (r).c[0]; \
  (s).c[3] = (r).c[1]


//s = r . gamma_5 . gamma_1 
#define _spinmatrix_g5g1(s, r) \
  _spinor_g5g1_f((s).c[0],(r).c[0]); \
  _spinor_g5g1_f((s).c[1],(r).c[1]); \
  _spinor_g5g1_f((s).c[2],(r).c[2]); \
  _spinor_g5g1_f((s).c[3],(r).c[3])

//s = gamma_5 . gamma_1 . r
#define _g5g1_spinmatrix(s, r) \
  _spinor_i_minus_f((s).c[0],(r).c[3]); \
  _spinor_i_minus_f((s).c[1],(r).c[2]); \
  _spinor_i_minus_f((s).c[2],(r).c[1]); \
  _spinor_i_minus_f((s).c[3],(r).c[0])

//s = r . gamma_5 .gamma_2 
#define _spinmatrix_g5g2(s, r) \
  _spinor_g5g2_f((s).c[0],(r).c[0]); \
  _spinor_g5g2_f((s).c[1],(r).c[1]); \
  _spinor_g5g2_f((s).c[2],(r).c[2]); \
  _spinor_g5g2_f((s).c[3],(r).c[3])

//s = gamma_5.gamma_2 . r
#define _g5g2_spinmatrix(s, r) \
  _spinor_minus_f((s).c[0],(r).c[3]); \
  (s).c[1] = (r).c[2]; \
  _spinor_minus_f((s).c[2],(r).c[1]); \
  (s).c[3] = (r).c[0]

//s = r . gamma_5 gamma_3
#define _spinmatrix_g5g3(s, r) \
  _spinor_g5g3_f((s).c[0],(r).c[0]); \
  _spinor_g5g3_f((s).c[1],(r).c[1]); \
  _spinor_g5g3_f((s).c[2],(r).c[2]); \
  _spinor_g5g3_f((s).c[3],(r).c[3])

//s = gamma_5 . gamma_3 . r
#define _g5g3_spinmatrix(s, r) \
  _spinor_i_minus_f((s).c[0],(r).c[2]); \
  _spinor_i_plus_f((s).c[1],(r).c[3]); \
  _spinor_i_minus_f((s).c[2],(r).c[0]);	   \
  _spinor_i_plus_f((s).c[3],(r).c[1])



//s = r . gamma_0 . gamma_1 
#define _spinmatrix_g0g1(s, r) \
  _spinor_g0g1_f((s).c[0],(r).c[0]); \
  _spinor_g0g1_f((s).c[1],(r).c[1]); \
  _spinor_g0g1_f((s).c[2],(r).c[2]); \
  _spinor_g0g1_f((s).c[3],(r).c[3])

//s =  gamma_0 . gamma_1 . r
#define _g0g1_spinmatrix(s, r) \
  _spinor_i_minus_f((s).c[0],(r).c[1]); \
  _spinor_i_minus_f((s).c[1],(r).c[0]); \
  _spinor_i_plus_f((s).c[2],(r).c[3]); \
  _spinor_i_plus_f((s).c[3],(r).c[2])

//s = r . gamma_0 . gamma_2 
#define _spinmatrix_g0g2(s, r) \
  _spinor_g0g2_f((s).c[0],(r).c[0]); \
  _spinor_g0g2_f((s).c[1],(r).c[1]); \
  _spinor_g0g2_f((s).c[2],(r).c[2]); \
  _spinor_g0g2_f((s).c[3],(r).c[3])

//s = gamma0 . gamma_2 . r
#define _g0g2_spinmatrix(s, r) \
  _spinor_minus_f((s).c[0],(r).c[1]); \
  _spinor_plus_f((s).c[1],(r).c[0]); \
  _spinor_plus_f((s).c[2],(r).c[3]); \
  _spinor_minus_f((s).c[3],(r).c[2])

//s = r . gamma0 .gamma_3
#define _spinmatrix_g0g3(s, r) \
  _spinor_g0g3_f((s).c[0],(r).c[0]); \
  _spinor_g0g3_f((s).c[1],(r).c[1]); \
  _spinor_g0g3_f((s).c[2],(r).c[2]); \
  _spinor_g0g3_f((s).c[3],(r).c[3])

//s = gamma_5 . gamma_3 . r
#define _g0g3_spinmatrix(s, r) \
  _spinor_i_minus_f((s).c[0],(r).c[0]); \
  _spinor_i_plus_f((s).c[1],(r).c[1]); \
  _spinor_i_plus_f((s).c[2],(r).c[2]);	   \
  _spinor_i_minus_f((s).c[3],(r).c[3])


//s = r . gamma_5 . gamma_0 . gamma_1 
#define _spinmatrix_g5g0g1(s, r) \
  _spinor_g5g0g1_f((s).c[0],(r).c[0]); \
  _spinor_g5g0g1_f((s).c[1],(r).c[1]); \
  _spinor_g5g0g1_f((s).c[2],(r).c[2]); \
  _spinor_g5g0g1_f((s).c[3],(r).c[3])

//s = gamma_5 . gamma_0 . gamma_1 . r
#define _g5g0g1_spinmatrix(s, r) \
  _spinor_i_minus_f((s).c[0],(r).c[1]); \
  _spinor_i_minus_f((s).c[1],(r).c[0]); \
  _spinor_i_minus_f((s).c[2],(r).c[3]); \
  _spinor_i_minus_f((s).c[3],(r).c[2])

//s = r . gamma_5 . gamma_0 . gamma_2 
#define _spinmatrix_g5g0g2(s, r) \
  _spinor_g5g0g2_f((s).c[0],(r).c[0]); \
  _spinor_g5g0g2_f((s).c[1],(r).c[1]); \
  _spinor_g5g0g2_f((s).c[2],(r).c[2]); \
  _spinor_g5g0g2_f((s).c[3],(r).c[3])

//s = gamma_5. gamma0 . gamma_2 . r
#define _g5g0g2_spinmatrix(s, r) \
  _spinor_minus_f((s).c[0],(r).c[1]); \
  (s).c[1] = (r).c[0]; \
  _spinor_minus_f((s).c[2],(r).c[3]); \
  (s).c[3] = (r).c[2]

//s = r . gamma_5 . gamma0 .gamma_3
#define _spinmatrix_g5g0g3(s, r) \
  _spinor_g5g0g3_f((s).c[0],(r).c[0]); \
  _spinor_g5g0g3_f((s).c[1],(r).c[1]); \
  _spinor_g5g0g3_f((s).c[2],(r).c[2]); \
  _spinor_g5g0g3_f((s).c[3],(r).c[3])

//s = gamma_5 . gamma_3 . r
#define _g5g0g3_spinmatrix(s, r) \
  _spinor_i_minus_f((s).c[0],(r).c[0]); \
  _spinor_i_plus_f((s).c[1],(r).c[1]); \
  _spinor_i_minus_f((s).c[2],(r).c[2]);	   \
  _spinor_i_plus_f((s).c[3],(r).c[3])


//r spinmatrix, s spinmatrix, k result; Tr [ r^dag . s]
#define _spinmatrix_mul_trace(k, r, s) \
   do { \
      (k).re=0.;(k).im=0.; \
      _spinor_prod_assign_f((k),(r).c[0],(s).c[0]); \
      _spinor_prod_assign_f((k),(r).c[1],(s).c[1]); \
      _spinor_prod_assign_f((k),(r).c[2],(s).c[2]); \
      _spinor_prod_assign_f((k),(r).c[3],(s).c[3]); \
   } while(0) 

//r spinmatrix, s spinmatrix, k += Tr [ r^dag . s]
#define _spinmatrix_mul_trace_assign(k, r, s) \
   do { \
       _spinor_prod_assign_f((k),(r).c[0],(s).c[0]); \
       _spinor_prod_assign_f((k),(r).c[1],(s).c[1]); \
       _spinor_prod_assign_f((k),(r).c[2],(s).c[2]); \
       _spinor_prod_assign_f((k),(r).c[3],(s).c[3]); \
   } while(0) 

#define _spinmatrix_mul_trace_re(k, r, s) \
    do { \
      double _tmpVAR;							\
     (k)=0;								\
     _spinor_prod_re_f((_tmpVAR),(r).c[0],(s).c[0]); (k)+=_tmpVAR;		\
     _spinor_prod_re_f((_tmpVAR),(r).c[1],(s).c[1]); (k)+=_tmpVAR;	\
     _spinor_prod_re_f((_tmpVAR),(r).c[2],(s).c[2]); (k)+=_tmpVAR;	\
     _spinor_prod_re_f((_tmpVAR),(r).c[3],(s).c[3]); (k)+=_tmpVAR;	\
    } while(0) 
	    
	    
#endif
