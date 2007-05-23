#ifndef FUNDAMENTAL_FERMIONS_H
#define FUNDAMENTAL_FERMIONS_H

#include "suN.h"

#define INVSQRT3 .57735026918962576451

#define _group_represent(ru,u) _suNg_matrix_copy((ru),(u))

/*
 * the resulting m is HERMITEAN TRACELESS
 * m = h_a \lambda_a ; \lambda = Gell-Mann matrices
 */
/*
#define _algebra_represent(m,h) \
        (m).c1_1.re = (h).c3 + INVSQRT3*(h).c8 ;\
        (m).c1_1.im = 0. ; \
        (m).c1_2.re = (h).c1 ; \
        (m).c1_2.im = -(h).c2 ; \
        (m).c1_3.re = (h).c4 ; \
        (m).c1_3.im = -(h).c5 ; \
        (m).c2_1.re = (h).c1 ; \
        (m).c2_1.im = (h).c2 ; \
        (m).c2_2.re = INVSQRT3*(h).c8-(h).c3 ; \
        (m).c2_2.im = 0. ; \
        (m).c2_3.re = (h).c6 ; \
        (m).c2_3.im = -(h).c7 ; \
        (m).c3_1.re = (h).c4 ; \
        (m).c3_1.im = (h).c5 ; \
        (m).c3_2.re = (h).c6 ; \
        (m).c3_2.im = (h).c7 ; \
        (m).c3_3.re = -2.*INVSQRT3*(h).c8 ; \
        (m).c3_3.im = 0.
*/

/*
 * Martin basis
 */
#define _algebra_represent(m,h) \
        (m).c1_1.re = 0. ;\
        (m).c1_1.im = (h).c1 + (h).c2 ; \
        (m).c1_2.re = (h).c3 ; \
        (m).c1_2.im = (h).c4 ; \
        (m).c1_3.re = (h).c5 ; \
        (m).c1_3.im = (h).c6 ; \
        (m).c2_1.re = -(h).c3 ; \
        (m).c2_1.im = (h).c4 ; \
        (m).c2_2.re = 0. ; \
        (m).c2_2.im = (h).c2-2.*(h).c1 ; \
        (m).c2_3.re = (h).c7 ; \
        (m).c2_3.im = (h).c8 ; \
        (m).c3_1.re = -(h).c5 ; \
        (m).c3_1.im = (h).c6 ; \
        (m).c3_2.re = -(h).c7 ; \
        (m).c3_2.im = (h).c8 ; \
        (m).c3_3.re = 0. ; \
        (m).c3_3.im = (h).c1-2.*(h).c2

#define _algebra_project(h,m) \
        (h).c1 = 0.33333333333333333333*((m).c1_1.im-(m).c2_2.im); \
        (h).c2 = 0.33333333333333333333*((m).c1_1.im-(m).c3_3.im); \
        (h).c3 = (m).c1_2.re; \
        (h).c4 = (m).c1_2.im; \
        (h).c5 = (m).c1_3.re; \
        (h).c6 = (m).c1_3.im; \
        (h).c7 = (m).c2_3.re; \
        (h).c8 = (m).c2_3.im

#endif
