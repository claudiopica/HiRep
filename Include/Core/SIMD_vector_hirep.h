#ifndef SIMD_VECTOR_HIREP_H
#define SIMD_VECTOR_HIREP_H

#ifdef SIMD_VECTOR_HIREP

#if (NF == 3) || (NG == 3)
#include "IO/logger.h"
typedef double suNg_vector_V __attribute__((vector_size(2 * NG * sizeof(double))));

#define _alt_sign_vector ((suNg_vector_V){ -1, +1, -1, +1, -1, +1 })

#define _conj_sign_vector ((suNg_vector_V){ +1, -1, +1, -1, +1, -1 })

typedef struct suNg_V {
    suNg_vector_V c[NG];
} suNg_V;

#define _mul_add(a, b, c) ((a) * (b) + (c))

#define _decompose_suNg_vector_V(reV, imV, inputV)                                            \
    const suNg_vector_V(reV) = __builtin_shufflevector((inputV), (inputV), 0, 0, 2, 2, 4, 4); \
    const suNg_vector_V(imV) = __builtin_shufflevector((inputV), (inputV), 1, 1, 3, 3, 5, 5)

#define _re_suNg_vector_V(inputV) __builtin_shufflevector((inputV), (inputV), 0, 0, 2, 2, 4, 4)

#define _im_suNg_vector_V(inputV) __builtin_shufflevector((inputV), (inputV), 1, 1, 3, 3, 5, 5)

#define _select_suNg_vector_V(inputV, i, j, sl, isl)                    \
    sl = __builtin_shufflevector((inputV), (inputV), i, j, i, j, i, j); \
    Isl = _conj_sign_vector * __builtin_shufflevector((inputV), (inputV), j, i, j, i, j, i)

#define _shuffle_suNg_vector_V(resV, inputV) \
    const suNg_vector_V(resV) = __builtin_shufflevector((inputV), (inputV), 1, 0, 3, 2, 5, 4);

#define _alt_mul(U, V) ((_alt_sign_vector * U) * V)

#define _vector_reduce(res, V1, V2, V3)                                               \
    do {                                                                              \
        suNg_vector_V C1 = __builtin_shufflevector((V1), (V2), 0, 1, 6, 7, -1, -1);   \
        C1[4] = (V3)[0];                                                              \
        C1[5] = (V3)[1];                                                              \
        suNg_vector_V C2 = __builtin_shufflevector((V1), (V2), 2, 3, 8, 9, -1, -1);   \
        C2[4] = (V3)[2];                                                              \
        C2[5] = (V3)[3];                                                              \
        suNg_vector_V C3 = __builtin_shufflevector((V1), (V2), 4, 5, 10, 11, -1, -1); \
        C3[4] = (V3)[4];                                                              \
        C3[5] = (V3)[5];                                                              \
        res = C1 + C2 + C3;                                                           \
    } while (0)

#define _prod_V(res, reU, imU, V, invV) suNg_vector_V(res) = _mul_add(reU, V, _alt_mul(imU, invV))
#define _prod_dag_V(res, U, V, IV) suNg_vector_V(res) = _mul_add(_re_suNg_vector_V(U), V, _im_suNg_vector_V(U) * IV)

#define _MVM_3x3C_SIMD_VEC(r, u, s)                               \
    do {                                                          \
        const suNg_vector_V s_V = *(suNg_vector_V *)(&s);         \
        suNg_vector_V *r_V = (suNg_vector_V *)(&r);               \
        const suNg_vector_V u_c0 = *(suNg_vector_V *)&((u).c[0]); \
        const suNg_vector_V u_c1 = *(suNg_vector_V *)&((u).c[3]); \
        const suNg_vector_V u_c2 = *(suNg_vector_V *)&((u).c[6]); \
        _decompose_suNg_vector_V(reS, imS, s_V);                  \
        _shuffle_suNg_vector_V(invC0, u_c0);                      \
        _prod_V(r1, reS, imS, u_c0, invC0);                       \
        _shuffle_suNg_vector_V(invC1, u_c1);                      \
        _prod_V(r2, reS, imS, u_c1, invC1);                       \
        _shuffle_suNg_vector_V(invC2, u_c2);                      \
        _prod_V(r3, reS, imS, u_c2, invC2);                       \
        _vector_reduce(*r_V, r1, r2, r3);                         \
    } while (0)

#define _MTVM_3x3C_SIMD_VEC(r, u, s)                              \
    do {                                                          \
        const suNg_vector_V s_V = *(suNg_vector_V *)(&s);         \
        suNg_vector_V *r_V = (suNg_vector_V *)(&r);               \
        const suNg_vector_V u_c0 = *(suNg_vector_V *)&((u).c[0]); \
        const suNg_vector_V u_c1 = *(suNg_vector_V *)&((u).c[3]); \
        const suNg_vector_V u_c2 = *(suNg_vector_V *)&((u).c[6]); \
        suNg_vector_V(sl);                                        \
        suNg_vector_V(Isl);                                       \
        _select_suNg_vector_V(s_V, 0, 1, sl, isl);                \
        _prod_dag_V(r1, (u_c0), sl, Isl);                         \
        _select_suNg_vector_V(s_V, 2, 3, sl, isl);                \
        _prod_dag_V(r2, (u_c1), sl, Isl);                         \
        _select_suNg_vector_V(s_V, 4, 5, sl, isl);                \
        _prod_dag_V(r3, (u_c2), sl, Isl);                         \
        *r_V = r1 + r2 + r3;                                      \
    } while (0)

#define _double_MVM_3x3C_SIMD_VEC(r1, r2, u, s1, s2)              \
    do {                                                          \
        const suNg_vector_V s1_V = *(suNg_vector_V *)(&s1);       \
        const suNg_vector_V s2_V = *(suNg_vector_V *)(&s2);       \
        suNg_vector_V *r1_V = (suNg_vector_V *)(&r1);             \
        suNg_vector_V *r2_V = (suNg_vector_V *)(&r2);             \
        const suNg_vector_V u_c0 = *(suNg_vector_V *)&((u).c[0]); \
        const suNg_vector_V u_c1 = *(suNg_vector_V *)&((u).c[3]); \
        const suNg_vector_V u_c2 = *(suNg_vector_V *)&((u).c[6]); \
        _decompose_suNg_vector_V(reS1, imS1, (s1_V));             \
        _decompose_suNg_vector_V(reS2, imS2, (s2_V));             \
        _shuffle_suNg_vector_V(invC0, (u_c0));                    \
        _prod_V(r11, reS1, imS1, (u_c0), invC0);                  \
        _prod_V(r21, reS2, imS2, (u_c0), invC0);                  \
        _shuffle_suNg_vector_V(invC1, (u_c1));                    \
        _prod_V(r12, reS1, imS1, (u_c1), invC1);                  \
        _prod_V(r22, reS2, imS2, (u_c1), invC1);                  \
        _shuffle_suNg_vector_V(invC2, (u_c2));                    \
        _prod_V(r13, reS1, imS1, (u_c2), invC2);                  \
        _prod_V(r23, reS2, imS2, (u_c2), invC2);                  \
        _vector_reduce(*r1_V, r11, r12, r13);                     \
        _vector_reduce(*r2_V, r21, r22, r23);                     \
    } while (0)

#define _double_MTVM_3x3C_SIMD_VEC(r1, r2, u, s1, s2)             \
    do {                                                          \
        const suNg_vector_V s1_V = *(suNg_vector_V *)(&s1);       \
        const suNg_vector_V s2_V = *(suNg_vector_V *)(&s2);       \
        suNg_vector_V *r1_V = (suNg_vector_V *)(&r1);             \
        suNg_vector_V *r2_V = (suNg_vector_V *)(&r2);             \
        const suNg_vector_V u_c0 = *(suNg_vector_V *)&((u).c[0]); \
        const suNg_vector_V u_c1 = *(suNg_vector_V *)&((u).c[3]); \
        const suNg_vector_V u_c2 = *(suNg_vector_V *)&((u).c[6]); \
        suNg_vector_V(sl);                                        \
        suNg_vector_V(Isl);                                       \
        _select_suNg_vector_V(s1_V, 0, 1, sl, isl);               \
        _prod_dag_V(r11, (u_c0), sl, Isl);                        \
        _select_suNg_vector_V(s2_V, 0, 1, sl, isl);               \
        _prod_dag_V(r21, (u_c0), sl, Isl);                        \
        _select_suNg_vector_V(s1_V, 2, 3, sl, isl);               \
        _prod_dag_V(r12, (u_c1), sl, Isl);                        \
        _select_suNg_vector_V(s2_V, 2, 3, sl, isl);               \
        _prod_dag_V(r22, (u_c1), sl, Isl);                        \
        _select_suNg_vector_V(s1_V, 4, 5, sl, isl);               \
        _prod_dag_V(r13, (u_c2), sl, Isl);                        \
        _select_suNg_vector_V(s2_V, 4, 5, sl, isl);               \
        _prod_dag_V(r23, (u_c2), sl, Isl);                        \
        *r1_V = r11 + r12 + r13;                                  \
        *r2_V = r21 + r22 + r23;                                  \
    } while (0)

#endif //(NF == 3) || (NG == 3)

#if (NF == 3) && !defined(REPR_IS_REAL)

#undef _suNf_multiply
#define _suNf_multiply(mc, mu, mp)                                  \
    _Generic((mc), suNf_vector                                      \
             : ({ _MVM_3x3C_SIMD_VEC((mc), (mu), (mp)); }), default \
             : ({ _suNf_multiply_default(mc, mu, mp); }))
#undef _suNf_inverse_multiply
#define _suNf_inverse_multiply(mc, mu, mp)                           \
    _Generic((mc), suNf_vector                                       \
             : ({ _MTVM_3x3C_SIMD_VEC((mc), (mu), (mp)); }), default \
             : ({ _suNf_inverse_multiply_default(mc, mu, mp); }))

#undef _suNf_double_multiply
#define _suNf_double_multiply(mc, mc2, mu, mp, mp2)                                      \
    _Generic((mc), suNf_vector                                                           \
             : ({ _double_MVM_3x3C_SIMD_VEC((mc), (mc2), (mu), (mp), (mp2)); }), default \
             : ({ _suNf_double_multiply_default(mc, mc2, mu, mp, mp2); }))

#undef _suNf_double_inverse_multiply
#define _suNf_double_inverse_multiply(mc, mc2, mu, mp, mp2)                               \
    _Generic((mc), suNf_vector                                                            \
             : ({ _double_MTVM_3x3C_SIMD_VEC((mc), (mc2), (mu), (mp), (mp2)); }), default \
             : ({ _suNf_double_inverse_multiply_default(mc, mc2, mu, mp, mp2); }))

#endif //(NF == 3) && !defined(REPR_IS_REAL)

#if (NG == 3)
#undef _suNg_multiply
#define _suNg_multiply(mc, mu, mp)                                  \
    _Generic((mc), suNg_vector                                      \
             : ({ _MVM_3x3C_SIMD_VEC((mc), (mu), (mp)); }), default \
             : ({ _suNg_multiply_default(mc, mu, mp); }))

#undef _suNg_inverse_multiply
#define _suNg_inverse_multiply(mc, mu, mp)                           \
    _Generic((mc), suNg_vector                                       \
             : ({ _MTVM_3x3C_SIMD_VEC((mc), (mu), (mp)); }), default \
             : ({ _suNg_inverse_multiply_default(mc, mu, mp); }))

#undef _suNg_double_multiply
#define _suNg_double_multiply(mc, mc2, mu, mp, mp2)                                      \
    _Generic((mc), suNg_vector                                                           \
             : ({ _double_MVM_3x3C_SIMD_VEC((mc), (mc2), (mu), (mp), (mp2)); }), default \
             : ({ _suNg_double_multiply_default(mc, mc2, mu, mp, mp2); }))

#undef _suNg_double_inverse_multiply
#define _suNg_double_inverse_multiply(mc, mc2, mu, mp, mp2)                               \
    _Generic((mc), suNg_vector                                                            \
             : ({ _double_MTVM_3x3C_SIMD_VEC((mc), (mc2), (mu), (mp), (mp2)); }), default \
             : ({ _suNg_double_inverse_multiply_default(mc, mc2, mu, mp, mp2); }))

#endif //(NG == 3)

#endif // SIMD_VECTOR_HIREP

#endif // SIMD_VECTOR_HIREP_H
