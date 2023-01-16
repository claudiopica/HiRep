#ifndef AVX2_HIREP_H
#define AVX2_HIREP_H

#ifdef AVX2_HIREP
#include <immintrin.h>
#include "io.h"
#if (NF == 3) || (NG == 3)

#define _MVM_3x3C_AVX2(mc, mu, mp)                                                                                                 \
  do                                                                                                                               \
  {                                                                                                                                \
    __m256d temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16; \
    __m128d chi_3rd;                                                                                                               \
    temp1 = _mm256_loadu_pd((double *)(&mu));                                                                                      \
    temp2 = _mm256_shuffle_pd(temp1, temp1, 0b0000);                                                                               \
    temp3 = _mm256_shuffle_pd(temp1, temp1, 0b1111);                                                                               \
    temp8 = _mm256_loadu_pd((double *)(&mu) + 2);                                                                                  \
    temp1 = _mm256_loadu_pd((double *)(&mu) + 6);                                                                                  \
    temp4 = _mm256_shuffle_pd(temp1, temp1, 0b0000);                                                                               \
    temp5 = _mm256_shuffle_pd(temp1, temp1, 0b1111);                                                                               \
    temp9 = _mm256_loadu_pd((double *)(&mu) + 8);                                                                                  \
    temp1 = _mm256_loadu_pd((double *)(&mu) + 12);                                                                                 \
    temp6 = _mm256_shuffle_pd(temp1, temp1, 0b0000);                                                                               \
    temp7 = _mm256_shuffle_pd(temp1, temp1, 0b1111);                                                                               \
    temp1 = _mm256_loadu_pd((double *)(&mp));                                                                                      \
    temp14 = _mm256_shuffle_pd(temp1, temp1, 0b0101);                                                                              \
    temp16 = _mm256_loadu_pd((double *)(&mp) + 2);                                                                                 \
    temp3 = _mm256_mul_pd(temp3, temp14);                                                                                          \
    temp2 = _mm256_fmaddsub_pd(temp2, temp1, temp3);                                                                               \
    temp5 = _mm256_mul_pd(temp5, temp14);                                                                                          \
    temp4 = _mm256_fmaddsub_pd(temp4, temp1, temp5);                                                                               \
    temp3 = _mm256_permute2f128_pd(temp4, temp2, 2);                                                                               \
    temp2 = _mm256_permute2f128_pd(temp2, temp2, 1);                                                                               \
    temp2 = _mm256_blend_pd(temp2, temp4, 12);                                                                                     \
    temp2 = _mm256_add_pd(temp3, temp2);                                                                                           \
    temp8 = _mm256_permute2f128_pd(temp8, temp8, 1);                                                                               \
    temp8 = _mm256_blend_pd(temp8, temp9, 12);                                                                                     \
    temp9 = _mm256_permute2f128_pd(temp16, temp16, 1);                                                                             \
    temp4 = _mm256_blend_pd(temp9, temp16, 12);                                                                                    \
    temp10 = _mm256_shuffle_pd(temp8, temp8, 0b0000);                                                                              \
    temp13 = _mm256_shuffle_pd(temp8, temp8, 0b1111);                                                                              \
    temp15 = _mm256_shuffle_pd(temp4, temp4, 0b0101);                                                                              \
    temp13 = _mm256_mul_pd(temp13, temp15);                                                                                        \
    temp10 = _mm256_fmaddsub_pd(temp10, temp4, temp13);                                                                            \
    temp2 = _mm256_add_pd(temp2, temp10);                                                                                          \
    temp7 = _mm256_mul_pd(temp7, temp14);                                                                                          \
    temp6 = _mm256_fmaddsub_pd(temp6, temp1, temp7);                                                                               \
    temp10 = _mm256_permute2f128_pd(temp6, temp6, 1);                                                                              \
    temp1 = _mm256_add_pd(temp10, temp6);                                                                                          \
    temp6 = _mm256_loadu_pd((double *)(&mu) + 14);                                                                                 \
    temp10 = _mm256_permute2f128_pd(temp6, temp6, 1);                                                                              \
    temp10 = _mm256_blend_pd(temp10, temp6, 12);                                                                                   \
    temp11 = _mm256_shuffle_pd(temp10, temp10, 0b0000);                                                                            \
    temp12 = _mm256_shuffle_pd(temp10, temp10, 0b1111);                                                                            \
    temp12 = _mm256_mul_pd(temp12, temp15);                                                                                        \
    temp11 = _mm256_fmaddsub_pd(temp11, temp4, temp12);                                                                            \
    temp11 = _mm256_add_pd(temp1, temp11);                                                                                         \
    chi_3rd = _mm256_castpd256_pd128(temp11);                                                                                      \
    _mm256_store_pd((double *)(&mc), temp2);                                                                                       \
    _mm_store_pd((double *)(&mc) + 4, chi_3rd);                                                                                    \
  } while (0)

#define _MTVM_3x3C_AVX2(mc, mu, mp)                                \
  do                                                               \
  {                                                                \
    __m256d temp1, temp2, temp3, temp4, temp5, temp6, temp7;       \
    const __m256d simd_mask = _mm256_set_pd(-1.0, 1.0, -1.0, 1.0); \
    __m128d chi_3rd;                                               \
    temp1 = _mm256_loadu_pd((double *)(&mu));                      \
    temp1 = _mm256_mul_pd(temp1, simd_mask);                       \
    temp2 = _mm256_loadu_pd((double *)(&mu) + 6);                  \
    temp2 = _mm256_mul_pd(temp2, simd_mask);                       \
    temp3 = _mm256_permute2f128_pd(temp2, temp1, 2);               \
    temp1 = _mm256_permute2f128_pd(temp1, temp1, 1);               \
    temp1 = _mm256_blend_pd(temp1, temp2, 12);                     \
    temp2 = _mm256_shuffle_pd(temp3, temp3, 0b0101);               \
    temp4 = _mm256_shuffle_pd(temp1, temp1, 0b0101);               \
    temp5 = _mm256_loadu_pd((double *)(&mp));                      \
    temp6 = _mm256_shuffle_pd(temp5, temp5, 0b0000);               \
    temp5 = _mm256_shuffle_pd(temp5, temp5, 0b1111);               \
    temp2 = _mm256_mul_pd(temp5, temp2);                           \
    temp2 = _mm256_fmaddsub_pd(temp6, temp3, temp2);               \
    temp3 = _mm256_mul_pd(temp5, temp4);                           \
    temp1 = _mm256_fmaddsub_pd(temp6, temp1, temp3);               \
    temp3 = _mm256_permute2f128_pd(temp1, temp2, 2);               \
    temp2 = _mm256_permute2f128_pd(temp2, temp2, 1);               \
    temp1 = _mm256_blend_pd(temp2, temp1, 12);                     \
    temp1 = _mm256_add_pd(temp3, temp1);                           \
    temp3 = _mm256_loadu_pd((double *)(&mp) + 2);                  \
    temp2 = _mm256_permute2f128_pd(temp3, temp3, 1);               \
    temp2 = _mm256_blend_pd(temp2, temp3, 12);                     \
    temp3 = _mm256_shuffle_pd(temp2, temp2, 0b0000);               \
    temp2 = _mm256_shuffle_pd(temp2, temp2, 0b1111);               \
    temp4 = _mm256_loadu_pd((double *)(&mu) + 12);                 \
    temp4 = _mm256_mul_pd(temp4, simd_mask);                       \
    temp7 = _mm256_shuffle_pd(temp4, temp4, 0b0101);               \
    temp7 = _mm256_mul_pd(temp2, temp7);                           \
    temp4 = _mm256_fmaddsub_pd(temp3, temp4, temp7);               \
    temp1 = _mm256_add_pd(temp1, temp4);                           \
    temp4 = _mm256_loadu_pd((double *)(&mu) + 2);                  \
    temp4 = _mm256_mul_pd(temp4, simd_mask);                       \
    temp7 = _mm256_loadu_pd((double *)(&mu) + 8);                  \
    temp7 = _mm256_mul_pd(temp7, simd_mask);                       \
    temp4 = _mm256_permute2f128_pd(temp4, temp4, 1);               \
    temp4 = _mm256_blend_pd(temp4, temp7, 12);                     \
    temp7 = _mm256_shuffle_pd(temp4, temp4, 0b0101);               \
    temp5 = _mm256_mul_pd(temp5, temp7);                           \
    temp4 = _mm256_fmaddsub_pd(temp6, temp4, temp5);               \
    temp5 = _mm256_permute2f128_pd(temp4, temp4, 1);               \
    temp5 = _mm256_add_pd(temp5, temp4);                           \
    temp4 = _mm256_loadu_pd((double *)(&mu) + 14);                 \
    temp4 = _mm256_mul_pd(temp4, simd_mask);                       \
    temp6 = _mm256_permute2f128_pd(temp4, temp4, 1);               \
    temp4 = _mm256_blend_pd(temp6, temp4, 12);                     \
    temp6 = _mm256_shuffle_pd(temp4, temp4, 0b0101);               \
    temp2 = _mm256_mul_pd(temp2, temp6);                           \
    temp2 = _mm256_fmaddsub_pd(temp3, temp4, temp2);               \
    temp2 = _mm256_add_pd(temp5, temp2);                           \
    chi_3rd = _mm256_castpd256_pd128(temp2);                       \
    _mm256_store_pd((double *)(&mc), temp1);                       \
    _mm_store_pd((double *)(&mc) + 4, chi_3rd);                    \
  } while (0)

#define _double_MVM_3x3C_AVX2(mc, mc2, mu, mp, mp2)                                                                                        \
  do                                                                                                                                       \
  {                                                                                                                                        \
    __m256d temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17; \
    __m128d chi_3rd, chi2_3rd;                                                                                                             \
    temp1 = _mm256_loadu_pd((double *)(&mu));                                                                                              \
    temp6 = _mm256_shuffle_pd(temp1, temp1, 0b0000);                                                                                       \
    temp1 = _mm256_shuffle_pd(temp1, temp1, 0b1111);                                                                                       \
    temp2 = _mm256_loadu_pd((double *)(&mu) + 6);                                                                                          \
    temp7 = _mm256_shuffle_pd(temp2, temp2, 0b0000);                                                                                       \
    temp2 = _mm256_shuffle_pd(temp2, temp2, 0b1111);                                                                                       \
    temp3 = _mm256_loadu_pd((double *)(&mu) + 12);                                                                                         \
    temp8 = _mm256_shuffle_pd(temp3, temp3, 0b0000);                                                                                       \
    temp3 = _mm256_shuffle_pd(temp3, temp3, 0b1111);                                                                                       \
    temp4 = _mm256_loadu_pd((double *)(&mp));                                                                                              \
    temp9 = _mm256_shuffle_pd(temp4, temp4, 0b0101);                                                                                       \
    temp5 = _mm256_loadu_pd((double *)(&mp2));                                                                                             \
    temp10 = _mm256_shuffle_pd(temp5, temp5, 0b0101);                                                                                      \
    temp12 = _mm256_mul_pd(temp1, temp9);                                                                                                  \
    temp11 = _mm256_fmaddsub_pd(temp6, temp4, temp12);                                                                                     \
    temp13 = _mm256_mul_pd(temp2, temp9);                                                                                                  \
    temp12 = _mm256_fmaddsub_pd(temp7, temp4, temp13);                                                                                     \
    temp13 = _mm256_permute2f128_pd(temp12, temp11, 2);                                                                                    \
    temp11 = _mm256_permute2f128_pd(temp11, temp11, 1);                                                                                    \
    temp11 = _mm256_blend_pd(temp11, temp12, 12);                                                                                          \
    temp11 = _mm256_add_pd(temp13, temp11);                                                                                                \
    temp12 = _mm256_loadu_pd((double *)(&mu) + 2);                                                                                         \
    temp12 = _mm256_permute2f128_pd(temp12, temp12, 1);                                                                                    \
    temp13 = _mm256_loadu_pd((double *)(&mu) + 8);                                                                                         \
    temp12 = _mm256_blend_pd(temp12, temp13, 12);                                                                                          \
    temp13 = _mm256_loadu_pd((double *)(&mp) + 2);                                                                                         \
    temp16 = _mm256_permute2f128_pd(temp13, temp13, 1);                                                                                    \
    temp13 = _mm256_blend_pd(temp16, temp13, 12);                                                                                          \
    temp15 = _mm256_shuffle_pd(temp12, temp12, 0b0000);                                                                                    \
    temp12 = _mm256_shuffle_pd(temp12, temp12, 0b1111);                                                                                    \
    temp14 = _mm256_shuffle_pd(temp13, temp13, 0b0101);                                                                                    \
    temp14 = _mm256_mul_pd(temp12, temp14);                                                                                                \
    temp13 = _mm256_fmaddsub_pd(temp15, temp13, temp14);                                                                                   \
    temp11 = _mm256_add_pd(temp11, temp13);                                                                                                \
    temp1 = _mm256_mul_pd(temp1, temp10);                                                                                                  \
    temp1 = _mm256_fmaddsub_pd(temp6, temp5, temp1);                                                                                       \
    temp2 = _mm256_mul_pd(temp2, temp10);                                                                                                  \
    temp7 = _mm256_fmaddsub_pd(temp7, temp5, temp2);                                                                                       \
    temp13 = _mm256_permute2f128_pd(temp7, temp1, 2);                                                                                      \
    temp1 = _mm256_permute2f128_pd(temp1, temp1, 1);                                                                                       \
    temp1 = _mm256_blend_pd(temp1, temp7, 12);                                                                                             \
    temp1 = _mm256_add_pd(temp13, temp1);                                                                                                  \
    temp13 = _mm256_loadu_pd((double *)(&mp2) + 2);                                                                                        \
    temp14 = _mm256_permute2f128_pd(temp13, temp13, 1);                                                                                    \
    temp14 = _mm256_blend_pd(temp14, temp13, 12);                                                                                          \
    temp17 = _mm256_shuffle_pd(temp14, temp14, 0b0101);                                                                                    \
    temp12 = _mm256_mul_pd(temp12, temp17);                                                                                                \
    temp12 = _mm256_fmaddsub_pd(temp15, temp14, temp12);                                                                                   \
    temp1 = _mm256_add_pd(temp1, temp12);                                                                                                  \
    temp12 = _mm256_mul_pd(temp3, temp9);                                                                                                  \
    temp4 = _mm256_fmaddsub_pd(temp8, temp4, temp12);                                                                                      \
    temp3 = _mm256_mul_pd(temp3, temp10);                                                                                                  \
    temp3 = _mm256_fmaddsub_pd(temp8, temp5, temp3);                                                                                       \
    temp5 = _mm256_permute2f128_pd(temp3, temp4, 2);                                                                                       \
    temp4 = _mm256_permute2f128_pd(temp4, temp4, 1);                                                                                       \
    temp3 = _mm256_blend_pd(temp4, temp3, 12);                                                                                             \
    temp3 = _mm256_add_pd(temp5, temp3);                                                                                                   \
    temp9 = _mm256_loadu_pd((double *)(&mu) + 14);                                                                                         \
    temp10 = _mm256_permute2f128_pd(temp9, temp9, 1);                                                                                      \
    temp9 = _mm256_blend_pd(temp10, temp9, 12);                                                                                            \
    temp10 = _mm256_shuffle_pd(temp9, temp9, 0b0000);                                                                                      \
    temp12 = _mm256_shuffle_pd(temp9, temp9, 0b1111);                                                                                      \
    temp9 = _mm256_blend_pd(temp16, temp13, 12);                                                                                           \
    temp13 = _mm256_shuffle_pd(temp9, temp9, 0b0101);                                                                                      \
    temp2 = _mm256_mul_pd(temp12, temp13);                                                                                                 \
    temp7 = _mm256_fmaddsub_pd(temp10, temp9, temp2);                                                                                      \
    temp2 = _mm256_add_pd(temp3, temp7);                                                                                                   \
    chi_3rd = _mm256_castpd256_pd128(temp2);                                                                                               \
    chi2_3rd = _mm256_extractf128_pd(temp2, 1);                                                                                            \
    _mm256_store_pd((double *)(&mc), temp11);                                                                                              \
    _mm_store_pd((double *)(&mc) + 4, chi_3rd);                                                                                            \
    _mm256_store_pd((double *)(&mc2), temp1);                                                                                              \
    _mm_store_pd((double *)(&mc2) + 4, chi2_3rd);                                                                                          \
  } while (0)

#define _double_MTVM_3x3C_AVX2(mc, mc2, mu, mp, mp2)                                                                       \
  do                                                                                                                       \
  {                                                                                                                        \
    __m256d temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15; \
    __m128d chi_3rd, chi2_3rd;                                                                                             \
    const __m256d simd_mask = _mm256_set_pd(-1.0, 1.0, -1.0, 1.0);                                                         \
    temp1 = _mm256_loadu_pd((double *)(&mu));                                                                              \
    temp1 = _mm256_mul_pd(temp1, simd_mask);                                                                               \
    temp2 = _mm256_loadu_pd((double *)(&mu) + 6);                                                                          \
    temp2 = _mm256_mul_pd(temp2, simd_mask);                                                                               \
    temp3 = _mm256_permute2f128_pd(temp2, temp1, 2);                                                                       \
    temp1 = _mm256_permute2f128_pd(temp1, temp1, 1);                                                                       \
    temp1 = _mm256_blend_pd(temp1, temp2, 12);                                                                             \
    temp2 = _mm256_shuffle_pd(temp3, temp3, 0b0101);                                                                       \
    temp4 = _mm256_shuffle_pd(temp1, temp1, 0b0101);                                                                       \
    temp5 = _mm256_loadu_pd((double *)(&mp));                                                                              \
    temp6 = _mm256_shuffle_pd(temp5, temp5, 0b0000);                                                                       \
    temp5 = _mm256_shuffle_pd(temp5, temp5, 0b1111);                                                                       \
    temp7 = _mm256_loadu_pd((double *)(&mp2));                                                                             \
    temp8 = _mm256_shuffle_pd(temp7, temp7, 0b0000);                                                                       \
    temp7 = _mm256_shuffle_pd(temp7, temp7, 0b1111);                                                                       \
    temp10 = _mm256_mul_pd(temp5, temp2);                                                                                  \
    temp9 = _mm256_fmaddsub_pd(temp6, temp3, temp10);                                                                      \
    temp11 = _mm256_mul_pd(temp5, temp4);                                                                                  \
    temp10 = _mm256_fmaddsub_pd(temp6, temp1, temp11);                                                                     \
    temp11 = _mm256_permute2f128_pd(temp10, temp9, 2);                                                                     \
    temp9 = _mm256_permute2f128_pd(temp9, temp9, 1);                                                                       \
    temp9 = _mm256_blend_pd(temp9, temp10, 12);                                                                            \
    temp9 = _mm256_add_pd(temp11, temp9);                                                                                  \
    temp10 = _mm256_loadu_pd((double *)(&mp) + 2);                                                                         \
    temp11 = _mm256_permute2f128_pd(temp10, temp10, 1);                                                                    \
    temp10 = _mm256_blend_pd(temp11, temp10, 12);                                                                          \
    temp12 = _mm256_shuffle_pd(temp10, temp10, 0b0000);                                                                    \
    temp10 = _mm256_shuffle_pd(temp10, temp10, 0b1111);                                                                    \
    temp13 = _mm256_loadu_pd((double *)(&mu) + 12);                                                                        \
    temp13 = _mm256_mul_pd(temp13, simd_mask);                                                                             \
    temp15 = _mm256_shuffle_pd(temp13, temp13, 0b0101);                                                                    \
    temp10 = _mm256_mul_pd(temp10, temp15);                                                                                \
    temp10 = _mm256_fmaddsub_pd(temp12, temp13, temp10);                                                                   \
    temp9 = _mm256_add_pd(temp9, temp10);                                                                                  \
    temp2 = _mm256_mul_pd(temp7, temp2);                                                                                   \
    temp2 = _mm256_fmaddsub_pd(temp8, temp3, temp2);                                                                       \
    temp3 = _mm256_mul_pd(temp7, temp4);                                                                                   \
    temp1 = _mm256_fmaddsub_pd(temp8, temp1, temp3);                                                                       \
    temp3 = _mm256_permute2f128_pd(temp1, temp2, 2);                                                                       \
    temp2 = _mm256_permute2f128_pd(temp2, temp2, 1);                                                                       \
    temp1 = _mm256_blend_pd(temp2, temp1, 12);                                                                             \
    temp1 = _mm256_add_pd(temp3, temp1);                                                                                   \
    temp2 = _mm256_loadu_pd((double *)(&mp2) + 2);                                                                         \
    temp3 = _mm256_permute2f128_pd(temp2, temp2, 1);                                                                       \
    temp3 = _mm256_blend_pd(temp3, temp2, 12);                                                                             \
    temp4 = _mm256_shuffle_pd(temp3, temp3, 0b0000);                                                                       \
    temp3 = _mm256_shuffle_pd(temp3, temp3, 0b1111);                                                                       \
    temp3 = _mm256_mul_pd(temp3, temp15);                                                                                  \
    temp3 = _mm256_fmaddsub_pd(temp4, temp13, temp3);                                                                      \
    temp1 = _mm256_add_pd(temp1, temp3);                                                                                   \
    temp3 = _mm256_loadu_pd((double *)(&mu) + 2);                                                                          \
    temp3 = _mm256_mul_pd(temp3, simd_mask);                                                                               \
    temp4 = _mm256_loadu_pd((double *)(&mu) + 8);                                                                          \
    temp4 = _mm256_mul_pd(temp4, simd_mask);                                                                               \
    temp3 = _mm256_permute2f128_pd(temp3, temp3, 1);                                                                       \
    temp3 = _mm256_blend_pd(temp3, temp4, 12);                                                                             \
    temp4 = _mm256_shuffle_pd(temp3, temp3, 0b0101);                                                                       \
    temp5 = _mm256_mul_pd(temp5, temp4);                                                                                   \
    temp5 = _mm256_fmaddsub_pd(temp6, temp3, temp5);                                                                       \
    temp4 = _mm256_mul_pd(temp7, temp4);                                                                                   \
    temp3 = _mm256_fmaddsub_pd(temp8, temp3, temp4);                                                                       \
    temp4 = _mm256_permute2f128_pd(temp3, temp5, 2);                                                                       \
    temp14 = _mm256_permute2f128_pd(temp5, temp5, 1);                                                                      \
    temp3 = _mm256_blend_pd(temp14, temp3, 12);                                                                            \
    temp3 = _mm256_add_pd(temp4, temp3);                                                                                   \
    temp2 = _mm256_blend_pd(temp11, temp2, 12);                                                                            \
    temp4 = _mm256_shuffle_pd(temp2, temp2, 0b0000);                                                                       \
    temp2 = _mm256_shuffle_pd(temp2, temp2, 0b1111);                                                                       \
    temp5 = _mm256_loadu_pd((double *)(&mu) + 14);                                                                         \
    temp5 = _mm256_mul_pd(temp5, simd_mask);                                                                               \
    temp6 = _mm256_permute2f128_pd(temp5, temp5, 1);                                                                       \
    temp5 = _mm256_blend_pd(temp6, temp5, 12);                                                                             \
    temp6 = _mm256_shuffle_pd(temp5, temp5, 0b0101);                                                                       \
    temp2 = _mm256_mul_pd(temp2, temp6);                                                                                   \
    temp2 = _mm256_fmaddsub_pd(temp4, temp5, temp2);                                                                       \
    temp2 = _mm256_add_pd(temp3, temp2);                                                                                   \
    chi_3rd = _mm256_castpd256_pd128(temp2);                                                                               \
    chi2_3rd = _mm256_extractf128_pd(temp2, 1);                                                                            \
    _mm256_store_pd((double *)(&mc), temp9);                                                                               \
    _mm_storeu_pd((double *)(&mc) + 4, chi_3rd);                                                                           \
    _mm256_store_pd((double *)(&mc2), temp1);                                                                              \
    _mm_storeu_pd((double *)(&mc2) + 4, chi2_3rd);                                                                         \
  } while (0)

#endif //(NF == 3) || (NG == 3)

#if (NF == 2) || (NG == 2)

#define _MVM_2x2C_AVX2(mc, mu, mp)                    \
  do                                                  \
  {                                                   \
    __m256d temp1, temp2, temp3, temp4, temp5, temp6; \
    temp1 = _mm256_loadu_pd((double *)(&mu));         \
    temp4 = _mm256_shuffle_pd(temp1, temp1, 0b0000);  \
    temp1 = _mm256_shuffle_pd(temp1, temp1, 0b1111);  \
    temp2 = _mm256_loadu_pd((double *)(&mu) + 4);     \
    temp5 = _mm256_shuffle_pd(temp2, temp2, 0b0000);  \
    temp2 = _mm256_shuffle_pd(temp2, temp2, 0b1111);  \
    temp3 = _mm256_loadu_pd((double *)(&mp));         \
    temp6 = _mm256_shuffle_pd(temp3, temp3, 0b0101);  \
    temp1 = _mm256_mul_pd(temp1, temp6);              \
    temp1 = _mm256_fmaddsub_pd(temp4, temp3, temp1);  \
    temp2 = _mm256_mul_pd(temp2, temp6);              \
    temp2 = _mm256_fmaddsub_pd(temp5, temp3, temp2);  \
    temp4 = _mm256_permute2f128_pd(temp2, temp1, 2);  \
    temp1 = _mm256_permute2f128_pd(temp1, temp1, 1);  \
    temp1 = _mm256_blend_pd(temp1, temp2, 12);        \
    temp1 = _mm256_add_pd(temp4, temp1);              \
    _mm256_store_pd((double *)(&mc), temp1);          \
  } while (0)

#define _MTVM_2x2C_AVX2(mc, mu, mp)                                \
  do                                                               \
  {                                                                \
    __m256d temp1, temp2, temp3, temp4, temp5, temp6;              \
    const __m256d simd_mask = _mm256_set_pd(-1.0, 1.0, -1.0, 1.0); \
    temp1 = _mm256_loadu_pd((double *)(&mu));                      \
    temp1 = _mm256_mul_pd(temp1, simd_mask);                       \
    temp2 = _mm256_loadu_pd((double *)(&mu) + 4);                  \
    temp2 = _mm256_mul_pd(temp2, simd_mask);                       \
    temp4 = _mm256_permute2f128_pd(temp2, temp1, 2);               \
    temp1 = _mm256_permute2f128_pd(temp1, temp1, 1);               \
    temp1 = _mm256_blend_pd(temp1, temp2, 12);                     \
    temp2 = _mm256_shuffle_pd(temp4, temp4, 0b0101);               \
    temp5 = _mm256_shuffle_pd(temp1, temp1, 0b0101);               \
    temp3 = _mm256_loadu_pd((double *)(&mp));                      \
    temp6 = _mm256_shuffle_pd(temp3, temp3, 0b0000);               \
    temp3 = _mm256_shuffle_pd(temp3, temp3, 0b1111);               \
    temp2 = _mm256_mul_pd(temp3, temp2);                           \
    temp2 = _mm256_fmaddsub_pd(temp6, temp4, temp2);               \
    temp3 = _mm256_mul_pd(temp3, temp5);                           \
    temp1 = _mm256_fmaddsub_pd(temp6, temp1, temp3);               \
    temp3 = _mm256_permute2f128_pd(temp1, temp2, 2);               \
    temp2 = _mm256_permute2f128_pd(temp2, temp2, 1);               \
    temp1 = _mm256_blend_pd(temp2, temp1, 12);                     \
    temp1 = _mm256_add_pd(temp3, temp1);                           \
    _mm256_store_pd((double *)(&mc), temp1);                       \
  } while (0)

#define _double_MVM_2x2C_AVX2(mc, mc2, mu, mp, mp2)                        \
  do                                                                       \
  {                                                                        \
    __m256d temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9; \
    temp1 = _mm256_loadu_pd((double *)(&mu));                              \
    temp2 = _mm256_loadu_pd((double *)(&mu) + 4);                          \
    temp3 = _mm256_loadu_pd((double *)(&mp));                              \
    temp5 = _mm256_shuffle_pd(temp3, temp3, 0b0101);                       \
    temp4 = _mm256_loadu_pd((double *)(&mp2));                             \
    temp6 = _mm256_shuffle_pd(temp4, temp4, 0b0101);                       \
    temp7 = _mm256_shuffle_pd(temp1, temp1, 0b0000);                       \
    temp1 = _mm256_shuffle_pd(temp1, temp1, 0b1111);                       \
    temp9 = _mm256_mul_pd(temp1, temp5);                                   \
    temp8 = _mm256_fmaddsub_pd(temp7, temp3, temp9);                       \
    temp9 = _mm256_shuffle_pd(temp2, temp2, 0b0000);                       \
    temp2 = _mm256_shuffle_pd(temp2, temp2, 0b1111);                       \
    temp5 = _mm256_mul_pd(temp2, temp5);                                   \
    temp3 = _mm256_fmaddsub_pd(temp9, temp3, temp5);                       \
    temp5 = _mm256_permute2f128_pd(temp3, temp8, 2);                       \
    temp8 = _mm256_permute2f128_pd(temp8, temp8, 1);                       \
    temp3 = _mm256_blend_pd(temp8, temp3, 12);                             \
    temp3 = _mm256_add_pd(temp5, temp3);                                   \
    temp1 = _mm256_mul_pd(temp1, temp6);                                   \
    temp1 = _mm256_fmaddsub_pd(temp7, temp4, temp1);                       \
    temp2 = _mm256_mul_pd(temp2, temp6);                                   \
    temp2 = _mm256_fmaddsub_pd(temp9, temp4, temp2);                       \
    temp4 = _mm256_permute2f128_pd(temp2, temp1, 2);                       \
    temp1 = _mm256_permute2f128_pd(temp1, temp1, 1);                       \
    temp1 = _mm256_blend_pd(temp1, temp2, 12);                             \
    temp1 = _mm256_add_pd(temp4, temp1);                                   \
    _mm256_store_pd((double *)(&mc), temp3);                               \
    _mm256_store_pd((double *)(&mc), temp1);                               \
  } while (0)

#define _double_MTVM_2x2C_AVX2(mc, mc2, mu, mp, mp2)                                       \
  do                                                                                       \
  {                                                                                        \
    __m256d temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11; \
    const __m256d simd_mask = _mm256_set_pd(-1.0, 1.0, -1.0, 1.0);                         \
    temp1 = _mm256_loadu_pd((double *)(&mu));                                              \
    temp1 = _mm256_mul_pd(temp1, simd_mask);                                               \
    temp2 = _mm256_loadu_pd((double *)(&mu) + 4);                                          \
    temp2 = _mm256_mul_pd(temp2, simd_mask);                                               \
    temp4 = _mm256_permute2f128_pd(temp2, temp1, 2);                                       \
    temp1 = _mm256_permute2f128_pd(temp1, temp1, 1);                                       \
    temp1 = _mm256_blend_pd(temp1, temp2, 12);                                             \
    temp2 = _mm256_shuffle_pd(temp4, temp4, 0b0101);                                       \
    temp5 = _mm256_shuffle_pd(temp1, temp1, 0b0101);                                       \
    temp3 = _mm256_loadu_pd((double *)(&mp));                                              \
    temp6 = _mm256_shuffle_pd(temp3, temp3, 0b0000);                                       \
    temp3 = _mm256_shuffle_pd(temp3, temp3, 0b1111);                                       \
    temp7 = _mm256_loadu_pd((double *)(&mp2));                                             \
    temp8 = _mm256_shuffle_pd(temp7, temp7, 0b0000);                                       \
    temp7 = _mm256_shuffle_pd(temp7, temp7, 0b1111);                                       \
    temp9 = _mm256_mul_pd(temp6, temp4);                                                   \
    temp10 = _mm256_mul_pd(temp3, temp2);                                                  \
    temp9 = _mm256_addsub_pd(temp9, temp10);                                               \
    temp3 = _mm256_mul_pd(temp3, temp5);                                                   \
    temp3 = _mm256_addsub_pd(temp6, temp3);                                                \
    temp6 = _mm256_permute2f128_pd(temp3, temp9, 2);                                       \
    temp11 = _mm256_permute2f128_pd(temp9, temp9, 1);                                      \
    temp3 = _mm256_blend_pd(temp11, temp3, 12);                                            \
    temp3 = _mm256_add_pd(temp6, temp3);                                                   \
    temp4 = _mm256_mul_pd(temp8, temp4);                                                   \
    temp2 = _mm256_mul_pd(temp7, temp2);                                                   \
    temp2 = _mm256_addsub_pd(temp4, temp2);                                                \
    temp1 = _mm256_mul_pd(temp8, temp1);                                                   \
    temp4 = _mm256_mul_pd(temp7, temp5);                                                   \
    temp1 = _mm256_addsub_pd(temp1, temp4);                                                \
    temp4 = _mm256_permute2f128_pd(temp1, temp2, 2);                                       \
    temp2 = _mm256_permute2f128_pd(temp2, temp2, 1);                                       \
    temp1 = _mm256_blend_pd(temp2, temp1, 12);                                             \
    temp2 = _mm256_add_pd(temp4, temp1);                                                   \
    _mm256_store_pd((double *)(&mc), temp3);                                               \
    _mm256_store_pd((double *)(&mc2), temp2);                                              \
  } while (0)

#endif //(NF == 2) || (NG == 2)

#if (NF == 3) && !defined(REPR_IS_REAL)

#undef _suNf_multiply
#define _suNf_multiply(mc, mu, mp) _Generic((mc),                                      \
                                            suNf_vector                                \
                                            : ({ _MVM_3x3C_AVX2((mc), (mu), (mp)); }), \
                                              default                                  \
                                            : ({ _suNf_multiply_default(mc, mu, mp); }))
#undef _suNf_inverse_multiply
#define _suNf_inverse_multiply(mc, mu, mp) _Generic((mc),                                       \
                                                    suNf_vector                                 \
                                                    : ({ _MTVM_3x3C_AVX2((mc), (mu), (mp)); }), \
                                                      default                                   \
                                                    : ({ _suNf_inverse_multiply_default(mc, mu, mp); }))

#undef _suNf_double_multiply
#define _suNf_double_multiply(mc, mc2, mu, mp, mp2) _Generic((mc),                                                           \
                                                             suNf_vector                                                     \
                                                             : ({ _double_MVM_3x3C_AVX2((mc), (mc2), (mu), (mp), (mp2)); }), \
                                                               default                                                       \
                                                             : ({ _suNf_double_multiply_default(mc, mc2, mu, mp, mp2); }))

#undef _suNf_double_inverse_multiply
#define _suNf_double_inverse_multiply(mc, mc2, mu, mp, mp2) _Generic((mc),                                                            \
                                                                     suNf_vector                                                      \
                                                                     : ({ _double_MTVM_3x3C_AVX2((mc), (mc2), (mu), (mp), (mp2)); }), \
                                                                       default                                                        \
                                                                     : ({ _suNf_double_inverse_multiply_default(mc, mc2, mu, mp, mp2); }))

#endif //(NF == 3) && !defined(REPR_IS_REAL)

#if (NG == 3)
#undef _suNg_multiply
#define _suNg_multiply(mc, mu, mp) _Generic((mc),                                      \
                                            suNg_vector                                \
                                            : ({ _MVM_3x3C_AVX2((mc), (mu), (mp)); }), \
                                              default                                  \
                                            : ({ _suNg_multiply_default(mc, mu, mp); }))

#undef _suNg_inverse_multiply
#define _suNg_inverse_multiply(mc, mu, mp) _Generic((mc),                                       \
                                                    suNg_vector                                 \
                                                    : ({ _MTVM_3x3C_AVX2((mc), (mu), (mp)); }), \
                                                      default                                   \
                                                    : ({ _suNg_inverse_multiply_default(mc, mu, mp); }))

#undef _suNg_double_multiply
#define _suNg_double_multiply(mc, mc2, mu, mp, mp2) _Generic((mc),                                                           \
                                                             suNg_vector                                                     \
                                                             : ({ _double_MVM_3x3C_AVX2((mc), (mc2), (mu), (mp), (mp2)); }), \
                                                               default                                                       \
                                                             : ({ _suNg_double_multiply_default(mc, mc2, mu, mp, mp2); }))

#undef _suNg_double_inverse_multiply
#define _suNg_double_inverse_multiply(mc, mc2, mu, mp, mp2) _Generic((mc),                                                            \
                                                                     suNg_vector                                                      \
                                                                     : ({ _double_MTVM_3x3C_AVX2((mc), (mc2), (mu), (mp), (mp2)); }), \
                                                                       default                                                        \
                                                                     : ({ _suNg_double_inverse_multiply_default(mc, mc2, mu, mp, mp2); }))

#endif //(NG == 3)

#if (NF == 2) && !defined(REPR_IS_REAL)

#undef _suNf_multiply
#define _suNf_multiply(mc, mu, mp) _Generic((mc),                                      \
                                            suNf_vector                                \
                                            : ({ _MVM_2x2C_AVX2((mc), (mu), (mp)); }), \
                                              default                                  \
                                            : ({ _suNf_multiply_default(mc, mu, mp); }))
#undef _suNf_inverse_multiply
#define _suNf_inverse_multiply(mc, mu, mp) _Generic((mc),                                       \
                                                    suNf_vector                                 \
                                                    : ({ _MTVM_2x2C_AVX2((mc), (mu), (mp)); }), \
                                                      default                                   \
                                                    : ({ _suNf_inverse_multiply_default(mc, mu, mp); }))

#undef _suNf_double_multiply
#define _suNf_double_multiply(mc, mc2, mu, mp, mp2) _Generic((mc),                                                           \
                                                             suNf_vector                                                     \
                                                             : ({ _double_MVM_2x2C_AVX2((mc), (mc2), (mu), (mp), (mp2)); }), \
                                                               default                                                       \
                                                             : ({ _suNf_double_multiply_default(mc, mc2, mu, mp, mp2); }))

#undef _suNf_double_inverse_multiply
#define _suNf_double_inverse_multiply(mc, mc2, mu, mp, mp2) _Generic((mc),                                                            \
                                                                     suNf_vector                                                      \
                                                                     : ({ _double_MTVM_2x2C_AVX2((mc), (mc2), (mu), (mp), (mp2)); }), \
                                                                       default                                                        \
                                                                     : ({ _suNf_double_inverse_multiply_default(mc, mc2, mu, mp, mp2); }))

#endif //(NF == 2) && !defined(REPR_IS_REAL)

#if (NG == 2)
#undef _suNg_multiply
#define _suNg_multiply(mc, mu, mp) _Generic((mc),                                      \
                                            suNg_vector                                \
                                            : ({ _MVM_2x2C_AVX2((mc), (mu), (mp)); }), \
                                              default                                  \
                                            : ({ _suNg_multiply_default(mc, mu, mp); }))

#undef _suNg_inverse_multiply
#define _suNg_inverse_multiply(mc, mu, mp) _Generic((mc),                                       \
                                                    suNg_vector                                 \
                                                    : ({ _MTVM_2x2C_AVX2((mc), (mu), (mp)); }), \
                                                      default                                   \
                                                    : ({ _suNg_inverse_multiply_default(mc, mu, mp); }))

#undef _suNg_double_multiply
#define _suNg_double_multiply(mc, mc2, mu, mp, mp2) _Generic((mc),                                                           \
                                                             suNg_vector                                                     \
                                                             : ({ _double_MVM_2x2C_AVX2((mc), (mc2), (mu), (mp), (mp2)); }), \
                                                               default                                                       \
                                                             : ({ _suNg_double_multiply_default(mc, mc2, mu, mp, mp2); }))

#undef _suNg_double_inverse_multiply
#define _suNg_double_inverse_multiply(mc, mc2, mu, mp, mp2) _Generic((mc),                                                            \
                                                                     suNg_vector                                                      \
                                                                     : ({ _double_MTVM_2x2C_AVX2((mc), (mc2), (mu), (mp), (mp2)); }), \
                                                                       default                                                        \
                                                                     : ({ _suNg_double_inverse_multiply_default(mc, mc2, mu, mp, mp2); }))

#endif //(NG == 2)

#endif // AVX2_HIREP

#endif // AVX2_HIREP_H
