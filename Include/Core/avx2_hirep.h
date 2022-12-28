#ifdef AVX2_HIREP
#define AVX2_HIREP

#include <immintrin.h>
#undef _suNf_double_multiply 
#define _suNf_double_multiply(mc, mc2, mu, mp, mp2)     \
do                                                 \
{                                                  \
__m256d temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17;\
__m128d chi_3rd, chi2_3rd;                        \
temp1 = _mm256_loadu_pd((double *)(&mu));           \
temp6 = _mm256_shuffle_pd(temp1, temp1, 0b0000);  \
temp1 = _mm256_shuffle_pd(temp1, temp1, 0b1111);  \
temp2 = _mm256_loadu_pd((double *)(&mu) + 6);      \
temp7 = _mm256_shuffle_pd(temp2, temp2, 0b0000);  \
temp2 = _mm256_shuffle_pd(temp2, temp2, 0b1111);  \
temp3 = _mm256_loadu_pd((double *)(&mu) + 12);     \
temp8 = _mm256_shuffle_pd(temp3, temp3, 0b0000);  \
temp3 = _mm256_shuffle_pd(temp3, temp3, 0b1111);  \
temp4 = _mm256_loadu_pd((double *)(&mp));          \
temp9 = _mm256_shuffle_pd(temp4, temp4, 0b0101);  \
temp5 = _mm256_loadu_pd((double *)(&mp2));         \
temp10 = _mm256_shuffle_pd(temp5, temp5, 0b0101); \
temp12 = _mm256_mul_pd(temp1, temp9);             \
temp11 = _mm256_fmaddsub_pd(temp6, temp4, temp12);\
temp13 = _mm256_mul_pd(temp2, temp9);             \
temp12 = _mm256_fmaddsub_pd(temp7, temp4, temp13);\
temp13 = _mm256_permute2f128_pd(temp12, temp11, 2);\
temp11 = _mm256_permute2f128_pd(temp11, temp11, 1);\
temp11 = _mm256_blend_pd(temp11, temp12, 12);      \
temp11 = _mm256_add_pd(temp13, temp11);\
temp12 = _mm256_loadu_pd((double *)(&mu) + 2);\
temp12 = _mm256_permute2f128_pd(temp12, temp12, 1);\
temp13 = _mm256_loadu_pd((double *)(&mu) + 8);\
temp12 = _mm256_blend_pd(temp12, temp13, 12);\
temp13 = _mm256_loadu_pd((double *)(&mp) + 2);\
temp16 = _mm256_permute2f128_pd(temp13, temp13, 1);\
temp13 = _mm256_blend_pd(temp16, temp13, 12);\
temp15 = _mm256_shuffle_pd(temp12, temp12, 0b0000);\
temp12 = _mm256_shuffle_pd(temp12, temp12, 0b1111);\
temp14 = _mm256_shuffle_pd(temp13, temp13, 0b0101);\
temp14 = _mm256_mul_pd(temp12, temp14);\
temp13 = _mm256_fmaddsub_pd(temp15, temp13, temp14);\
temp11 = _mm256_add_pd(temp11, temp13);\
temp1 = _mm256_mul_pd(temp1, temp10);\
temp1 = _mm256_fmaddsub_pd(temp6, temp5, temp1);\
temp2 = _mm256_mul_pd(temp2, temp10);\
temp7 = _mm256_fmaddsub_pd(temp7, temp5, temp2);\
temp13 = _mm256_permute2f128_pd(temp7, temp1, 2);\
temp1 = _mm256_permute2f128_pd(temp1, temp1, 1);\
temp1 = _mm256_blend_pd(temp1, temp7, 12);\
temp1 = _mm256_add_pd(temp13, temp1);\
temp13 = _mm256_loadu_pd((double *)(&mp2) + 2);\
temp14 = _mm256_permute2f128_pd(temp13, temp13, 1);\
temp14 = _mm256_blend_pd(temp14, temp13, 12);\
temp17 = _mm256_shuffle_pd(temp14, temp14, 0b0101);\
temp12 = _mm256_mul_pd(temp12, temp17);\
temp12 = _mm256_fmaddsub_pd(temp15, temp14, temp12);\
temp1 = _mm256_add_pd(temp1, temp12);\
temp12 = _mm256_mul_pd(temp3, temp9);\
temp4 = _mm256_fmaddsub_pd(temp8, temp4, temp12);\
temp3 = _mm256_mul_pd(temp3, temp10);\
temp3 = _mm256_fmaddsub_pd(temp8, temp5, temp3);\
temp5 = _mm256_permute2f128_pd(temp3, temp4, 2);\
temp4 = _mm256_permute2f128_pd(temp4, temp4, 1);\
temp3 = _mm256_blend_pd(temp4, temp3, 12);\
temp3 = _mm256_add_pd(temp5, temp3);\
temp9 = _mm256_loadu_pd((double *)(&mu) + 14);\
temp10 = _mm256_permute2f128_pd(temp9, temp9, 1);\
temp9 = _mm256_blend_pd(temp10, temp9, 12);\
temp10 = _mm256_shuffle_pd(temp9, temp9, 0b0000);\
temp12 = _mm256_shuffle_pd(temp9, temp9, 0b1111);\
temp9 = _mm256_blend_pd(temp16, temp13, 12);\
temp13 = _mm256_shuffle_pd(temp9, temp9, 0b0101);\
temp2 = _mm256_mul_pd(temp12, temp13);\
temp7 = _mm256_fmaddsub_pd(temp10, temp9, temp2);\
temp2 = _mm256_add_pd(temp3, temp7);\
chi_3rd = _mm256_castpd256_pd128(temp2);\
chi2_3rd = _mm256_extractf128_pd(temp2, 1);\
_mm256_store_pd((double *)(&mc), temp11);\
_mm_store_pd((double *)(&mc) + 4, chi_3rd);\
_mm256_store_pd((double *)(&mc2), temp1);\
_mm_store_pd((double *)(&mc2) + 4, chi2_3rd);\
} while (0)

#endif

