/*******************************************************************************
*
* Consistency in linear algebra
*
*******************************************************************************/

#include "libhr.h"

int errors = 0;

static int n_geometries = 3;
static geometry_descriptor *geometries[3] = { &glattice, &glat_even, &glat_odd };
static char *desc[3] = { "glattice", "glat_even", "glat_odd" };

static hr_complex c = 2.0 + 1.0 * I;
static double d = 2.43782;

static const int niter = 1;
static const int ninputs = 3;

#define _TEST_BASE(_FIELD_TYPE, _max_geom, _prec)                                                                                                                 \
    do {                                                                                                                                                          \
        for (int geom_idx = 0; geom_idx < _max_geom; ++geom_idx) {                                                                                                \
            _FIELD_TYPE *in = alloc(in, 2, geometries[geom_idx]);                                                                                                 \
            _FIELD_TYPE *out = in + 1;                                                                                                                            \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST BASE", 0,                                                                                                                               \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                /* copy */                                                                                                                                        \
                _TEST_CPU_INV_OP(errors, "copy", 1, in, out, copy(out, in);, "TEST BASE", _prec);                                                                 \
                                                                                                                                                                  \
                /* zero */                                                                                                                                        \
                _TEST_CPU_INV_OP(errors, "zero", 1, in, out, zero(out); zero(in);, "TEST BASE", _prec);                                                           \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

#define _TEST_BASE_OP(_FIELD_TYPE, _max_geom, _prec)                                                                                                              \
    do {                                                                                                                                                          \
        for (int geom_idx = 0; geom_idx < _max_geom; ++geom_idx) {                                                                                                \
            _FIELD_TYPE *in = alloc(in, ninputs + 2, geometries[geom_idx]);                                                                                       \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST BASE OP", 0,                                                                                                                            \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                /* add and sub */                                                                                                                                 \
                _TEST_CPU_INV_OP(errors, "s1 + s2 - s2 = s1", 2, in, in + 3, add(&in[2], &in[0], &in[1]);                                                         \
                                 sub(&in[3], &in[2], &in[1]);, "TEST BASE OP", _prec);                                                                            \
                                                                                                                                                                  \
                /* mul */                                                                                                                                         \
                _TEST_CPU_INV_OP(errors, "s1 * rho / rho = s1", 1, in, in + 2, mul(&in[1], d, &in[0]);                                                            \
                                 mul(&in[2], 1. / d, &in[1]);, "TEST BASE OP", _prec);                                                                            \
                                                                                                                                                                  \
                /* mulc */                                                                                                                                        \
                _TEST_CPU_INV_OP(errors, "s1 * c / c = s1", 1, in, in + 2, mulc(&in[1], c, &in[0]);                                                               \
                                 mulc(&in[2], 1. / c, &in[1]);, "TEST BASE OP", _prec);                                                                           \
                                                                                                                                                                  \
                /* add and sub_assign */                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "s3 = s1 + s2, s3 -= s2", 2, in, in + 2, add(&in[2], &in[0], &in[1]);                                                    \
                                 sub_assign(&in[2], &in[1]);, "TEST BASE OP", _prec);                                                                             \
                                                                                                                                                                  \
                /* sub and add_assign */                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "s3 = s1 - s2, s3 += s2", 2, in, in + 2, sub(&in[2], &in[0], &in[1]);                                                    \
                                 add_assign(&in[2], &in[1]);, "TEST BASE OP", _prec);                                                                             \
                                                                                                                                                                  \
                /* minus */                                                                                                                                       \
                _TEST_CPU_INV_OP(errors, "s1 = -(-s1)", 2, in, in + 2, minus(&in[1], &in[0]); minus(&in[2], &in[1]);                                              \
                                 , "TEST BASE OP", _prec);                                                                                                        \
                                                                                                                                                                  \
                /* mul_add_assign */                                                                                                                              \
                _TEST_CPU_INV_OP(errors, "s2 += s1 * rho, s1 -= s1 * rho", 1, in, in + 2, copy(&in[2], &in[0]);                                                   \
                                 mul_add_assign(&in[2], d, &in[0]); mul(&in[1], d, &in[0]); sub_assign(&in[2], &in[1]);                                           \
                                 , "TEST BASE OP", _prec);                                                                                                        \
                                                                                                                                                                  \
                /* mulc_add_assign */                                                                                                                             \
                _TEST_CPU_INV_OP(errors, "s2 += s1 * c, s1 -= s1 * c", 1, in, in + 2, copy(&in[2], &in[0]);                                                       \
                                 mulc_add_assign(&in[2], c, &in[0]); mulc(&in[1], c, &in[0]); sub_assign(&in[2], &in[1]);                                         \
                                 , "TEST BASE OP", _prec);                                                                                                        \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

#define _TEST_RED(_FIELD_TYPE, _max_geom, _prec)                                                                                                                  \
    do {                                                                                                                                                          \
        for (int geom_idx = 0; geom_idx < _max_geom; ++geom_idx) {                                                                                                \
            _FIELD_TYPE *in = alloc(in, 3, geometries[geom_idx]);                                                                                                 \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST RED", 0,                                                                                                                                \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                /* sqnorm */                                                                                                                                      \
                _TEST_RED_INV_OP(errors, "sqnorm", 1, in, in + 1, double abs1 = sqnorm(&in[0]); mul(&in[1], d, &in[0]);                                           \
                                 double abs2 = sqnorm(&in[1]) / d / d;, "TEST RED", _prec);                                                                       \
                                                                                                                                                                  \
                /* prod */                                                                                                                                        \
                _TEST_RED_INV_OP(errors, "prod", 2, in, in + 1, double abs1 = prod(&in[0], &in[1]);                                                               \
                                 double abs2 = prod(&in[1], &in[0]);, "TEST RED", _prec);                                                                         \
                _TEST_RED_INV_OP(errors, "prod", 2, in, in + 1, double abs1 = prod(&in[0], &in[1]); mul(&in[2], d, &in[1]);                                       \
                                 double abs2 = prod(&in[0], &in[2]) / d;, "TEST RED", _prec);                                                                     \
                                                                                                                                                                  \
                /* max */                                                                                                                                         \
                _TEST_RED_INV_OP(errors, "max", 2, in, in + 1, zero(&in[0]); double abs1 = max(&in[0]) + 0.1;                                                     \
                                 double abs2 = 0.1;, "TEST RED", _prec);                                                                                          \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

#define _TEST_BASE_OP_REAL(_FIELD_TYPE, _max_geom, _prec)                                                                                                         \
    do {                                                                                                                                                          \
        for (int geom_idx = 0; geom_idx < _max_geom; ++geom_idx) {                                                                                                \
            _FIELD_TYPE *in = alloc(in, ninputs + 2, geometries[geom_idx]);                                                                                       \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST BASE OP", 0,                                                                                                                            \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                /* add and sub */                                                                                                                                 \
                _TEST_CPU_INV_OP(errors, "s1 + s2 - s2 = s1", 2, in, in + 3, add(&in[2], &in[0], &in[1]);                                                         \
                                 sub(&in[3], &in[2], &in[1]);, "TEST BASE OP", _prec);                                                                            \
                                                                                                                                                                  \
                /* mul */                                                                                                                                         \
                _TEST_CPU_INV_OP(errors, "s1 * rho / rho = s1", 1, in, in + 2, mul(&in[1], d, &in[0]);                                                            \
                                 mul(&in[2], 1. / d, &in[1]);, "TEST BASE OP", _prec);                                                                            \
                                                                                                                                                                  \
                /* add and sub_assign */                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "s3 = s1 + s2, s3 -= s2", 2, in, in + 2, add(&in[2], &in[0], &in[1]);                                                    \
                                 sub_assign(&in[2], &in[1]);, "TEST BASE OP", _prec);                                                                             \
                                                                                                                                                                  \
                /* sub and add_assign */                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "s3 = s1 - s2, s3 += s2", 2, in, in + 2, sub(&in[2], &in[0], &in[1]);                                                    \
                                 add_assign(&in[2], &in[1]);, "TEST BASE OP", _prec);                                                                             \
                                                                                                                                                                  \
                /* minus */                                                                                                                                       \
                _TEST_CPU_INV_OP(errors, "s1 = -(-s1)", 2, in, in + 2, minus(&in[1], &in[0]); minus(&in[2], &in[1]);                                              \
                                 , "TEST BASE OP", _prec);                                                                                                        \
                                                                                                                                                                  \
                /* mul_add_assign */                                                                                                                              \
                _TEST_CPU_INV_OP(errors, "s2 += s1 * rho, s1 -= s1 * rho", 1, in, in + 2, copy(&in[2], &in[0]);                                                   \
                                 mul_add_assign(&in[2], d, &in[0]); mul(&in[1], d, &in[0]); sub_assign(&in[2], &in[1]);                                           \
                                 , "TEST BASE OP", _prec);                                                                                                        \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

#define _TEST_GAMMA_AND_LC(_FIELD_TYPE, _prec)                                                                                                                    \
    do {                                                                                                                                                          \
        _FIELD_TYPE *in;                                                                                                                                          \
                                                                                                                                                                  \
        for (int geom_idx = 0; geom_idx < n_geometries; geom_idx++) {                                                                                             \
            hr_complex c = 2.0 + 1.0 * I;                                                                                                                         \
            double d = 2.43782;                                                                                                                                   \
            in = alloc(in, ninputs + 2, geometries[geom_idx]);                                                                                                    \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST GAMMA", 0,                                                                                                                              \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                /* g5 */                                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "g5*g5", 1, in, in + 2, g5(&in[1], &in[0]); g5(&in[2], &in[1]);                                                          \
                                 , "TEST GAMMA", _prec);                                                                                                          \
                                                                                                                                                                  \
                /* g1 */                                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "g0*g0", 1, in, in + 2, g0(&in[1], &in[0]); g0(&in[2], &in[1]);                                                          \
                                 , "TEST GAMMA", _prec);                                                                                                          \
                                                                                                                                                                  \
                /* g2 */                                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "g1*g1", 1, in, in + 2, g1(&in[1], &in[0]); g1(&in[2], &in[1]);                                                          \
                                 , "TEST GAMMA", _prec);                                                                                                          \
                                                                                                                                                                  \
                /* g3 */                                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "g3*g3", 1, in, in + 2, g3(&in[1], &in[0]); g3(&in[2], &in[1]);                                                          \
                                 , "TEST GAMMA", _prec);                                                                                                          \
                                                                                                                                                                  \
                /* g5 assign */                                                                                                                                   \
                _TEST_CPU_INV_OP(errors, "g5_assign", 1, in, in + 1, g5(&in[1], &in[0]); g5_assign(&in[0]);                                                       \
                                 , "TEST GAMMA", _prec);                                                                                                          \
                                                                                                                                                                  \
                /* g5 = g0*g1*g2*g3 */                                                                                                                            \
                _TEST_CPU_INV_OP(errors, "g5=g0*g1*g2*g3", 1, in, in + 2, g3(&in[1], &in[0]); g2(&in[2], &in[1]);                                                 \
                                 g1(&in[1], &in[2]); g0(&in[2], &in[1]); g5_assign(&in[0]);, "TEST GAMMA", _prec);                                                \
                                                                                                                                                                  \
                /* lc */                                                                                                                                          \
                _TEST_CPU_INV_OP(errors, "lc", 2, in, in + 3, lc(&in[2], d, &in[0], d, &in[1]);                                                                   \
                                 mul_add_assign(&in[2], -d, &in[1]); mul(&in[3], 1. / d, &in[2]);, "TEST GAMMA", _prec);                                          \
                                                                                                                                                                  \
                /* lc add assign*/                                                                                                                                \
                _TEST_CPU_INV_OP(errors, "lc_add_assign", 3, in, in + 3, copy(&in[3], &in[0]);                                                                    \
                                 lc_add_assign(&in[3], d, &in[1], d, &in[2]); mul_add_assign(&in[3], -d, &in[1]);                                                 \
                                 mul_add_assign(&in[3], -d, &in[2]);, "TEST GAMMA", _prec);                                                                       \
                                                                                                                                                                  \
                /* clc */                                                                                                                                         \
                _TEST_CPU_INV_OP(errors, "clc", 2, in, in + 3, clc(&in[2], c, &in[0], c, &in[1]);                                                                 \
                                 mulc_add_assign(&in[2], -c, &in[1]); mulc(&in[3], 1. / c, &in[2]);, "TEST GAMMA", _prec);                                        \
                                                                                                                                                                  \
                /* clc add assign */                                                                                                                              \
                _TEST_CPU_INV_OP(errors, "clc_add_assign", 3, in, in + 3, copy(&in[3], &in[0]);                                                                   \
                                 clc_add_assign(&in[3], c, &in[1], c, &in[2]); mulc_add_assign(&in[3], -c, &in[1]);                                               \
                                 mulc_add_assign(&in[3], -c, &in[2]);, "TEST GAMMA", _prec);                                                                      \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

int main(int argc, char *argv[]) {
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    rlxs_init(1, 1234);

    /* Base */
    _TEST_BASE(spinor_field, 3, EPSILON_TEST);
    _TEST_BASE(spinor_field_flt, 3, EPSILON_FLT_TEST);
    _TEST_BASE(scalar_field, 1, EPSILON_TEST); //why not 3?(SAM)
    _TEST_BASE(suNg_field, 1, EPSILON_TEST);
    _TEST_BASE(suNf_field, 1, EPSILON_TEST);
    _TEST_BASE(suNg_field_flt, 1, EPSILON_FLT_TEST);
    _TEST_BASE(suNf_field_flt, 1, EPSILON_FLT_TEST);
#ifdef REPR_IS_REAL
    _TEST_BASE(suNfc_field, 1, EPSILON_TEST);
#endif
    _TEST_BASE(suNg_av_field, 1, EPSILON_TEST);
    _TEST_BASE(gtransf, 1, EPSILON_TEST);
    _TEST_BASE(clover_term, 1, EPSILON_TEST);
    _TEST_BASE(clover_force, 1, EPSILON_TEST);
    _TEST_BASE(staple_field, 1, EPSILON_TEST);

    /* Base operations */
    _TEST_BASE_OP(spinor_field, 3, EPSILON_TEST);
    _TEST_BASE_OP(spinor_field_flt, 3, EPSILON_FLT_TEST);
    _TEST_BASE_OP_REAL(scalar_field, 1, EPSILON_TEST);
    _TEST_BASE_OP(suNg_field, 1, EPSILON_TEST);
    _TEST_BASE_OP(suNg_field_flt, 1, EPSILON_FLT_TEST);
#ifdef REPR_IS_REAL
    _TEST_BASE_OP(suNfc_field, 1, EPSILON_TEST);
    _TEST_BASE_OP_REAL(suNf_field, 1, EPSILON_TEST);
    _TEST_BASE_OP_REAL(suNf_field_flt, 1, EPSILON_FLT_TEST);
    _TEST_BASE_OP_REAL(clover_force, 1, EPSILON_TEST);
#else
    _TEST_BASE_OP(suNf_field, 1, EPSILON_TEST);
    _TEST_BASE_OP(suNf_field_flt, 1, EPSILON_TEST);
    _TEST_BASE_OP(clover_force, 1, EPSILON_TEST);
#endif
    _TEST_BASE_OP_REAL(suNg_av_field, 1, EPSILON_TEST);
    _TEST_BASE_OP(gtransf, 1, EPSILON_TEST);
    _TEST_BASE_OP(clover_term, 1, EPSILON_TEST);
    _TEST_BASE_OP(staple_field, 1, EPSILON_TEST);

    /* Reductions */
    _TEST_RED(spinor_field, 3, EPSILON_TEST);
    _TEST_RED(spinor_field_flt, 3, EPSILON_FLT_TEST);
    _TEST_RED(scalar_field, 1, EPSILON_TEST);
    _TEST_RED(suNg_field, 1, EPSILON_TEST);
    _TEST_RED(suNf_field, 1, EPSILON_TEST);
    _TEST_RED(suNg_field_flt, 1, EPSILON_FLT_TEST);
    _TEST_RED(suNf_field_flt, 1, EPSILON_FLT_TEST);
#ifdef REPR_IS_REAL
    _TEST_RED(suNfc_field, 1, EPSILON_TEST);
#endif
    _TEST_RED(suNg_av_field, 1, EPSILON_TEST);
    _TEST_RED(gtransf, 1, EPSILON_TEST);
    _TEST_RED(clover_term, 1, EPSILON_TEST);
    _TEST_RED(clover_force, 1, EPSILON_TEST);
    _TEST_RED(staple_field, 1, EPSILON_TEST);

    /* Gamma and linear combinations */
    _TEST_GAMMA_AND_LC(spinor_field, EPSILON_TEST);
    _TEST_GAMMA_AND_LC(spinor_field_flt, EPSILON_FLT_TEST);

    finalize_process();
    return errors;
}
