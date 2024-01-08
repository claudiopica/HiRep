/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Test of modules
*
*******************************************************************************/

#include "libhr.h"

int errors = 0; // count the number of errors during this test unit

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
                _TEST_GPU_OP(errors, "s2+=r*s1", 2, in, in + 1, double r = 10.0; mul_add_assign(&in[1], r, &in[0]);                                               \
                             mul_add_assign_cpu(&in[1], r, &in[0]);, "GPU TEST", _prec);                                                                          \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2+=c*s1", 2, in, in + 1, hr_complex c = 2.0 + 1.0 * I;                                                                     \
                             mulc_add_assign(&in[1], c, &in[0]); mulc_add_assign_cpu(&in[1], c, &in[0]);, "GPU TEST", _prec);                                     \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2=r*s1", 1, in, in + 1, double r = 10.0; mul(&in[1], r, &in[0]);                                                           \
                             mul_cpu(&in[1], r, &in[0]);, "GPU TEST", _prec);                                                                                     \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2=c*s1", 1, in, in + 1, hr_complex c = 2.0 + 1.0 * I; mulc(&in[1], c, &in[0]);                                             \
                             mulc_cpu(&in[1], c, &in[0]);, "GPU TEST", _prec);                                                                                    \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s3=s1+s2", 2, in, in + 2, add(&in[2], &in[0], &in[1]); add_cpu(&in[2], &in[0], &in[1]);                                     \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s3=s1-s2", 2, in, in + 2, sub(&in[2], &in[0], &in[1]); sub_cpu(&in[2], &in[0], &in[1]);                                     \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2+=s1", 2, in, in + 1, add_assign(&in[1], &in[0]); add_assign_cpu(&in[1], &in[0]);                                         \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2-=s1", 2, in, in + 1, sub_assign(&in[1], &in[0]); sub_assign_cpu(&in[1], &in[0]);                                         \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2=-s1", 1, in, in + 1, minus(&in[1], &in[0]); minus_cpu(&in[1], &in[0]);                                                   \
                             , "GPU TEST", _prec);                                                                                                                \
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
                _TEST_GPU_OP(errors, "s2+=r*s1", 2, in, in + 1, double r = 10.0; mul_add_assign(&in[1], r, &in[0]);                                               \
                             mul_add_assign_cpu(&in[1], r, &in[0]);, "GPU TEST", _prec);                                                                          \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2=r*s1", 1, in, in + 1, double r = 10.0; mul(&in[1], r, &in[0]);                                                           \
                             mul_cpu(&in[1], r, &in[0]);, "GPU TEST", _prec);                                                                                     \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s3=s1+s2", 2, in, in + 2, add(&in[2], &in[0], &in[1]); add_cpu(&in[2], &in[0], &in[1]);                                     \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s3=s1-s2", 2, in, in + 2, sub(&in[2], &in[0], &in[1]); sub_cpu(&in[2], &in[0], &in[1]);                                     \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2+=s1", 2, in, in + 1, add_assign(&in[1], &in[0]); add_assign_cpu(&in[1], &in[0]);                                         \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2-=s1", 2, in, in + 1, sub_assign(&in[1], &in[0]); sub_assign_cpu(&in[1], &in[0]);                                         \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2=-s1", 1, in, in + 1, minus(&in[1], &in[0]); minus_cpu(&in[1], &in[0]);                                                   \
                             , "GPU TEST", _prec);                                                                                                                \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

#define _TEST_BASE(_FIELD_TYPE, _max_geom, _prec)                                                                                                                 \
    do {                                                                                                                                                          \
        for (int geom_idx = 0; geom_idx < _max_geom; ++geom_idx) {                                                                                                \
            _FIELD_TYPE *in = alloc(in, ninputs + 2, geometries[geom_idx]);                                                                                       \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST BASE", 0,                                                                                                                               \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s1=s2", 1, in, in + 1, copy(&in[1], &in[0]); copy_cpu(&in[1], &in[0]);                                                      \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s1=0", 0, in, in, zero(in); zero_cpu(in);, "GPU TEST", _prec);                                                              \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

#define _TEST_RED(_FIELD_TYPE, _max_geom, _prec)                                                                                                                  \
    do {                                                                                                                                                          \
        for (int geom_idx = 0; geom_idx < _max_geom; ++geom_idx) {                                                                                                \
            _FIELD_TYPE *in = alloc(in, ninputs + 2, geometries[geom_idx]);                                                                                       \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST RED", 0,                                                                                                                                \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "Re<s1,s2>", 2, in, double abs1 = prod_re(&in[0], &in[1]);                                                                   \
                             double abs2 = prod_re_cpu(&in[0], &in[1]);, "GPU TEST", _prec);                                                                      \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "Im<s1,s2>", 2, in, double abs1 = prod_im(&in[0], &in[1]);                                                                   \
                             double abs2 = prod_im_cpu(&in[0], &in[1]);, "GPU TEST", _prec);                                                                      \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "<s1,s2>", 2, in, hr_complex c1 = prod(&in[0], &in[1]);                                                                      \
                             hr_complex c2 = prod_cpu(&in[0], &in[1]); double abs1 = _complex_prod_re(c1, c1);                                                    \
                             double abs2 = _complex_prod_re(c2, c2);, "GPU TEST", _prec);                                                                         \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "|s1|^2", 1, in, double abs1 = sqnorm(&in[0]); double abs2 = sqnorm_cpu(&in[0]);                                             \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "|s1|_infty", 1, in, double abs1 = max(&in[0]); double abs2 = max_cpu(&in[0]);                                               \
                             , "GPU TEST", _prec);                                                                                                                \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

#define _TEST_RED_REAL(_FIELD_TYPE, _max_geom, _prec)                                                                                                             \
    do {                                                                                                                                                          \
        for (int geom_idx = 0; geom_idx < _max_geom; ++geom_idx) {                                                                                                \
            _FIELD_TYPE *in = alloc(in, ninputs + 2, geometries[geom_idx]);                                                                                       \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST RED", 0,                                                                                                                                \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "Re<s1,s2>", 2, in, double abs1 = prod_re(&in[0], &in[1]);                                                                   \
                             double abs2 = prod_re_cpu(&in[0], &in[1]);, "GPU TEST", _prec);                                                                      \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "<s1,s2>", 2, in, hr_complex c1 = prod(&in[0], &in[1]);                                                                      \
                             hr_complex c2 = prod_cpu(&in[0], &in[1]); double abs1 = _complex_prod_re(c1, c1);                                                    \
                             double abs2 = _complex_prod_re(c2, c2);, "GPU TEST", _prec);                                                                         \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "|s1|^2", 1, in, double abs1 = sqnorm(&in[0]); double abs2 = sqnorm_cpu(&in[0]);                                             \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "|s1|_infty", 1, in, double abs1 = max(&in[0]); double abs2 = max_cpu(&in[0]);                                               \
                             , "GPU TEST", _prec);                                                                                                                \
            }                                                                                                                                                     \
            free_field(in);                                                                                                                                       \
        }                                                                                                                                                         \
    } while (0)

#define _TEST_GAMMA_AND_LC(_FIELD_TYPE, _prec)                                                                                                                    \
    do {                                                                                                                                                          \
        for (int geom_idx = 0; geom_idx < 3; ++geom_idx) {                                                                                                        \
            _FIELD_TYPE *in = alloc(in, ninputs + 2, geometries[geom_idx]);                                                                                       \
                                                                                                                                                                  \
            for (int k = 0; k < niter; k++) {                                                                                                                     \
                lprintf(                                                                                                                                          \
                    "TEST GAMMA", 0,                                                                                                                              \
                    "Loop #%d geometry: %s, field: %s\n======================================================================================================\n", \
                    k, desc[geom_idx], _FIELD_DESC(in));                                                                                                          \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s3=r*s1+s*s2", 2, in, in + 2, double r = 2.0; double s = 3.0;                                                               \
                             lc(&in[2], r, &in[0], s, &in[1]); lc_cpu(&in[2], r, &in[0], s, &in[1]);, "GPU TEST", _prec);                                         \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s3+=r*s1+s*s2", 3, in, in + 2, double r = 2.0; double s = 3.0;                                                              \
                             lc_add_assign(&in[2], r, &in[0], s, &in[1]); lc_add_assign_cpu(&in[2], r, &in[0], s, &in[1]);                                        \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s3=c*s1+d*s2", 2, in, in + 2, hr_complex c = 2.0 + 3.0 * I;                                                                 \
                             hr_complex d = 3.0 + 4.0 * I; clc(&in[2], c, &in[0], d, &in[1]);                                                                     \
                             clc_cpu(&in[2], c, &in[0], d, &in[1]);, "GPU TEST", _prec);                                                                          \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s3+=c*s1+d*s2", 3, in, in + 2, hr_complex c = 2.0 + 3.0 * I;                                                                \
                             hr_complex d = 3.0 + 4.0 * I; clc_add_assign(&in[2], c, &in[0], d, &in[1]);                                                          \
                             clc_add_assign_cpu(&in[2], c, &in[0], d, &in[1]);, "GPU TEST", _prec);                                                               \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s2=g5*s1", 1, in, in + 1, g5(&in[1], &in[0]); g5_cpu(&in[1], &in[0]);                                                       \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s1=g5*s1", 1, in, in, g5_assign(&in[0]); g5_assign_cpu(&in[0]);, "GPU TEST", _prec);                                        \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s1=g0*s2", 1, in, in + 1, g0(&in[1], &in[0]); g0_cpu(&in[1], &in[0]);                                                       \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s1=g1*s2", 1, in, in + 1, g1(&in[1], &in[0]); g1_cpu(&in[1], &in[0]);                                                       \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s1=g2*s2", 1, in, in + 1, g2(&in[1], &in[0]); g2_cpu(&in[1], &in[0]);                                                       \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_GPU_OP(errors, "s1=g3*s2", 1, in, in + 1, g3(&in[1], &in[0]); g3_cpu(&in[1], &in[0]);                                                       \
                             , "GPU TEST", _prec);                                                                                                                \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "Im<g5*s1,s2>", 2, in, double abs1 = g5_prod_im(&in[0], &in[1]);                                                             \
                             double abs2 = g5_prod_im_cpu(&in[0], &in[1]);, "GPU TEST", _prec);                                                                   \
                                                                                                                                                                  \
                _TEST_RED_OP(errors, "Re<g5*s1,s2>", 2, in, double abs1 = g5_prod_re(&in[0], &in[1]);                                                             \
                             double abs2 = g5_prod_re_cpu(&in[0], &in[1]);, "GPU TEST", _prec);                                                                   \
                free_field(in);                                                                                                                                   \
            }                                                                                                                                                     \
        }                                                                                                                                                         \
    } while (0)

#define _TEST

static geometry_descriptor *geometries[3] = { &glattice, &glat_even, &glat_odd };
static char *desc[3] = { "glattice", "glat_even", "glat_odd" };

int main(int argc, char *argv[]) {
    /* setup process id and communications */
    //logger_setlevel(0,10000);

    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    rlxs_init(1, 1234);

    const int niter = 1;
    // Allocate memory for CPU and GPU spinor fields
    // add 2 for the output results used in the macro TEST
    int ninputs = 3; //max number of inputs

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
    _TEST_RED_REAL(scalar_field, 1, EPSILON_TEST);
    _TEST_RED_REAL(suNg_field, 1, EPSILON_TEST);
    _TEST_RED_REAL(suNg_field_flt, 1, EPSILON_FLT_TEST);
#ifdef REPR_IS_REAL
    _TEST_RED_REAL(suNfc_field, 1, EPSILON_TEST);
#endif
    _TEST_RED_REAL(suNf_field, 1, EPSILON_TEST);
    _TEST_RED_REAL(suNf_field_flt, 1, EPSILON_FLT_TEST);
    _TEST_RED_REAL(suNg_av_field, 1, EPSILON_TEST);
    _TEST_RED_REAL(gtransf, 1, EPSILON_TEST);
    _TEST_RED_REAL(clover_term, 1, EPSILON_TEST);
    _TEST_RED_REAL(clover_force, 1, EPSILON_TEST);
    _TEST_RED_REAL(staple_field, 1, EPSILON_TEST);

    /* Gamma and linear combinations */
    _TEST_GAMMA_AND_LC(spinor_field, EPSILON_TEST);
    _TEST_GAMMA_AND_LC(spinor_field_flt, EPSILON_FLT_TEST);

    finalize_process();

    return errors;
}
