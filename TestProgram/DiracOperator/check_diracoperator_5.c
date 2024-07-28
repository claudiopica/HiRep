/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that CPU and GPU copies are identical
*
********************************************************************************/

#include "libhr.h"

static int errors = 0; // count the number of errors during this test uni
static double hmass = 0.1;

#define EPSILON_SQ_TEST 1.e-12
#define EPSILON_SQ_FLT_TEST 1.e-3

spinor_field *tmp, *in;
spinor_field in_oe[3], in_eo[3], in_ee[3], in_oo[3];

#ifdef DPHI_FLT
spinor_field_flt *tmp_flt, *in_flt;
spinor_field_flt in_eo_flt[3], in_oe_flt[3], in_ee_flt[3], in_oo_flt[3];
#endif

void reset_fields() {
    in = tmp;
#ifdef DPHI_FLT
    in_flt = tmp_flt;
#endif
}

void setup_fields_oe() {
    in_oe[0] = get_even_part_spinor_field(&in[0]);
    in_oe[1] = get_odd_part_spinor_field(&in[1]); // out
    in_oe[2] = get_odd_part_spinor_field(&in[2]); // diff

    tmp = in;
    in = in_oe;
#ifdef DPHI_FLT
    in_oe_flt[0] = get_even_part_spinor_field_flt(&in_flt[0]);
    in_oe_flt[1] = get_odd_part_spinor_field_flt(&in_flt[1]); // out
    in_oe_flt[2] = get_odd_part_spinor_field_flt(&in_flt[2]); // diff

    tmp_flt = in_flt;
    in_flt = in_oe_flt;
#endif
}

void setup_fields_eo() {
    in_eo[0] = get_odd_part_spinor_field(&in[0]);
    in_eo[1] = get_even_part_spinor_field(&in[1]); // out
    in_eo[2] = get_even_part_spinor_field(&in[2]); // diff

    tmp = in;
    in = in_eo;

#ifdef DPHI_FLT
    in_eo_flt[0] = get_odd_part_spinor_field_flt(&in_flt[0]);
    in_eo_flt[1] = get_even_part_spinor_field_flt(&in_flt[1]); // out
    in_eo_flt[2] = get_even_part_spinor_field_flt(&in_flt[2]); // diff

    tmp_flt = in_flt;
    in_flt = in_eo_flt;
#endif
}

void setup_fields_ee() {
    in_ee[0] = get_even_part_spinor_field(&in[0]);
    in_ee[1] = get_even_part_spinor_field(&in[1]); // out
    in_ee[2] = get_even_part_spinor_field(&in[2]); // diff

    tmp = in;
    in = in_ee;

#ifdef DPHI_FLT
    in_ee_flt[0] = get_even_part_spinor_field_flt(&in_flt[0]);
    in_ee_flt[1] = get_even_part_spinor_field_flt(&in_flt[1]); // out
    in_ee_flt[2] = get_even_part_spinor_field_flt(&in_flt[2]); // diff

    tmp_flt = in_flt;
    in_flt = in_ee_flt;
#endif
}

void setup_fields_oo() {
    in_oo[0] = get_odd_part_spinor_field(&in[0]);
    in_oo[1] = get_odd_part_spinor_field(&in[1]); // out
    in_oo[2] = get_odd_part_spinor_field(&in[2]); // diff

    tmp = in;
    in = in_oo;

#ifdef DPHI_FLT
    in_oo_flt[0] = get_odd_part_spinor_field_flt(&in_flt[0]);
    in_oo_flt[1] = get_odd_part_spinor_field_flt(&in_flt[1]); // out
    in_oo_flt[2] = get_odd_part_spinor_field_flt(&in_flt[2]); // diff

    tmp_flt = in_flt;
    in_flt = in_oo_flt;
#endif
}

int main(int argc, char *argv[]) {
    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary
    setup_process(&argc, &argv);
    int ninputs = 1;
    int noutputs = 1;

    in = alloc_spinor_field(ninputs + noutputs + 1, &glattice);

#ifdef DPHI_FLT
    in_flt = alloc_spinor_field_flt(ninputs + noutputs + 1, &glattice);
#endif

    setup_random_gauge_fields();

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    compute_clover_term();
    double csw_check = get_csw(); // Query GPU setting of csw
    set_csw_cpu(&csw_check); // Set CPU setting to the same value
    setup_clover();
#endif

    _TEST_GPU_OP(errors, "Unit", ninputs + noutputs + 1, in, in + 1, mul_spinor_field(&in[1], 1, in);
                 mul_spinor_field_cpu(&in[1], 1, in);, "GPU TEST", 0.0);

    _TEST_GPU_OP(errors, "Dphi", ninputs + noutputs + 1, in, in + 1, Dphi(-hmass, &in[1], in); Dphi_cpu(-hmass, &in[1], in);
                 , "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Dphi_", ninputs + noutputs + 1, in, in + 1, Dphi_(&in[1], in); Dphi_cpu_(&in[1], in);
                 , "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "g5Dphi", ninputs + noutputs + 1, in, in + 1, g5Dphi(-hmass, &in[1], in);
                 g5Dphi_cpu(-hmass, &in[1], in);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Q^2", ninputs + noutputs + 1, in, in + 1, g5Dphi_sq(-hmass, &in[1], in);
                 g5Dphi_sq_cpu(-hmass, &in[1], in);, "GPU_TEST", EPSILON_SQ_TEST);

#ifdef DPHI_FLT

    _TEST_GPU_OP(errors, "Dphi_flt", ninputs + noutputs + 1, in_flt, in_flt + 1, Dphi_flt(-hmass, &in_flt[1], in_flt);
                 Dphi_flt_cpu(-hmass, &in_flt[1], in_flt);, "GPU_TEST", EPSILON_FLT_TEST);

    _TEST_GPU_OP(errors, "Dphi_flt_", ninputs + noutputs + 1, in_flt, in_flt + 1, Dphi_flt_(&in_flt[1], in_flt);
                 Dphi_flt_cpu_(&in_flt[1], in_flt);, "GPU_TEST", EPSILON_FLT_TEST);

    _TEST_GPU_OP(errors, "g5Dphi_flt", ninputs + noutputs + 1, in_flt, in_flt + 1, g5Dphi_flt(-hmass, &in_flt[1], in_flt);
                 g5Dphi_flt_cpu(-hmass, &in_flt[1], in_flt);, "GPU_TEST", EPSILON_FLT_TEST);

    _TEST_GPU_OP(errors, "Q^2 flt", ninputs + noutputs + 1, in_flt, in_flt + 1, g5Dphi_sq_flt(-hmass, &in_flt[1], in_flt);
                 g5Dphi_sq_flt_cpu(-hmass, &in_flt[1], in_flt);, "GPU_TEST", EPSILON_SQ_FLT_TEST);
#endif

#ifdef WITH_CLOVER

    _TEST_GPU_OP(errors, "Cphi_", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_ (+=)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 1);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_", ninputs + noutputs + 1, in, in + 1, Cphi_inv_(-hmass, &in[1], in, 0);
                 Cphi_inv_cpu_(-hmass, &in[1], in, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_ (+=)", ninputs + noutputs + 1, in, in + 1, Cphi_inv_(-hmass, &in[1], in, 1);
                 Cphi_inv_cpu_(-hmass, &in[1], in, 1);, "GPU_TEST", EPSILON_TEST);

#endif

#ifdef WITH_EXPCLOVER

    _TEST_GPU_OP(errors, "Cphi_", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 0, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_ (+=)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 1, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 0, 1);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_ (+=)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 1, 1);, "GPU_TEST", EPSILON_TEST);

#endif

    setup_fields_oe();

    _TEST_GPU_OP(errors, "Dphi_ (OE)", ninputs + noutputs + 1, in, in + 1, Dphi_(&in[1], in); Dphi_cpu_(&in[1], in);
                 , "GPU_TEST", EPSILON_TEST);

#ifdef DPHI_FLT

    _TEST_GPU_OP(errors, "Dphi_flt_ (OE)", ninputs + noutputs + 1, in_flt, in_flt + 1, Dphi_flt_(&in_flt[1], in_flt);
                 Dphi_flt_cpu_(&in_flt[1], in_flt);, "GPU_TEST", EPSILON_FLT_TEST);

#endif

    reset_fields();
    setup_fields_eo();

    _TEST_GPU_OP(errors, "Dphi_ (EO)", ninputs + noutputs + 1, in, in + 1, Dphi_(&in[1], in); Dphi_cpu_(&in[1], in);
                 , "GPU_TEST", EPSILON_TEST);

#ifdef DPHI_FLT

    _TEST_GPU_OP(errors, "Dphi_flt_ (EO)", ninputs + noutputs + 1, in_flt, in_flt + 1, Dphi_flt_(&in_flt[1], in_flt);
                 Dphi_flt_cpu_(&in_flt[1], in_flt);, "GPU_TEST", EPSILON_TEST);

#endif

    reset_fields();
    setup_fields_ee();

    _TEST_GPU_OP(errors, "Dphi_eopre", ninputs + noutputs + 1, in, in + 1, Dphi_eopre(-hmass, &in[1], in);
                 Dphi_eopre_cpu(-hmass, &in[1], in);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "g5Dphi_oepre", ninputs + noutputs + 1, in, in + 1, g5Dphi_eopre(-hmass, &in[1], in);
                 g5Dphi_eopre_cpu(-hmass, &in[1], in);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "g5Dphi_oepre_sq", ninputs + noutputs + 1, in, in + 1, g5Dphi_eopre_sq(-hmass, &in[1], in);
                 g5Dphi_eopre_sq_cpu(-hmass, &in[1], in);, "GPU_TEST", EPSILON_SQ_TEST);

    _TEST_GPU_OP(errors, "Q_eopre", ninputs + noutputs + 1, in, in + 1, g5Dphi_eopre(-hmass, &in[1], in);
                 g5Dphi_eopre_cpu(-hmass, &in[1], in);, "GPU_TEST", EPSILON_TEST);

#ifdef DPHI_FLT

    _TEST_GPU_OP(errors, "Dphi_eopre_flt", ninputs + noutputs + 1, in_flt, in_flt + 1,
                 Dphi_eopre_flt(-hmass, &in_flt[1], in_flt);
                 Dphi_eopre_flt_cpu(-hmass, &in_flt[1], in_flt);, "GPU_TEST", EPSILON_FLT_TEST);

    _TEST_GPU_OP(errors, "Q_eopre_flt", ninputs + noutputs + 1, in_flt, in_flt + 1,
                 g5Dphi_eopre_flt(-hmass, &in_flt[1], in_flt);
                 g5Dphi_eopre_flt_cpu(-hmass, &in_flt[1], in_flt);, "GPU_TEST", EPSILON_FLT_TEST);

#endif

#ifdef WITH_CLOVER

    _TEST_GPU_OP(errors, "Cphi_ (EE)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_ (EE,+=)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 1);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_ (EE)", ninputs + noutputs + 1, in, in + 1, Cphi_inv_(-hmass, &in[1], in, 0);
                 Cphi_inv_cpu_(-hmass, &in[1], in, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_,EE,+=", ninputs + noutputs + 1, in, in + 1, Cphi_inv_(-hmass, &in[1], in, 1);
                 Cphi_inv_cpu_(-hmass, &in[1], in, 1);, "GPU_TEST", EPSILON_TEST);

#endif

#ifdef WITH_EXPCLOVER

    _TEST_GPU_OP(errors, "Cphi_ (EE)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 0, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_ (EE,+=)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 1, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_ (EE)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 0, 1);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_,EE,+=", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 1, 1);, "GPU_TEST", EPSILON_TEST);

#endif

    reset_fields();
    setup_fields_oo();

    _TEST_GPU_OP(errors, "Dphi_oepre", ninputs + noutputs + 1, in, in + 1, Dphi_oepre(-hmass, &in[1], in);
                 Dphi_oepre_cpu(-hmass, &in[1], in);, "GPU_TEST", EPSILON_TEST);

#ifdef WITH_CLOVER

    _TEST_GPU_OP(errors, "Cphi_ (OO)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_ (OO,+=)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 1);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_ (OO)", ninputs + noutputs + 1, in, in + 1, Cphi_inv_(-hmass, &in[1], in, 0);
                 Cphi_inv_cpu_(-hmass, &in[1], in, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_,OO,+=", ninputs + noutputs + 1, in, in + 1, Cphi_inv_(-hmass, &in[1], in, 1);
                 Cphi_inv_cpu_(-hmass, &in[1], in, 1);, "GPU_TEST", EPSILON_TEST);

#endif

#ifdef WITH_EXPCLOVER

    _TEST_GPU_OP(errors, "Cphi_ (OO)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 0, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_ (OO,+=)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1, 0);
                 Cphi_cpu_(-hmass, &in[1], in, 1, 0);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_ (OO)", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 0, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 0, 1);, "GPU_TEST", EPSILON_TEST);

    _TEST_GPU_OP(errors, "Cphi_inv_,OO,+=", ninputs + noutputs + 1, in, in + 1, Cphi_(-hmass, &in[1], in, 1, 1);
                 Cphi_cpu_(-hmass, &in[1], in, 1, 1);, "GPU_TEST", EPSILON_TEST);

#endif

    reset_fields();
    finalize_process();
    return errors;
}