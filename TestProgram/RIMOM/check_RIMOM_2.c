/*******************************************************************************
 * NOCOMPILE= BASIC_SF
 * NOCOMPILE= BC_T_SF
 * NOCOMPILE= ROTATED_SF
 * NOCOMPILE= BC_T_THETA
 * NOCOMPILE= FERMION_THETA && !BC_X_THETA
 * NOCOMPILE= FERMION_THETA && !BC_Y_THETA
 * NOCOMPILE= FERMION_THETA && !BC_Z_THETA
 *
 * RI'-MOM self-test 2: gauge covariance (CHECK_GT).
 *
 * On a random configuration, U -> g U g^dag together with eta -> g eta must
 * leave every vertex channel unchanged (Sin/Sout excluded: the propagator is
 * gauge variant). This covers the conserved current cg0..cg3 in all four
 * directions, the T direction included. No twist.
 *
 * Run with MPI (check_RIMOM_2.in: NP_X = 2) so that the gauge transformation
 * and the conserved current cross rank boundaries.
 *
 * Test by Vincent Drach
 *******************************************************************************/

#define MAIN_PROGRAM

#include "check_RIMOM_utils.h"

static const double TOL_GAUGE_COVARIANCE = 1e-8;

static int check_gauge_covariance(rimom_ctx *ctx, double p[4], double tol) {
    /* random links: not gauge-equivalent to U=1, so U != U^dag matters */
    random_u(u_gauge);
    represent_gauge_field(); /* also completes random_u's halo exchange */
    double plaq_before = avr_plaquette();
    lprintf("CHECK_GT", 0, "random configuration, plaquette %1.6f\n", plaq_before);

    data_storage_array *outA =
        measure_all_channels(ctx->source, ctx->prop_in, ctx->prop_out, ctx->nm, ctx->m, p, 0, "GT check A");

    random_gauge_transform(u_gauge); /* U -> g U g^dag, g left in gaugefix.c's global g */
    represent_gauge_field();
    lprintf("CHECK_GT", 0, "after random gauge transformation, plaquette %1.6f (must be unchanged)\n", avr_plaquette());
    if (fabs(avr_plaquette() - plaq_before) > 1e-12) {
        lprintf("CHECK_GT", 0, "plaquette changed: the gauge transformation is broken FAILED\n");
        free_data_storage(outA);
        return 0;
    }

    data_storage_array *outB =
        measure_all_channels(ctx->source, ctx->prop_in, ctx->prop_out, ctx->nm, ctx->m, p, 1, "GT check B");

    int ok = compare_channels(outA, outB, vertex_channels, N_VERTEX_CHANNELS, tol, "CHECK_GT");
    free_data_storage(outB);
    free_data_storage(outA);
    return ok;
}

int main(int argc, char *argv[]) {
    rimom_ctx ctx;
    rimom_setup(&argc, &argv, "check_RIMOM_2", &ctx);
    double p[4] = { mes_var.pt_in, mes_var.px_in, mes_var.py_in, mes_var.pz_in };
    int return_value = check_gauge_covariance(&ctx, p, TOL_GAUGE_COVARIANCE) ? 0 : 1;
    return rimom_finish(&ctx, return_value, "check_RIMOM_2", "gauge covariance of all vertex channels");
}
