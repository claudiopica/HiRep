/*******************************************************************************
 * NOCOMPILE= BASIC_SF
 * NOCOMPILE= BC_T_SF
 * NOCOMPILE= ROTATED_SF
 * NOCOMPILE= BC_T_THETA
 * NOCOMPILE= FERMION_THETA && !BC_X_THETA
 * NOCOMPILE= FERMION_THETA && !BC_Y_THETA
 * NOCOMPILE= FERMION_THETA && !BC_Z_THETA
 *
 * RI'-MOM self-test 3: twist equivalence (CHECK_TW).
 *
 * On a random configuration, momentum n with twist theta = +n and momentum 2n
 * with theta = -n reach the same p + theta/2 = 3n/2 (units of 2*pi/L). The
 * two setups differ by the U(1) gauge transformation exp(2*pi*i n.x/L),
 * which is periodic, so every channel, Sin/Sout included, must be identical.
 * Tests the twist (links + MPI buffers, or eitheta, as in the build: see
 * check_RIMOM_utils.h) and the Fourier phases of the vertices, the conserved
 * current included (twist on the point-split link). Only twisted directions
 * are tested: n must have no T component, and every spatial component should
 * be non-zero (check_RIMOM_3.in: n = (0,1,1,-1)) so that cg1..cg3 are all
 * covered. cg0 is covered by check_RIMOM_2.
 *
 * Run with MPI (check_RIMOM_3.in: NP_X = 2).
 *
 * Test by Vincent Drach
 *******************************************************************************/

#define MAIN_PROGRAM

#include "check_RIMOM_utils.h"

static const double TOL_TWIST = 1e-8;

static int check_twist_equivalence(rimom_ctx *ctx, double tol) {
    const double n[4] = { mes_var.pt_in, mes_var.px_in, mes_var.py_in, mes_var.pz_in };
    if (n[0] != 0. || (n[1] == 0. && n[2] == 0. && n[3] == 0.)) {
        lprintf("CHECK_TW", 0, "FAILED: needs a momentum with no T component and a non-zero spatial one\n");
        return 0;
    }
    if (n[1] == 0. || n[2] == 0. || n[3] == 0.) {
        lprintf("CHECK_TW", 0, "WARNING: a spatial component of n is zero, that direction of cg is not tested\n");
    }

    random_u(u_gauge);
    represent_gauge_field(); /* also completes random_u's halo exchange */
    lprintf("CHECK_TW", 0, "random configuration, plaquette %1.6f\n", avr_plaquette());
    suNf_field *u_untwisted = alloc_suNf_field(&glattice);
    copy_suNf_field(u_untwisted, u_gauge_f);

    double pA[4], pB[4];
    for (int mu = 0; mu < 4; mu++) {
        pA[mu] = n[mu];
        pB[mu] = 2. * n[mu];
    }

    apply_twist(u_untwisted, +n[1], +n[2], +n[3]);
    data_storage_array *outA = measure_all_channels(ctx->source, ctx->prop_in, ctx->prop_out, ctx->nm, ctx->m, pA, 0,
                                                    "TW check A (p=n, theta=+n)");

    apply_twist(u_untwisted, -n[1], -n[2], -n[3]);
    data_storage_array *outB = measure_all_channels(ctx->source, ctx->prop_in, ctx->prop_out, ctx->nm, ctx->m, pB, 0,
                                                    "TW check B (p=2n, theta=-n)");

    static const char *propagators[2] = { "Sin", "Sout" };
    int ok = compare_channels(outA, outB, propagators, 2, tol, "CHECK_TW");
    if (!compare_channels(outA, outB, vertex_channels, N_VERTEX_CHANNELS, tol, "CHECK_TW")) { ok = 0; }

    free_data_storage(outB);
    free_data_storage(outA);
    apply_twist(u_untwisted, 0., 0., 0.);
    free_suNf_field(u_untwisted);
    return ok;
}

int main(int argc, char *argv[]) {
    rimom_ctx ctx;
    rimom_setup(&argc, &argv, "check_RIMOM_3", &ctx);
    int return_value = check_twist_equivalence(&ctx, TOL_TWIST) ? 0 : 1;
    return rimom_finish(&ctx, return_value, "check_RIMOM_3", "twist equivalence of all channels");
}
