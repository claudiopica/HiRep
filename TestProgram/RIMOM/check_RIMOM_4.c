/*******************************************************************************
 * NOCOMPILE= !BC_X_THETA
 * NOCOMPILE= !BC_Y_THETA
 * NOCOMPILE= !BC_Z_THETA
 * NOCOMPILE= WITH_GPU
 * NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
 *
 * RI'-MOM self-test 4: twisting the links vs eitheta (theta builds only).
 *
 * Old method (former measure_Z_mom.c; non-theta path of check_RIMOM_1..3):
 *   multiply every spatial link
 *   of u_gauge_f by exp(i pi theta_mu / L_mu)            [theta in units of pi]
 * New method (FERMION_THETA): leave the links alone, set
 *   eitheta[mu] = exp(i pi theta_mu / L_mu), which Dphi applies to each
 *   forward hop (and its conjugate to each backward hop).
 *
 * Both give the same Dirac operator, so on the SAME random gauge
 * configuration, with the clover term on, the propagators from a momentum
 * source must agree to solver precision. In the old-method run the clover
 * term is recomputed from the twisted links (closed loops, so the phases
 * should cancel): agreement therefore also checks that the clover term is
 * insensitive to the twist.
 *
 * Extra check: for theta != 0 the new-method propagator must differ from the
 * untwisted one, i.e. eitheta is really used by the operator.
 *
 * Gauge-field integrity: after every propagator computation the plaquette,
 * the saved untwisted links and the links in use are checked (see
 * check_gauge). This looks for the gauge-field corruption seen in
 * measure_Z_mom with expclover + MPI.
 *
 * Gauge field: random (default) or read from mes:configfile, optionally
 * Landau gauge fixed (mes:gaugefix = 1) as in measure_Z_mom.c. The lattice
 * size in the input file must match the configuration.
 *
 * MPI: copy_suNf_field only copies local sites, so after every restore of
 * the untwisted links the boundary buffers are refreshed by a sendrecv
 * (restore_links). Field comparisons are done on local sites and summed
 * over processes.
 *
 * Test by Vincent Drach
 *******************************************************************************/

#include "libhr.h"
#include <string.h>

#ifndef FERMION_THETA
#error check_RIMOM_4 needs FERMION_THETA (compile with BC_X_THETA, BC_Y_THETA, BC_Z_THETA)
#endif

/* Same function as in measure_Z_mom.c / check_RIMOM.c (old method). */
static void twist_XYZ_bc(double theta_x, double theta_y, double theta_z) {
    int index;
    int ix, iy, iz, it;
    suNf *u;
    suNf utmp;
    hr_complex eith_x, eith_y, eith_z;
    eith_x = cexp(I * PI * theta_x / (double)GLB_X);
    eith_y = cexp(I * PI * theta_y / (double)GLB_Y);
    eith_z = cexp(I * PI * theta_z / (double)GLB_Z);

    for (it = 0; it < T_EXT; ++it) {
        for (ix = 0; ix < X_EXT; ++ix) {
            for (iy = 0; iy < Y_EXT; ++iy) {
                for (iz = 0; iz < Z_EXT; ++iz) {
                    index = ipt_ext(it, ix, iy, iz);
                    u = pu_gauge_f(index, 1);
                    _suNf_mulc(utmp, eith_x, *u);
                    *u = utmp;
                    u = pu_gauge_f(index, 2);
                    _suNf_mulc(utmp, eith_y, *u);
                    *u = utmp;
                    u = pu_gauge_f(index, 3);
                    _suNf_mulc(utmp, eith_z, *u);
                    *u = utmp;
                }
            }
        }
    }
}

/* New method: same phases, put into eitheta (index 0 = T, 1..3 = X, Y, Z). */
static void set_eitheta(double theta_x, double theta_y, double theta_z) {
    eitheta[0] = 1.0;
    eitheta[1] = cexp(I * PI * theta_x / (double)GLB_X);
    eitheta[2] = cexp(I * PI * theta_y / (double)GLB_Y);
    eitheta[3] = cexp(I * PI * theta_z / (double)GLB_Z);
}

typedef struct _input_twist {
    char configfile[256]; /* empty: random gauge field */
    int gaugefix;         /* 1: Landau gauge fix the configuration first */
    double mass;
    double precision;
    double csw;
    double tol;
    int n_mom;
    int n_twist;
    int pt;
    int px;
    int py;
    int pz;

    input_record_t read[13];
} input_twist;

#define init_input_twist(varname)                                                         \
    {                                                                                     \
        .read = {                                                                         \
            { "gauge configuration (optional)", "mes:configfile = %s", STRING_T, (varname).configfile }, \
            { "gauge fix the configuration (optional)", "mes:gaugefix = %d", INT_T, &(varname).gaugefix }, \
            { "quark mass", "mes:mass = %lf", DOUBLE_T, &(varname).mass },                \
            { "inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision }, \
            { "csw coefficient", "mes:csw = %lf", DOUBLE_T, &(varname).csw },             \
            { "relative tolerance", "mes:tol = %lf", DOUBLE_T, &(varname).tol },          \
            { "number of momenta", "mes:n_mom = %d", INT_T, &(varname).n_mom },           \
            { "number of twists", "mes:n_twist = %d", INT_T, &(varname).n_twist },        \
            { "momentum t", "mes:pt = %d", INT_T, &(varname).pt },                        \
            { "momentum x", "mes:px = %d", INT_T, &(varname).px },                        \
            { "momentum y", "mes:py = %d", INT_T, &(varname).py },                        \
            { "momentum z", "mes:pz = %d", INT_T, &(varname).pz },                        \
            { NULL, NULL, INT_T, NULL }                                                   \
        }                                                                                 \
    }

input_twist mes_var = init_input_twist(mes_var);

/* Propagator (4 spins x NF colours) from a momentum source, with whatever
 * links / eitheta are currently set. */
static void momentum_propagator(spinor_field *prop, spinor_field *source, int pt, int px, int py, int pz) {
    for (int k = 0; k < NF; ++k) {
        create_gauge_fixed_momentum_source(source, pt, px, py, pz, k);
        calc_propagator(prop + 4 * k, source, 4);
    }
}

/* Exact (bitwise) comparison of two gauge fields on the local sites (the
 * boundary buffers are not copied by copy_suNf_field), summed over all
 * processes. Returns 1 if identical everywhere. */
static int same_suNf_field(suNf_field *a, suNf_field *b) {
    int ndiff = 0;
    for (int it = 0; it < T; ++it) {
        for (int ix = 0; ix < X; ++ix) {
            for (int iy = 0; iy < Y; ++iy) {
                for (int iz = 0; iz < Z; ++iz) {
                    int idx = ipt(it, ix, iy, iz);
                    if (memcmp(_4FIELD_AT(a, idx, 0), _4FIELD_AT(b, idx, 0), 4 * sizeof(suNf)) != 0) { ndiff++; }
                }
            }
        }
    }
    global_sum_int(&ndiff, 1);
    return ndiff == 0;
}

/* Put the untwisted links back into u_gauge_f, including the boundary
 * buffers (needed with MPI, harmless without), and mark the clover term
 * for recomputation. */
static void restore_links(suNf_field *u_save) {
    copy_suNf_field(u_gauge_f, u_save);
#ifdef WITH_MPI
    start_sendrecv_suNf_field(u_gauge_f);
    complete_sendrecv_suNf_field(u_gauge_f);
#endif
    stale_clover_cpu = 1;
}

/* Gauge-field integrity check, run after every propagator computation:
 *  - the plaquette must still be plaq0 (twisting does not change it);
 *  - the saved untwisted copy u_save must still equal the pristine copy;
 *  - the links currently in use must still equal the ones we set up.
 * Returns 1 if everything is intact. */
static double plaq0;
static int check_gauge(const char *label, suNf_field *u_save, suNf_field *u_pristine, suNf_field *expected) {
    double plaq = avr_plaquette();
    int plaq_ok = fabs(plaq - plaq0) < 1e-10;
    int save_ok = same_suNf_field(u_save, u_pristine);
    int links_ok = same_suNf_field(u_gauge_f, expected);
    if (!plaq_ok || !save_ok || !links_ok) {
        lprintf("GAUGE_CHECK", 0, "%s: CORRUPTED  <p> = %.12f (expected %.12f)%s%s\n", label, plaq, plaq0,
                save_ok ? "" : ", saved untwisted copy changed", links_ok ? "" : ", links in use changed");
        return 0;
    }
    lprintf("GAUGE_CHECK", 10, "%s: intact, <p> = %.12f\n", label, plaq);
    return 1;
}

/* |a - b|^2 / |a|^2 over all 4*NF spinor fields; tmp is overwritten. */
static double rel_diff2(spinor_field *a, spinor_field *b, spinor_field *tmp) {
    double num = 0., den = 0.;
    for (int j = 0; j < 4 * NF; ++j) {
        copy_spinor_field(tmp, &a[j]);
        sub_assign_spinor_field(tmp, &b[j]);
        num += sqnorm_spinor_field(tmp);
        den += sqnorm_spinor_field(&a[j]);
    }
    return num / den;
}

int main(int argc, char *argv[]) {
    int return_value = 0;
    double m[1];

    setup_process(&argc, &argv);
    setup_random_gauge_fields();

    read_input(mes_var.read, get_input_filename());

    lprintf("MAIN", 0, "Compiled with macros: %s\n", MACROS);
    lprintf("MAIN", 0, "Gauge group: SU(%d), representation " REPR_NAME " [dim=%d]\n", NG, NF);
    lprintf("MAIN", 0, "mass = %f, precision = %e, tolerance = %e\n", mes_var.mass, mes_var.precision, mes_var.tol);
    lprintf("MAIN", 0, "base momentum (%d, %d, %d, %d), n_mom = %d, n_twist = %d\n", mes_var.pt, mes_var.px, mes_var.py,
            mes_var.pz, mes_var.n_mom, mes_var.n_twist);

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    set_csw(&mes_var.csw);
    lprintf("MAIN", 0, "csw = %f\n", mes_var.csw);
#else
    lprintf("MAIN", 0, "No clover term (plain Wilson): csw is ignored\n");
#endif

    if (mes_var.configfile[0] != '\0') {
        lprintf("MAIN", 0, "Reading gauge configuration %s\n", mes_var.configfile);
        read_gauge_field(mes_var.configfile);
        lprintf("MAIN", 0, "configuration read, <p> = %1.6f\n", avr_plaquette());
        if (mes_var.gaugefix) {
            double act = gaugefix(10, 1.8, 10000, 1e-10, u_gauge); /* Landau, as in measure_Z_mom.c */
            lprintf("MAIN", 0, "Landau gauge fixed, action %1.6f, <p> = %1.6f\n", act, avr_plaquette());
        }
        represent_gauge_field();
#ifdef WITH_MPI
        start_sendrecv_suNf_field(u_gauge_f);
        complete_sendrecv_suNf_field(u_gauge_f);
#endif
    } else {
        lprintf("MAIN", 0, "random gauge field, <p> = %1.6f\n", avr_plaquette());
    }

    /* Untwisted links, restored before every run.
     * NB: for REPR_FUNDAMENTAL u_gauge_f is the same memory as u_gauge. */
    suNf_field *u_save = alloc_suNf_field(&glattice);
    copy_suNf_field(u_save, u_gauge_f);
    /* second untwisted copy, never used for anything: detects writes into u_save */
    suNf_field *u_pristine = alloc_suNf_field(&glattice);
    copy_suNf_field(u_pristine, u_gauge_f);
    /* twisted links of the current old-method run, to detect writes into u_gauge_f */
    suNf_field *u_twisted = alloc_suNf_field(&glattice);
    plaq0 = avr_plaquette();
    char label[256];

    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop_new = alloc_spinor_field(4 * NF, &glattice);
    spinor_field *prop_old = alloc_spinor_field(4 * NF, &glattice);
    spinor_field *prop_ref = alloc_spinor_field(4 * NF, &glattice);
    spinor_field *tmp = alloc_spinor_field(1, &glattice);

    m[0] = mes_var.mass;
    init_propagator_eo(1, m, mes_var.precision);

    for (int l = 1; l <= mes_var.n_mom; ++l) {
        int p[4] = { mes_var.pt * l, mes_var.px * l, mes_var.py * l, mes_var.pz * l };

        /* untwisted reference: untwisted links, eitheta = 1 */
        restore_links(u_save);
        set_eitheta(0., 0., 0.);
        momentum_propagator(prop_ref, source, p[0], p[1], p[2], p[3]);
        sprintf(label, "p=(%d,%d,%d,%d) untwisted reference", p[0], p[1], p[2], p[3]);
        if (!check_gauge(label, u_save, u_pristine, u_save)) { return_value = 1; }

        for (int tw = -mes_var.n_twist; tw <= mes_var.n_twist; ++tw) {
            double twist = (double)tw * 0.5; /* same convention as measure_Z_mom.c */
            double th[3] = { twist * p[1], twist * p[2], twist * p[3] };

            /* New method: untwisted links, phases in eitheta */
            restore_links(u_save);
            set_eitheta(th[0], th[1], th[2]);
            momentum_propagator(prop_new, source, p[0], p[1], p[2], p[3]);
            sprintf(label, "p=(%d,%d,%d,%d) tw=%+d new method", p[0], p[1], p[2], p[3], tw);
            if (!check_gauge(label, u_save, u_pristine, u_save)) { return_value = 1; }

            /* Old method: phases in the links, eitheta = 1, clover term
             * recomputed from the twisted links */
            restore_links(u_save);
            set_eitheta(0., 0., 0.);
            twist_XYZ_bc(th[0], th[1], th[2]); /* local sites and buffers, as in measure_Z_mom.c */
            copy_suNf_field(u_twisted, u_gauge_f);
            stale_clover_cpu = 1;
            momentum_propagator(prop_old, source, p[0], p[1], p[2], p[3]);
            sprintf(label, "p=(%d,%d,%d,%d) tw=%+d old method", p[0], p[1], p[2], p[3], tw);
            if (!check_gauge(label, u_save, u_pristine, u_twisted)) { return_value = 1; }

            double d_old_new = sqrt(rel_diff2(prop_new, prop_old, tmp));
            double d_new_ref = sqrt(rel_diff2(prop_new, prop_ref, tmp));

            int ok = (d_old_new < mes_var.tol);
            /* for a non-zero twist the twisted propagator must differ from the untwisted one */
            int twisted = (th[0] != 0. || th[1] != 0. || th[2] != 0.);
            int effect_ok = twisted ? (d_new_ref > 1e3 * mes_var.tol) : (d_new_ref < mes_var.tol);

            lprintf("TWIST", 0,
                    "p=(%d,%d,%d,%d) tw=%+d theta/pi=(%g,%g,%g): |old-new|/|new| = %.3e %s, |new-untwisted|/|new| = %.3e %s\n",
                    p[0], p[1], p[2], p[3], tw, th[0], th[1], th[2], d_old_new, ok ? "OK" : "MISMATCH", d_new_ref,
                    effect_ok ? "OK" : (twisted ? "NO EFFECT" : "MISMATCH"));

            if (!ok || !effect_ok) { return_value = 1; }
        }
    }

    /* leave everything untwisted */
    restore_links(u_save);
    set_eitheta(0., 0., 0.);

    global_sum_int(&return_value, 1); /* fail on every process if any failed */
    if (return_value != 0) { return_value = 1; }

    lprintf("MAIN", 0, "%s\n", return_value == 0 ? "Old and new twisting agree: PASSED" : "Test FAILED");

    free_propagator_eo();
    free_spinor_field(tmp);
    free_spinor_field(prop_ref);
    free_spinor_field(prop_old);
    free_spinor_field(prop_new);
    free_spinor_field(source);
    free_suNf_field(u_twisted);
    free_suNf_field(u_pristine);
    free_suNf_field(u_save);

    finalize_process();
    return return_value;
}
