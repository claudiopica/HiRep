/*******************************************************************************
 * Shared code for the RI'-MOM self-tests check_RIMOM_1..3
 * (LibHR/Observables/measure_renormalization.c as used by
 * RenormalizationFactors/measure_Z_mom.c: gauge-fixed momentum sources,
 * twisted BCs, exceptional kinematics p_in = p_out).
 *
 * Header-only (static functions), included by exactly one .c per executable,
 * so no change to Make/build.pl is needed.
 *
 * Twisted BCs psi(x + L_mu) = exp(i*pi*theta_mu) psi(x) in the spatial
 * directions (a momentum shift of theta_mu/2 in units of 2*pi/L_mu) are
 * applied as in the build at hand (apply_twist):
 *  - theta build (-x T -y T -z T, FERMION_THETA): through eitheta, as in
 *    measure_Z_mom.c. This is the production path, including the eitheta
 *    factor on the conserved-current link (vertex_link).
 *  - otherwise: by multiplying the spatial links (former measure_Z_mom.c
 *    method).
 * Every test must pass in both builds.
 *
 * Each .c must carry the NOCOMPILE lines (they are read from the .c only):
 *   BASIC_SF, BC_T_SF, ROTATED_SF, BC_T_THETA,
 *   FERMION_THETA && !BC_X_THETA (and Y, Z)
 *
 * Originally written by Rudy Arthur
 * Tests by Vincent Drach
 *******************************************************************************/

#ifndef CHECK_RIMOM_UTILS_H
#define CHECK_RIMOM_UTILS_H

#include "libhr.h"
#include <string.h>

#ifdef FERMION_THETA
#if !defined(BC_X_THETA) || !defined(BC_Y_THETA) || !defined(BC_Z_THETA)
#error check_RIMOM with theta BCs needs them in X, Y and Z (write_mkflags.pl -x T -y T -z T)
#endif
#endif

#define PROP_N (4 * NF)

/* Vertex channels checked: the 16 local bilinears first, then the conserved
 * vector current. Names as in measure_renormalization.c. */
#define N_LOCAL_CHANNELS 16
#define N_VERTEX_CHANNELS 20
static const char *vertex_channels[N_VERTEX_CHANNELS] = { "id",     "g5",     "g0",     "g1",   "g2",   "g3",   "g0g1",
                                                          "g0g2",   "g0g3",   "g5g0",   "g5g1", "g5g2", "g5g3", "g5g0g1",
                                                          "g5g0g2", "g5g0g3", "cg0",    "cg1",  "cg2",  "cg3" };

/*******************************************************************************
 * Input
 *******************************************************************************/

typedef struct _input_renormalization {
    char mstring[256];
    double precision;
    int ne;
    int n_mom;
    int n_twist;
    int pt_in;
    int px_in;
    int py_in;
    int pz_in;

    /* for the reading function */
    input_record_t read[10];
} input_renormalization;

#define init_input_renormalization(varname)                                                         \
    {                                                                                               \
        .read = {                                                                                   \
            { "quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring },            \
            { "inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision },        \
            { "non-exceptional kinematics (not supported)", "mes:ne = %d", INT_T, &(varname).ne },  \
            { "number of momenta", "mes:n_mom = %d", INT_T, &(varname).n_mom },                     \
            { "number of twists at each momentum", "mes:n_twist = %d", INT_T, &(varname).n_twist }, \
            { "momentum t component", "mes:pt_in = %d", INT_T, &(varname).pt_in },                  \
            { "momentum x component", "mes:px_in = %d", INT_T, &(varname).px_in },                  \
            { "momentum y component", "mes:py_in = %d", INT_T, &(varname).py_in },                  \
            { "momentum z component", "mes:pz_in = %d", INT_T, &(varname).pz_in },                  \
            { NULL, NULL, INT_T, NULL }                                                             \
        }                                                                                           \
    }

char input_filename[256] = "input_file";
char output_filename[256] = "renormalization.out";

input_renormalization mes_var = init_input_renormalization(mes_var);

/*******************************************************************************
 * Twisted BCs
 *******************************************************************************/

#ifndef FERMION_THETA
/* Multiply the spatial links of u_gauge_f, local sites and MPI buffers, by
 * exp(i*pi*theta_mu/L_mu): a momentum shift of theta_mu/2 in units of 2*pi/L_mu. */
static inline void twist_XYZ_bc(double theta_x, double theta_y, double theta_z) {
    hr_complex eith_x = cexp(I * PI * theta_x / (double)GLB_X);
    hr_complex eith_y = cexp(I * PI * theta_y / (double)GLB_Y);
    hr_complex eith_z = cexp(I * PI * theta_z / (double)GLB_Z);
    suNf *u, utmp;

    for (int it = 0; it < T_EXT; ++it) {
        for (int ix = 0; ix < X_EXT; ++ix) {
            for (int iy = 0; iy < Y_EXT; ++iy) {
                for (int iz = 0; iz < Z_EXT; ++iz) {
                    int index = ipt_ext(it, ix, iy, iz);
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

/* Restore the untwisted links into u_gauge_f, including the MPI boundary
 * buffers: copy_suNf_field only copies local sites, while twist_XYZ_bc also
 * multiplies the buffers. Without the sendrecv, the buffers keep the previous
 * twist and the Dirac operator is no longer gamma5-hermitian across process
 * boundaries (wrong results, solver can hang). */
static inline void restore_untwisted_links(suNf_field *u_untwisted) {
    copy_suNf_field(u_gauge_f, u_untwisted);
#ifdef WITH_MPI
    start_sendrecv_suNf_field(u_gauge_f);
    complete_sendrecv_suNf_field(u_gauge_f);
#endif
}
#endif /* !FERMION_THETA */

/* Set the twist theta (units of pi) in X, Y, Z; T is never twisted.
 * theta build: eitheta, the links stay untwisted (u_untwisted unused).
 * otherwise:   restore the untwisted links, then multiply them. */
static inline void apply_twist(suNf_field *u_untwisted, double theta_x, double theta_y, double theta_z) {
#ifdef FERMION_THETA
    (void)u_untwisted;
    eitheta[0] = 1.0;
    eitheta[1] = cexp(I * PI * theta_x / (double)GLB_X);
    eitheta[2] = cexp(I * PI * theta_y / (double)GLB_Y);
    eitheta[3] = cexp(I * PI * theta_z / (double)GLB_Z);
#else
    restore_untwisted_links(u_untwisted);
    twist_XYZ_bc(theta_x, theta_y, theta_z);
#endif
}

/*******************************************************************************
 * Measurement
 *******************************************************************************/

/* Rotate n spinor fields by the gauge transformation g(x) left behind by
 * random_gauge_transform() (gaugefix.c keeps it in the global g, slot mu=0):
 * eta(x) -> R[g(x)] eta(x), R = fermion representation. Together with
 * U -> g U g^dag this makes every propagator column transform as
 * psi(x) -> R[g(x)] psi(x). */
static inline void gauge_rotate_spinor_fields(spinor_field *s, int n) {
    extern suNg_field *g;
    for (int k = 0; k < n; k++) {
        _MASTER_FOR(&glattice, ix) {
            suNf gf;
            suNf_vector v;
            _group_represent2(&gf, _4FIELD_AT(g, ix, 0));
            suNf_spinor *sp = _FIELD_AT(&s[k], ix);
            for (int alpha = 0; alpha < 4; alpha++) {
                _suNf_multiply(v, gf, sp->c[alpha]);
                sp->c[alpha] = v;
            }
        }
        start_sendrecv_spinor_field(&s[k]);
        complete_sendrecv_spinor_field(&s[k]);
    }
}

/* Momentum sources (optionally rotated by g), propagators and all vertex
 * channels at exceptional kinematics (prop_out = prop_in). Prints the
 * channels with `label` and returns the STORE output (caller frees it with
 * free_data_storage).
 * print_renormalization is also what resets the accumulators of
 * measure_renormalization.c: without it, the next call would add on top. */
static inline data_storage_array *measure_all_channels(spinor_field *source, spinor_field *prop_in,
                                                       spinor_field *prop_out, int nm, double *m, double p[4],
                                                       int rotate, const char *label) {
    data_storage_array *out = NULL;
    for (int k = 0; k < NF; ++k) {
        create_gauge_fixed_momentum_source(source, p[0], p[1], p[2], p[3], k);
        if (rotate) { gauge_rotate_spinor_fields(source, 4); }
        calc_propagator(prop_in + 4 * k, source, 4);
    }
    for (int j = 0; j < 4 * NF; j++) {
        copy_spinor_field(&prop_out[j], &prop_in[j]);
    }
    measure_renormalization(prop_in, prop_out, nm, p[0], p[1], p[2], p[3], p[0], p[1], p[2], p[3], STORE, &out);
    print_renormalization(0, nm, m, (char *)label, p[0], p[1], p[2], p[3], p[0], p[1], p[2], p[3]);
    return out;
}

/* Read channel `name`, mass index 0, out of a STORE-filled data_storage_array
 * (already global-summed and 1/GLB_VOLUME normalised). */
static inline void channel_from_storage(data_storage_array *dat, const char *name, suNf_propagator *out) {
    int channel = renorm_channel_index(name);
    error(channel < 0, 1, "check_RIMOM", "unknown channel name (mismatch with measure_renormalization.c?)");
    for (int r = 0; r < PROP_N; r++) {
        for (int c = 0; c < PROP_N; c++) {
            int idx_re[5] = { channel, 0, r, c, 0 };
            int idx_im[5] = { channel, 0, r, c, 1 };
            _PROP_IDX(*out, r, c) = *data_storage_element(dat, 0, idx_re) + I * (*data_storage_element(dat, 0, idx_im));
        }
    }
}

/* Compare channels of two STORE outputs: max|A-B|/max|A| <= tol for each. */
static inline int compare_channels(data_storage_array *outA, data_storage_array *outB, const char **names, int n_names,
                                   double tol, const char *tag) {
    int ok = 1;
    for (int c = 0; c < n_names; c++) {
        suNf_propagator A, B;
        channel_from_storage(outA, names[c], &A);
        channel_from_storage(outB, names[c], &B);
        double diff = 0., norm = 0.;
        for (int r = 0; r < PROP_N; r++) {
            for (int s = 0; s < PROP_N; s++) {
                diff = fmax(diff, cabs(_PROP_IDX(A, r, s) - _PROP_IDX(B, r, s)));
                norm = fmax(norm, cabs(_PROP_IDX(A, r, s)));
            }
        }
        int pass = (diff <= tol * norm);
        lprintf(tag, 0, "%-7s max|A-B|/max|A| = %.3e  %s\n", names[c], diff / norm, pass ? "ok" : "FAILED");
        if (!pass) { ok = 0; }
    }
    return ok;
}

/*******************************************************************************
 * Common main(): setup and teardown
 *******************************************************************************/

/* Fields shared by every test; one mass (the propagator layout prop + 4*k
 * assumes it). */
typedef struct {
    int nm;
    double m[1];
    spinor_field *source, *prop_in, *prop_out;
} rimom_ctx;

static inline void rimom_setup(int *argc, char ***argv, const char *name, rimom_ctx *ctx) {
    setup_process(argc, argv);
    setup_gauge_fields();
    read_input(mes_var.read, get_input_filename());

    lprintf("MAIN", 0, "%s\n", name);
    lprintf("MAIN", 0, "Compiled with macros: %s\n", MACROS);
    lprintf("MAIN", 0, "PId =  %d [world_size: %d]\n\n", PID, WORLD_SIZE);
    lprintf("MAIN", 0, "input file [%s]\n", input_filename);
    lprintf("MAIN", 0, "output file [%s]\n", output_filename);
#ifdef GAUGE_SON
    lprintf("MAIN", 0, "Gauge group: SO(%d)\n", NG);
#else
    lprintf("MAIN", 0, "Gauge group: SU(%d)\n", NG);
#endif
    lprintf("MAIN", 0, "Fermion representation: " REPR_NAME " [dim=%d]\n", NF);
#ifdef FERMION_THETA
    lprintf("MAIN", 0, "Twisted BCs via eitheta (FERMION_THETA)\n");
#else
    lprintf("MAIN", 0, "Twisted BCs via twisted links\n");
#endif

    error(mes_var.ne != 0, 1, name,
          "mes:ne = 1 (non-exceptional kinematics) is not supported by this test: the checks need p_in = p_out. "
          "Set mes:ne = 0, or use RenormalizationFactors/measure_Z_mom for non-exceptional measurements.");

    ctx->nm = 1;
    ctx->m[0] = -atof(mes_var.mstring); /* note: measure_Z_mom uses +masses */
    lprintf("MAIN", 0, "Inverter precision = %e\n", mes_var.precision);
    lprintf("MAIN", 0, "Mass[0] = %f\n", ctx->m[0]);
    lprintf("MAIN", 0, "Number of momenta n_mom = %d, twists per momentum: -%d..%d\n", mes_var.n_mom, mes_var.n_twist,
            mes_var.n_twist);
    lprintf("MAIN", 0, "Momentum (%d, %d, %d, %d)\n", mes_var.pt_in, mes_var.px_in, mes_var.py_in, mes_var.pz_in);

    ctx->source = alloc_spinor_field(4, &glattice);
    ctx->prop_in = alloc_spinor_field(4 * ctx->nm * NF, &glattice);
    ctx->prop_out = alloc_spinor_field(4 * ctx->nm * NF, &glattice);
    init_propagator_eo(ctx->nm, ctx->m, mes_var.precision);
}

/* Prints PASS/FAIL, frees everything, returns the (global) return value. */
static inline int rimom_finish(rimom_ctx *ctx, int return_value, const char *name, const char *what) {
    global_sum_int(&return_value, 1);
    if (return_value == 0) {
        lprintf("MAIN", 0, "%s: PASS (%s)\n", name, what);
    } else {
        lprintf("MAIN", 0, "%s: FAIL (see the lines marked FAILED above)\n", name);
    }

    free_spinor_field(ctx->source);
    free_spinor_field(ctx->prop_in);
    free_spinor_field(ctx->prop_out);
    free_propagator_eo();
    free_BCs();
    free_suNg_field(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
    free_suNf_field(u_gauge_f);
#endif

    finalize_process();
    lprintf("MAIN", 0, "return_value= %d\n ", return_value);
    return return_value;
}

#endif /* CHECK_RIMOM_UTILS_H */
