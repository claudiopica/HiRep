/*******************************************************************************
 * NOCOMPILE= BASIC_SF
 * NOCOMPILE= BC_T_SF
 * NOCOMPILE= ROTATED_SF
 * NOCOMPILE= FERMION_THETA
 *
 * Self-test of the RI'-MOM vertex code (LibHR/Observables/measure_renormalization.c)
 * as used by RenormalizationFactors/measure_Z_mom.c (gauge-fixed momentum
 * sources, twisted links, exceptional kinematics p_in = p_out).
 *
 * 1. Tree level (CHECK_Z). On the free field U = 1 (Landau gauge fixed), for
 *    each momentum l*p_in, l = 1..n_mom, and twist tw = -n_twist..n_twist:
 *    Zq = 1 and, for all 16 local bilinears (S, P, V, A, T), the amputated
 *    vertex equals Gamma exactly, so Z_Gamma = 1. Gamma is built from the
 *    single-gamma macros following the channel name, independently of the
 *    composite macros used in measure_renormalization.c. The twist enters
 *    through the effective momentum p + theta/2.
 *
 * 2. Gauge covariance (CHECK_GT). On a random configuration, U -> g U g^dag
 *    together with eta -> g eta must leave every vertex channel unchanged
 *    (Sin/Sout excluded: the propagator is gauge variant). This covers the
 *    conserved current cg0..cg3 in all four directions.
 *
 * 3. Twist equivalence (CHECK_TW). On a random configuration, momentum n with
 *    twist theta = +n and momentum 2n with theta = -n reach the same
 *    p + theta/2: every channel, Sin/Sout included, must be identical. Tests
 *    the twisted links and the Fourier phases of the vertices (twisted
 *    directions only, so T is covered by check 2 alone).
 *
 * Run with MPI (see check_RIMOM.in: NP_X = 2) so that the twisted links,
 * the gauge transformation and the conserved current cross rank boundaries.
 *
 * Originally written by Rudy Arthur
 * Test by Vincent Drach
 *******************************************************************************/

#define MAIN_PROGRAM

#include "libhr.h"
#include <string.h>

#define PROP_N (4 * NF)

static const double TOL_TREE_LEVEL = 1e-6;
static const double TOL_GAUGE_COVARIANCE = 1e-8;

/* Vertex channels checked here: the 16 local bilinears first, then the
 * conserved vector current. Names as in measure_renormalization.c. */
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

#define init_input_renormalization(varname)                                                      \
    {                                                                                            \
        .read = {                                                                                \
            { "quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring },         \
            { "inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision },     \
            { "non-exceptional kinematics (not supported)", "mes:ne = %d", INT_T, &(varname).ne }, \
            { "number of momenta", "mes:n_mom = %d", INT_T, &(varname).n_mom },                  \
            { "number of twists at each momentum", "mes:n_twist = %d", INT_T, &(varname).n_twist }, \
            { "momentum t component", "mes:pt_in = %d", INT_T, &(varname).pt_in },               \
            { "momentum x component", "mes:px_in = %d", INT_T, &(varname).px_in },               \
            { "momentum y component", "mes:py_in = %d", INT_T, &(varname).py_in },               \
            { "momentum z component", "mes:pz_in = %d", INT_T, &(varname).pz_in },               \
            { NULL, NULL, INT_T, NULL }                                                          \
        }                                                                                        \
    }

char input_filename[256] = "input_file";
char output_filename[256] = "renormalization.out";

input_renormalization mes_var = init_input_renormalization(mes_var);

/*******************************************************************************
 * Twisted links (same convention as measure_Z_mom.c)
 *******************************************************************************/

/* Multiply the spatial links of u_gauge_f, local sites and MPI buffers, by
 * exp(i*pi*theta_mu/L_mu): a momentum shift of theta_mu/2 in units of 2*pi/L_mu. */
static void twist_XYZ_bc(double theta_x, double theta_y, double theta_z) {
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
static void restore_untwisted_links(suNf_field *u_untwisted) {
    copy_suNf_field(u_gauge_f, u_untwisted);
#ifdef WITH_MPI
    start_sendrecv_suNf_field(u_gauge_f);
    complete_sendrecv_suNf_field(u_gauge_f);
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
static void gauge_rotate_spinor_fields(spinor_field *s, int n) {
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
static data_storage_array *measure_all_channels(spinor_field *source, spinor_field *prop_in, spinor_field *prop_out,
                                                int nm, double *m, double p[4], int rotate, const char *label) {
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

/*******************************************************************************
 * 12x12 (spin x colour) matrix helpers
 *******************************************************************************/

/* Generic complex Gauss-Jordan inverse with partial pivoting, operating on
 * the full spin(x)colour object via _PROP_IDX.
 *-- the accompanying check_matrix_inverse() residual
 * check is there specifically to catch a bug in this function  */
static int propagator_inverse(suNf_propagator *out, suNf_propagator *in)
{
  int i, j, k, piv;
  hr_complex a[PROP_N][PROP_N], id[PROP_N][PROP_N], tmp, factor;

  for (i = 0; i < PROP_N; i++)
    for (j = 0; j < PROP_N; j++)
    {
      a[i][j] = _PROP_IDX(*in, i, j);
      id[i][j] = (i == j) ? 1.0 : 0.0;
    }

  for (k = 0; k < PROP_N; k++)
  {
    piv = k;
    for (i = k + 1; i < PROP_N; i++)
    {
      if (cabs(a[i][k]) > cabs(a[piv][k])) piv = i;
    }
    if (cabs(a[piv][k]) < 1e-30) return -1; /* singular */
    if (piv != k)
    {
      for (j = 0; j < PROP_N; j++)
      {
        tmp = a[k][j]; a[k][j] = a[piv][j]; a[piv][j] = tmp;
        tmp = id[k][j]; id[k][j] = id[piv][j]; id[piv][j] = tmp;
      }
    }
    factor = a[k][k];
    for (j = 0; j < PROP_N; j++)
    {
      a[k][j] /= factor;
      id[k][j] /= factor;
    }
    for (i = 0; i < PROP_N; i++)
    {
      if (i == k) continue;
      factor = a[i][k];
      for (j = 0; j < PROP_N; j++)
      {
        a[i][j] -= factor * a[k][j];
        id[i][j] -= factor * id[k][j];
      }
    }
  }
  for (i = 0; i < PROP_N; i++)
    for (j = 0; j < PROP_N; j++)
      _PROP_IDX(*out, i, j) = id[i][j];
  return 0;
}

/* max|S @ Sinv - 1|, should sit near solver precision; a bad inverse would
 * show up here before it corrupts the Zq result below. */
static double check_matrix_inverse(suNf_propagator *S, suNf_propagator *Sinv)
{
  suNf_propagator prod;
  double maxdev = 0.0;
  int i, j;
  _propagator_mul(prod, (*S), (*Sinv));
  for (i = 0; i < PROP_N; i++)
    for (j = 0; j < PROP_N; j++)
    {
      hr_complex expect = (i == j) ? 1.0 : 0.0;
      double dev = cabs(_PROP_IDX(prod, i, j) - expect);
      if (dev > maxdev) maxdev = dev;
    }
  return maxdev;
}

/* Read channel `name`, mass index 0, out of a STORE-filled data_storage_array
 * (already global-summed and 1/GLB_VOLUME normalised). */
static void channel_from_storage(data_storage_array *dat, const char *name, suNf_propagator *out) {
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

/* Gamma (x 1_colour) for a channel name made of "g0".."g3","g5" tokens,
 * e.g. "g5g0g1" = g5*g0*g1, or "id" = 1. Built from the single-gamma macros
 * only (applied right to left), so it is independent of the composite
 * _g5g0g1_propagator etc. that measure_renormalization.c uses.
 * Returns 0 on success, -1 for an unknown name. */
static int build_gamma(suNf_propagator *G, const char *name) {
    suNf_propagator cur, tmp;
    _propagator_one(cur);
    if (strcmp(name, "id") != 0) {
        int n = (int)strlen(name);
        if (n == 0 || n % 2 != 0) { return -1; }
        for (int k = n - 2; k >= 0; k -= 2) {
            if (name[k] != 'g') { return -1; }
            switch (name[k + 1]) {
            case '0': _g0_propagator(tmp, cur); break;
            case '1': _g1_propagator(tmp, cur); break;
            case '2': _g2_propagator(tmp, cur); break;
            case '3': _g3_propagator(tmp, cur); break;
            case '5': _g5_propagator(tmp, cur); break;
            default: return -1;
            }
            cur = tmp;
        }
    }
    *G = cur;
    return 0;
}

/* Tr[ sum_mu gamma_mu sin(p_mu) * Sinv ] */
static hr_complex trace_gamma_sin_dot_Sinv(suNf_propagator *Sinv, double sp[4]) {
    suNf_propagator g0S, g1S, g2S, g3S;
    hr_complex tr0, tr1, tr2, tr3;
    _g0_propagator(g0S, *Sinv);
    _g1_propagator(g1S, *Sinv);
    _g2_propagator(g2S, *Sinv);
    _g3_propagator(g3S, *Sinv);
    _propagator_trace(tr0, g0S);
    _propagator_trace(tr1, g1S);
    _propagator_trace(tr2, g2S);
    _propagator_trace(tr3, g3S);
    return sp[0] * tr0 + sp[1] * tr1 + sp[2] * tr2 + sp[3] * tr3;
}

/*******************************************************************************
 * Check 1: tree level
 *******************************************************************************/

/* One local channel: the amputated vertex Lambda = Sout^-1 G Sin^-1 must
 * equal Gamma exactly, hence Z_Gamma = 4 NF Zq / Tr[Lambda Gamma^-1] = 1.
 * Also checks the full matrix (max|Lambda - Gamma|), which catches pieces
 * the trace projection misses. */
static int check_local_channel(suNf_propagator *Sin_inv, suNf_propagator *Sout_inv, suNf_propagator *vertex,
                               const char *name, hr_complex zq, double tol) {
    suNf_propagator G, Ginv, tmp, lambda, lg;
    hr_complex tr;
    if (build_gamma(&G, name) != 0 || propagator_inverse(&Ginv, &G) != 0) {
        lprintf("CHECK_Z", 0, "%s: cannot build Gamma from the channel name\n", name);
        return 0;
    }
    _propagator_mul(tmp, (*Sout_inv), (*vertex));
    _propagator_mul(lambda, tmp, (*Sin_inv));
    _propagator_mul(lg, lambda, Ginv);
    _propagator_trace(tr, lg);
    hr_complex z = 4.0 * (double)NF * zq / tr;

    double dev = 0.;
    for (int r = 0; r < PROP_N; r++) {
        for (int s = 0; s < PROP_N; s++) {
            dev = fmax(dev, cabs(_PROP_IDX(lambda, r, s) - _PROP_IDX(G, r, s)));
        }
    }
    int pass = (cabs(z - 1.0) < tol) && (dev < tol);
    lprintf("CHECK_Z", 0, "Z_%-7s = %.10f %+.10fi  max|Lambda-Gamma| = %.2e  %s\n", name, creal(z), cimag(z), dev,
            pass ? "ok" : "FAILED");
    return pass;
}

/* Zq and all 16 local channels at one momentum point. p_eff: momentum seen
 * by the Dirac operator in units of 2*pi/L_mu (twist folded in).
 * Returns 1 (pass) or 0 (fail). */
static int check_tree_level_z(data_storage_array *dat, double p_eff[4], double tol) {
    const double glb[4] = { GLB_T, GLB_X, GLB_Y, GLB_Z };
    double sp[4], p2 = 0.0;
    suNf_propagator Sin, Sout, Sin_inv, Sout_inv;
    int ok = 1;

    for (int mu = 0; mu < 4; mu++) {
        sp[mu] = sin(2.0 * PI * p_eff[mu] / glb[mu]);
        p2 += sp[mu] * sp[mu];
    }
    if (p2 < 1e-12) {
        lprintf("CHECK_Z", 0, "p=0 at this point, skipping (expected for some twist/p_ref combinations)\n");
        return 1;
    }

    channel_from_storage(dat, "Sin", &Sin);
    channel_from_storage(dat, "Sout", &Sout);
    if (propagator_inverse(&Sin_inv, &Sin) != 0 || propagator_inverse(&Sout_inv, &Sout) != 0) {
        lprintf("CHECK_Z", 0, "FAIL: Sin or Sout singular\n");
        return 0;
    }
    double residual = check_matrix_inverse(&Sin, &Sin_inv);
    if (residual > 1e-8) {
        lprintf("CHECK_Z", 0, "WARNING: Sin inverse residual = %e (expected ~solver precision)\n", residual);
    }

    /* Zq = (-i/(4*NF)) * Tr[ sum_mu gamma_mu sin(p_mu) Sinv ] / sum_mu sin(p_mu)^2 */
    hr_complex zq = -I * trace_gamma_sin_dot_Sinv(&Sin_inv, sp) / (4.0 * (double)NF * p2);
    lprintf("CHECK_Z", 0, "Zq = %.10f %+.10fi  (expect 1+0i)\n", creal(zq), cimag(zq));
    if (cabs(zq - 1.0) > tol) { ok = 0; }

    for (int c = 0; c < N_LOCAL_CHANNELS; c++) {
        suNf_propagator vertex;
        channel_from_storage(dat, vertex_channels[c], &vertex);
        if (!check_local_channel(&Sin_inv, &Sout_inv, &vertex, vertex_channels[c], zq, tol)) { ok = 0; }
    }
    return ok;
}

/* Free field U = 1 in Landau gauge; scan momenta l*p_in, l = 1..n_mom, and
 * twists tw = -n_twist..n_twist. Returns 1 if every point passes. */
static int run_tree_level_scan(spinor_field *source, spinor_field *prop_in, spinor_field *prop_out, int nm, double *m,
                               double tol) {
    const double mom[4] = { mes_var.pt_in, mes_var.px_in, mes_var.py_in, mes_var.pz_in };
    struct timeval start, end, etime;
    int ok = 1, num = 0;

    lprintf("MAIN", 0, "Generating a unit configuration\n");
    unit_gauge(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "initial plaq %1.6f\n", calc_plaq(u_gauge));

    gettimeofday(&start, 0);
    double act = gaugefix(10, /* Landau (0..3 would be Coulomb) */
                          1.8, /* overrelaxation */
                          10000, /* max iterations */
                          1e-10, /* tolerance */
                          u_gauge);
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MAIN", 0, "Unit configuration gauge fixed in [%ld sec %ld usec], action %1.6f, plaq %1.6f\n",
            etime.tv_sec, etime.tv_usec, act, calc_plaq(u_gauge));
    represent_gauge_field();

    suNf_field *u_untwisted = alloc_suNf_field(&glattice);
    copy_suNf_field(u_untwisted, u_gauge_f);

    gettimeofday(&start, 0);
    for (int l = 1; l <= mes_var.n_mom; ++l) {
        for (int tw = -mes_var.n_twist; tw <= mes_var.n_twist; tw++) {
            /* theta_mu = (tw/2) * p_in_mu on the spatial links (T untwisted):
             * the Dirac operator sees p_eff = l*p_in + theta/2 */
            double twist = 0.5 * tw;
            double p[4], p_eff[4];
            restore_untwisted_links(u_untwisted);
            twist_XYZ_bc(twist * mom[1], twist * mom[2], twist * mom[3]);
            lprintf("TEST", 0, "<p> after twist %1.6f (must stay 1)\n", avr_plaquette());
            for (int mu = 0; mu < 4; mu++) {
                p[mu] = l * mom[mu];
                p_eff[mu] = p[mu] + (mu == 0 ? 0.0 : 0.5 * twist * mom[mu]);
            }

            char label[256];
            sprintf(label, "NPR mom_idx %d twist %d ", num, tw);
            data_storage_array *out = measure_all_channels(source, prop_in, prop_out, nm, m, p, 0, label);
            int pass = check_tree_level_z(out, p_eff, tol);
            free_data_storage(out);

            lprintf("CHECK_Z", 0, "mom_idx %d twist %d: tree-level Z check %s\n", num, tw, pass ? "PASSED" : "FAILED");
            if (!pass) { ok = 0; }
            num++;
        }
    }
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MAIN", 0, "Unit configuration analysed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

    restore_untwisted_links(u_untwisted);
    free_suNf_field(u_untwisted);
    return ok;
}

/*******************************************************************************
 * Check 2: gauge covariance
 *******************************************************************************/

/* Compare channels of two STORE outputs: max|A-B|/max|A| <= tol for each. */
static int compare_channels(data_storage_array *outA, data_storage_array *outB, const char **names, int n_names,
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

static int check_gauge_covariance(spinor_field *source, spinor_field *prop_in, spinor_field *prop_out, int nm, double *m,
                                  double p[4], double tol) {
    /* random links: not gauge-equivalent to U=1, so U != U^dag matters */
    random_u(u_gauge);
    represent_gauge_field(); /* also completes random_u's halo exchange */
    double plaq_before = avr_plaquette();
    lprintf("CHECK_GT", 0, "random configuration, plaquette %1.6f\n", plaq_before);

    data_storage_array *outA = measure_all_channels(source, prop_in, prop_out, nm, m, p, 0, "GT check A");

    random_gauge_transform(u_gauge); /* U -> g U g^dag, g left in gaugefix.c's global g */
    represent_gauge_field();
    lprintf("CHECK_GT", 0, "after random gauge transformation, plaquette %1.6f (must be unchanged)\n", avr_plaquette());
    if (fabs(avr_plaquette() - plaq_before) > 1e-12) {
        lprintf("CHECK_GT", 0, "plaquette changed: the gauge transformation is broken\n");
        free_data_storage(outA);
        return 0;
    }

    data_storage_array *outB = measure_all_channels(source, prop_in, prop_out, nm, m, p, 1, "GT check B");

    int ok = compare_channels(outA, outB, vertex_channels, N_VERTEX_CHANNELS, tol, "CHECK_GT");
    free_data_storage(outB);
    free_data_storage(outA);
    return ok;
}

/*******************************************************************************
 * Check 3: twist equivalence
 *******************************************************************************/

/* (p = n, theta = +n) and (p = 2n, theta = -n) give the same momentum
 * p + theta/2 = 3n/2 (units of 2*pi/L), n = the input momentum. The two
 * setups differ by the U(1) gauge transformation exp(2*pi*i n.x/L), which is
 * periodic, so every channel, Sin/Sout included, must be identical on any
 * configuration. Tests the twisted links (and their MPI buffers) and the
 * Fourier phases of the vertices, the conserved current included. Only the
 * spatial links are twisted, so n must have no T component. */
static int check_twist_equivalence(spinor_field *source, spinor_field *prop_in, spinor_field *prop_out, int nm,
                                   double *m, double tol) {
    const double n[4] = { mes_var.pt_in, mes_var.px_in, mes_var.py_in, mes_var.pz_in };
    if (n[0] != 0. || (n[1] == 0. && n[2] == 0. && n[3] == 0.)) {
        lprintf("CHECK_TW", 0, "skipped: needs a momentum with no T component and a non-zero spatial one\n");
        return 1;
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

    restore_untwisted_links(u_untwisted);
    twist_XYZ_bc(+n[1], +n[2], +n[3]);
    data_storage_array *outA = measure_all_channels(source, prop_in, prop_out, nm, m, pA, 0, "TW check A (p=n, theta=+n)");

    restore_untwisted_links(u_untwisted);
    twist_XYZ_bc(-n[1], -n[2], -n[3]);
    data_storage_array *outB = measure_all_channels(source, prop_in, prop_out, nm, m, pB, 0, "TW check B (p=2n, theta=-n)");

    static const char *propagators[2] = { "Sin", "Sout" };
    int ok = compare_channels(outA, outB, propagators, 2, tol, "CHECK_TW");
    if (!compare_channels(outA, outB, vertex_channels, N_VERTEX_CHANNELS, tol, "CHECK_TW")) { ok = 0; }

    free_data_storage(outB);
    free_data_storage(outA);
    restore_untwisted_links(u_untwisted);
    free_suNf_field(u_untwisted);
    return ok;
}

/*******************************************************************************
 * Main
 *******************************************************************************/

int main(int argc, char *argv[]) {
    int return_value = 0;
    int nm = 1; /* one mass: the propagator layout prop + 4*k assumes it */
    double m[1];

    setup_process(&argc, &argv);
    setup_gauge_fields();
    read_input(mes_var.read, get_input_filename());

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

    error(mes_var.ne != 0, 1, "check_RIMOM",
          "mes:ne = 1 (non-exceptional kinematics) is not supported by this test: the checks need p_in = p_out. "
          "Set mes:ne = 0, or use RenormalizationFactors/measure_Z_mom for non-exceptional measurements.");

    m[0] = -atof(mes_var.mstring); /* note: measure_Z_mom uses +masses */
    lprintf("MAIN", 0, "Inverter precision = %e\n", mes_var.precision);
    lprintf("MAIN", 0, "Mass[0] = %f\n", m[0]);
    lprintf("MAIN", 0, "Number of momenta n_mom = %d, twists per momentum: -%d..%d\n", mes_var.n_mom, mes_var.n_twist,
            mes_var.n_twist);
    lprintf("MAIN", 0, "Momentum (%d, %d, %d, %d)\n", mes_var.pt_in, mes_var.px_in, mes_var.py_in, mes_var.pz_in);

    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop_in = alloc_spinor_field(4 * nm * NF, &glattice);
    spinor_field *prop_out = alloc_spinor_field(4 * nm * NF, &glattice);
    init_propagator_eo(nm, m, mes_var.precision);

    if (!run_tree_level_scan(source, prop_in, prop_out, nm, m, TOL_TREE_LEVEL)) { return_value = 1; }

    double p_gt[4] = { mes_var.pt_in, mes_var.px_in, mes_var.py_in, mes_var.pz_in };
    if (!check_gauge_covariance(source, prop_in, prop_out, nm, m, p_gt, TOL_GAUGE_COVARIANCE)) { return_value = 1; }
    if (!check_twist_equivalence(source, prop_in, prop_out, nm, m, TOL_GAUGE_COVARIANCE)) { return_value = 1; }

    global_sum_int(&return_value, 1);
    if (return_value == 0) {
        lprintf("MAIN", 0,
                "check_RIMOM: PASS (tree-level Z = 1 for all local channels, gauge covariance, twist equivalence)\n");
    } else {
        lprintf("MAIN", 0, "check_RIMOM: FAIL (see the CHECK_Z / CHECK_GT / CHECK_TW lines marked FAILED above)\n");
    }

    free_spinor_field(source);
    free_spinor_field(prop_in);
    free_spinor_field(prop_out);
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
