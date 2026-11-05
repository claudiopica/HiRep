/*******************************************************************************
 * NOCOMPILE= BASIC_SF
 * NOCOMPILE= BC_T_SF
 * NOCOMPILE= ROTATED_SF
 * NOCOMPILE= BC_T_THETA
 * NOCOMPILE= FERMION_THETA && !BC_X_THETA
 * NOCOMPILE= FERMION_THETA && !BC_Y_THETA
 * NOCOMPILE= FERMION_THETA && !BC_Z_THETA
 *
 * RI'-MOM self-test 1: tree level (CHECK_Z).
 *
 * On the free field U = 1 (Landau gauge fixed), for each momentum l*p_in,
 * l = 1..n_mom, and twist tw = -n_twist..n_twist: Zq = 1 and, for all 16
 * local bilinears (S, P, V, A, T), the amputated vertex equals Gamma exactly,
 * so Z_Gamma = 1. Gamma is built from the single-gamma macros following the
 * channel name, independently of the composite macros used in
 * measure_renormalization.c. The twist enters through the effective momentum
 * p + theta/2. Twisting as in the build (see check_RIMOM_utils.h).
 *
 * Run with MPI (check_RIMOM_1.in: NP_X = 2).
 *
 * Originally written by Rudy Arthur
 * Test by Vincent Drach
 *******************************************************************************/

#define MAIN_PROGRAM

#include "check_RIMOM_utils.h"

static const double TOL_TREE_LEVEL = 1e-6;

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
 * Tree-level checks
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
static int run_tree_level_scan(rimom_ctx *ctx, double tol) {
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
            apply_twist(u_untwisted, twist * mom[1], twist * mom[2], twist * mom[3]);
            lprintf("TEST", 0, "<p> after twist %1.6f (must stay 1)\n", avr_plaquette());
            for (int mu = 0; mu < 4; mu++) {
                p[mu] = l * mom[mu];
                p_eff[mu] = p[mu] + (mu == 0 ? 0.0 : 0.5 * twist * mom[mu]);
            }

            char label[256];
            sprintf(label, "NPR mom_idx %d twist %d ", num, tw);
            data_storage_array *out =
                measure_all_channels(ctx->source, ctx->prop_in, ctx->prop_out, ctx->nm, ctx->m, p, 0, label);
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

    apply_twist(u_untwisted, 0., 0., 0.);
    free_suNf_field(u_untwisted);
    return ok;
}

int main(int argc, char *argv[]) {
    rimom_ctx ctx;
    rimom_setup(&argc, &argv, "check_RIMOM_1", &ctx);
    int return_value = run_tree_level_scan(&ctx, TOL_TREE_LEVEL) ? 0 : 1;
    return rimom_finish(&ctx, return_value, "check_RIMOM_1", "tree-level Zq = 1 and Z = 1 for all local channels");
}
