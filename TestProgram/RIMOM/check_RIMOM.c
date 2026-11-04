/*******************************************************************************
 * NOCOMPILE= BASIC_SF
 * NOCOMPILE= BC_T_SF
 * NOCOMPILE= ROTATED_SF
 * NOCOMPILE= FERMION_THETA
 * Computation of Renomalization constants (Z_a,Z_q,Z_s,Z_ps,Z_t,Z_m,Z_v)  
 * factors with gauge fixed momentum sources. 
 *
 * Originally written by Rudy Arthur
 *
 * Test by Vincent Drach
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include "libhr.h"


static void twist_XYZ_bc(double theta_x, double theta_y, double theta_z)
{

  int index;
  int ix, iy, iz, it;
  suNf *u;
  suNf utmp;
  hr_complex eith_x, eith_y, eith_z;
  eith_x = cexp(I * PI * theta_x / (double)GLB_X);
  eith_y = cexp(I * PI * theta_y / (double)GLB_Y);
  eith_z = cexp(I * PI * theta_z / (double)GLB_Z);

  for (it = 0; it < T_EXT; ++it)
    for (ix = 0; ix < X_EXT; ++ix)
      for (iy = 0; iy < Y_EXT; ++iy)
        for (iz = 0; iz < Z_EXT; ++iz)
        {
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

/* Momentum sources (optionally rotated by g), propagators, all channels.
 * Exceptional kinematics: prop_out = prop_in. Returns the STORE output
 * (caller frees with free_data_storage). */
static data_storage_array *measure_all_channels(spinor_field *source, spinor_field *prop_in, spinor_field *prop_out,
                                                int nm, double *m, double p[4], int rotate, const char *label) {
    data_storage_array *out = NULL;
    for (int k = 0; k < NF; ++k) {
        create_gauge_fixed_momentum_source(source, p[0], p[1], p[2], p[3], k);
        if (rotate) gauge_rotate_spinor_fields(source, 4);
        calc_propagator(prop_in + 4 * k, source, 4);
    }
    for (int j = 0; j < 4 * NF; j++) {
        copy_spinor_field(&prop_out[j], &prop_in[j]);
    }
    measure_renormalization(prop_in, prop_out, nm, p[0], p[1], p[2], p[3], p[0], p[1], p[2], p[3], STORE, &out);
    // calling print renormalizatoin reset the accumulators.
    print_renormalization(0, nm, m, (char *)label, p[0], p[1], p[2], p[3], p[0], p[1], p[2], p[3]);
    return out;
}
/* Renormalization parameters */
typedef struct _input_renormalization
{
  char mstring[256];
  char configlist[256]; /* list of configuration */
  double precision;
  int ne;
  int n_mom;
  int n_twist;
  int pt_out;
  int px_out;
  int py_out;
  int pz_out;
  int pt_in;
  int px_in;
  int py_in;
  int pz_in;

  /* for the reading function */
  input_record_t read[15];
} input_renormalization;


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

#define init_input_renormalization(varname)                                                 \
  {                                                                                         \
    .read = {                                                                               \
      {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring},            \
      {"Configuration list:", "mes:configlist = %s", STRING_T, &(varname).configlist},      \
      {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},        \
      {"non-excepional configuration or not?", "mes:ne = %d", INT_T, &(varname).ne},        \
      {"number of momenta in each direction", "mes:n_mom = %d", INT_T, &(varname).n_mom},   \
      {"number of twists at each momentum", "mes:n_twist = %d", INT_T, &(varname).n_twist}, \
      {"mom1 t component", "mes:pt_out = %d", INT_T, &(varname).pt_out},                    \
      {"mom1 t component", "mes:px_out = %d", INT_T, &(varname).px_out},                    \
      {"mom1 t component", "mes:py_out = %d", INT_T, &(varname).py_out},                    \
      {"mom1 t component", "mes:pz_out = %d", INT_T, &(varname).pz_out},                    \
      {"mom2 t component", "mes:pt_in = %d", INT_T, &(varname).pt_in},                      \
      {"mom2 t component", "mes:px_in = %d", INT_T, &(varname).px_in},                      \
      {"mom2 t component", "mes:py_in = %d", INT_T, &(varname).py_in},                      \
      {"mom2 t component", "mes:pz_in = %d", INT_T, &(varname).pz_in},                      \
      {NULL, NULL, INT_T, NULL}                                                             \
    }                                                                                       \
  }

char cnfg_filename[256] = "";
char list_filename[256] = "";
char input_filename[256] = "input_file";
char output_filename[256] = "renormalization.out";
enum
{
  UNKNOWN_CNFG,
  DYNAMICAL_CNFG,
  QUENCHED_CNFG
};

input_renormalization mes_var = init_input_renormalization(mes_var);

typedef struct
{
  char string[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;


/*******************************************************************************
 * Tree-level RI'-MOM self-check (U=1 free field): Zq and Z_S/Z_P must equal
 * exactly 1 (up to solver precision) at tree level -- a defining property
 * of the RI'-MOM normalisation, not an approximation, independent of mass
 * or Wilson parameter (the mass term is proportional to the identity in
 * spin space and drops out of the trace).
 *
 * measure_renormalization(..., STORE, &out) now fills out with the
 * already global-summed, 1/GLB_VOLUME-normalised 3-point functions --
 * same storage_switch/data_storage_array convention as measure_spectrum_pt
 * (see meson_measurements.c / check_triplets_3.c). No manual
 * global_sum/normalise step needed here any more; just read the channels
 * we need out via data_storage_element, invert, amputate, trace.
 *
 * Only valid for EXCEPTIONAL kinematics (p_in == p_out); this check is
 * skipped automatically when mes_var.ne is set (see the call site below).
 ******************************************************************************/

#define PROP_N (4 * NF)

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

/* Read channel `channel`, mass index `mass_idx` out of a STORE-filled
 * data_storage_array into a plain suNf_propagator (already normalised --
 * see measure_renormalization's STORE branch). */
static void propagator_from_storage(data_storage_array *dat, int channel, int mass_idx, suNf_propagator *out)
{
  int r, c;
  for (r = 0; r < PROP_N; r++)
  {
    for (c = 0; c < PROP_N; c++)
    {
      int idx_re[5] = { channel, mass_idx, r, c, 0 };
      int idx_im[5] = { channel, mass_idx, r, c, 1 };
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
        if (n == 0 || n % 2 != 0) return -1;
        for (int k = n - 2; k >= 0; k -= 2) {
            if (name[k] != 'g') return -1;
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
static hr_complex trace_gamma_sin_dot_Sinv(suNf_propagator *Sinv, double sp[4])
{
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

/* Tree-level check of one local channel: the amputated vertex
 * Lambda = Sout^-1 G Sin^-1 must equal Gamma exactly, hence
 * Z_Gamma = 4 NF Zq / Tr[Lambda Gamma^-1] = 1. Also checks the full matrix
 * (max|Lambda - Gamma|), which catches pieces the trace projection misses. */
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


/* Sin/Sout/id_prop/g5_prop: already normalised (came straight out of
 * propagator_from_storage). p_lat: (pt,px,py,pz) for this momentum point,
 * twist already folded in. tol: e.g. 1e-6. Returns 1 (pass) or 0 (fail);
 * prints Zq/Z_S/Z_P either way. */
static int check_tree_level_z(suNf_propagator *Sin, suNf_propagator *Sout,
                              suNf_propagator *id_prop, suNf_propagator *g5_prop, data_storage_array *dat,
                              double p_lat[4], double tol)
{
  int mu, ok = 1;
  double sp[4], p2 = 0.0;
  double glb[4];
  suNf_propagator Sin_inv, Sout_inv, tmp1, tmp2, amputated;
  hr_complex zq, tr;
  double residual;

  glb[0] = GLB_T; glb[1] = GLB_X; glb[2] = GLB_Y; glb[3] = GLB_Z;
  for (mu = 0; mu < 4; mu++)
  {
    double p = 2.0 * PI * p_lat[mu] / glb[mu];
    sp[mu] = sin(p);
    p2 += sp[mu] * sp[mu];
  }
  if (p2 < 1e-12)
  {
    lprintf("CHECK_Z", 0, "p=0 at this point, skipping (expected for some twist/p_ref combinations)\n");
    return 1;
  }

  if (propagator_inverse(&Sin_inv, Sin) != 0 ||
      propagator_inverse(&Sout_inv, Sout) != 0)
  {
    lprintf("CHECK_Z", 0, "FAIL: Sin or Sout singular\n");
    return 0;
  }
  residual = check_matrix_inverse(Sin, &Sin_inv);
  if (residual > 1e-8)
  {
    lprintf("CHECK_Z", 0, "WARNING: Sin inverse residual = %e (expected ~solver precision)\n", residual);
  }

  /* Zq = (-i/(4*NF)) * Tr[ sum_mu gamma_mu sin(p_mu) Sinv ] / sum_mu sin(p_mu)^2 */
  tr = trace_gamma_sin_dot_Sinv(&Sin_inv, sp);
  zq = -I * tr / (4.0 * (double)NF * p2);
  lprintf("CHECK_Z", 0, "Zq = %.10f %+.10fi  (expect 1+0i)\n", creal(zq), cimag(zq));
  if (cabs(zq - 1.0) > tol) ok = 0;

  /* Z_S from the 'id' channel: amputate, trace against P_S = 1_spin.
   * NOTE: confirm this matches your P_Gamma normalisation convention
   * (eq 19) before trusting the overall constant here. */
  _propagator_mul(tmp1, Sout_inv, (*id_prop));
  _propagator_mul(amputated, tmp1, Sin_inv);
  _propagator_trace(tr, amputated);
  hr_complex z_s = 4.0 * (double)NF * zq / tr;
  lprintf("CHECK_Z", 0, "Z_S = %.10f %+.10fi  (expect 1+0i)\n", creal(z_s), cimag(z_s));
  if (cabs(z_s - 1.0) > tol) ok = 0;

  /* Z_P from the 'g5' channel, P_P = gamma5 */
  _propagator_mul(tmp1, Sout_inv, (*g5_prop));
  _propagator_mul(amputated, tmp1, Sin_inv);
  _g5_propagator(tmp2, amputated);
  _propagator_trace(tr, tmp2);
  hr_complex z_p = 4.0 * (double)NF * zq / tr;
  lprintf("CHECK_Z", 0, "Z_P = %.10f %+.10fi  (expect 1+0i)\n", creal(z_p), cimag(z_p));
  if (cabs(z_p - 1.0) > tol) ok = 0;

  /* all 16 local bilinears: tree-level Z_Gamma = 1 (V, A, T included) */
  static const char *local_names[] = { "id",   "g5",   "g0",   "g1",   "g2",     "g3",     "g0g1",  "g0g2",
                                       "g0g3", "g5g0", "g5g1", "g5g2", "g5g3",   "g5g0g1", "g5g0g2", "g5g0g3" };
  for (int c = 0; c < 16; c++) {
    suNf_propagator vertex;
    int ch = renorm_channel_index(local_names[c]);
    error(ch < 0, 1, "check_tree_level_z", "unknown channel name");
    propagator_from_storage(dat, ch, 0, &vertex);
    if (!check_local_channel(&Sin_inv, &Sout_inv, &vertex, local_names[c], zq, tol)) ok = 0;
  }

  return ok;
}

static int check_gauge_covariance(spinor_field *source, spinor_field *prop_in, spinor_field *prop_out, int nm, double *m,
                                  double p[4], double tol) {
    /* random links: not gauge-equivalent to U=1, so U != U^dag matters */
    random_u(u_gauge);
    represent_gauge_field();
    lprintf("CHECK_GT", 0, "random configuration, plaquette %1.6f\n", avr_plaquette());

    data_storage_array *outA = measure_all_channels(source, prop_in, prop_out, nm, m, p, 0,"GT check A");

    double plaq_before = avr_plaquette();
    random_gauge_transform(u_gauge); /* U -> g U g^dag, g left in gaugefix.c's global g */
    represent_gauge_field();
    lprintf("CHECK_GT", 0, "after random gauge transformation, plaquette %1.6f (must be unchanged)\n",
            avr_plaquette());
    if (fabs(avr_plaquette() - plaq_before) > 1e-12) {
        lprintf("CHECK_GT", 0, "plaquette changed: the gauge transformation is broken\n");
        free_data_storage(outA);
        return 0;
    }

    data_storage_array *outB = measure_all_channels(source, prop_in, prop_out, nm, m, p, 1,"GT check B");


        static const char *names[] = { "g5",   "id",     "g0",     "g1",     "g2",     "g3",     "g0g1", "g0g2",
                                   "g0g3", "g5g0",   "g5g1",   "g5g2",   "g5g3",   "g5g0g1", "g5g0g2",
                                   "g5g0g3", "cg0",  "cg1",    "cg2",    "cg3" }; /* not Sin/Sout: gauge variant */
    int ok = 1;
    for (int c = 0; c < (int)(sizeof(names) / sizeof(names[0])); c++) {
        int ch = renorm_channel_index(names[c]);
        error(ch < 0, 1, "check_gauge_covariance", "unknown channel name");
        suNf_propagator A, B;
        propagator_from_storage(outA, ch, 0, &A);
        propagator_from_storage(outB, ch, 0, &B);
        double diff = 0., norm = 0.;
        for (int r = 0; r < PROP_N; r++) {
            for (int s = 0; s < PROP_N; s++) {
                diff = fmax(diff, cabs(_PROP_IDX(A, r, s) - _PROP_IDX(B, r, s)));
                norm = fmax(norm, cabs(_PROP_IDX(A, r, s)));
            }
        }
        int pass = (diff <= tol * norm);
        lprintf("CHECK_GT", 0, "%-7s max|A-B|/max|A| = %.3e  %s\n", names[c], diff / norm, pass ? "ok" : "FAILED");
        if (!pass) ok = 0;
    }
    free_data_storage(outB);
    free_data_storage(outA);
    return ok;
}




int main(int argc, char *argv[])
{
  int return_value=0;
  int nm;
  int k;    /* used as the colour index in the loops further down --
             * undeclared in the originally uploaded file, a pre-existing
             * bug unrelated to the STORE changes, just now surfaced by a
             * stricter compiler default. (The Mass[%d] print used to read
             * k here too, uninitialised at that point -- fixed to print
             * mass index 0 directly, since nm is always 1 in this program.) */
  double m[256];

  spinor_field *source;
  spinor_field *prop_in;
  spinor_field *prop_out;

  /* setup process communications */
  setup_process(&argc, &argv);

  const char *t[] = { "id", "g5", "g0", "g5g0", "g0g1", "g5g0g1" };
    for (int k = 0; k < 6; k++) {
      suNf_propagator G, G2; hr_complex tr;
      build_gamma(&G, t[k]);
      _propagator_mul(G2, G, G);
      _propagator_trace(tr, G2);
      lprintf("CHECK_G", 0, "%-7s Tr[G G]/(4NF) = %+.3f\n", t[k], creal(tr) / (4.0 * NF));
    }

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


  nm = 1;
  m[0] = -atof(mes_var.mstring);
  lprintf("MAIN", 0, "Inverter precision = %e\n", mes_var.precision);
  lprintf("MAIN", 0, "Mass[%d] = %f\n", 0, m[0]);

  lprintf("MAIN", 0, "Number of maximum monentum component\n", mes_var.n_mom);
  lprintf("MAIN", 0, "Number of twists per momentum\n", mes_var.n_twist);
  if (mes_var.ne)
    lprintf("MAIN", 0, "Doing non-exceptional configuration\n");
  lprintf("MAIN", 0, "Momentum 1 (%d, %d, %d, %d)\n", mes_var.pt_out, mes_var.px_out, mes_var.py_out, mes_var.pz_out);
  if (mes_var.ne)
    lprintf("MAIN", 0, "Momentum 2 (%d, %d, %d, %d)\n", mes_var.pt_in, mes_var.px_in, mes_var.py_in, mes_var.pz_in);

  source = alloc_spinor_field(4, &glattice);
  prop_in = alloc_spinor_field(4 * nm * NF, &glattice);
  prop_out = alloc_spinor_field(4 * nm * NF, &glattice);

  suNf_field *u_gauge_old_f = alloc_suNf_field(&glattice);
  suNg_field *u_gauge_old = alloc_suNg_field(&glattice);

  struct timeval start, end, etime;

  lprintf("MAIN", 0, "Generating a unit configuration\n");

  unit_gauge(u_gauge);
  copy_suNg_field(u_gauge_old, u_gauge);

  represent_gauge_field();

  lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());
  full_plaquette();

  //Fix Gauge
  double p2 = calc_plaq(u_gauge);
  lprintf("MAIN", 0, "initial plaq %1.6f\n", p2);

  gettimeofday(&start, 0);
  double act = gaugefix(10,     //= 0, 1, 2, 3 for Coulomb guage else Landau
                        1.8,    //overrelax
                        10000,  //maxit
                        1e-10,  //tolerance
                        u_gauge //gauge
  );
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("MAIN", 0, "Unit configuration Gauge Fixed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);
  lprintf("MAIN", 0, "action  %1.6f\n", act);
  p2 = calc_plaq(u_gauge);
  lprintf("MAIN", 0, "fixed gauge plaq %1.6f\n", p2);

  represent_gauge_field();
  gettimeofday(&start, 0);

  copy_suNf_field(u_gauge_old_f, u_gauge_f);

  init_propagator_eo(nm, m, mes_var.precision);

  int l, j, tw, num;
  num = 0;

  /* Looked up once by name (see renorm_channel_index in measure_renormalization.c
   * -- the channel list itself stays private to that file) rather than
   * hardcoded indices. */
  int idx_Sin = renorm_channel_index("Sin");
  int idx_Sout = renorm_channel_index("Sout");
  int idx_id = renorm_channel_index("id");
  int idx_g5 = renorm_channel_index("g5");
  if (idx_Sin < 0 || idx_Sout < 0 || idx_id < 0 || idx_g5 < 0)
  {
    error(1, 1, "check_RIMOM " __FILE__, "renorm_channel_index couldn't find one of Sin/Sout/id/g5 -- "
                                        "channel name mismatch with measure_renormalization.c?");
  }

  double mom_in[4], mom_out[4];
  mom_out[0] = mes_var.pt_out;
  mom_out[1] = mes_var.px_out;
  mom_out[2] = mes_var.py_out;
  mom_out[3] = mes_var.pz_out;
  mom_in[0] = mes_var.pt_in;
  mom_in[1] = mes_var.px_in;
  mom_in[2] = mes_var.py_in;
  mom_in[3] = mes_var.pz_in;

  for (l = 1; l <= mes_var.n_mom; ++l)
  {
    double p_in[4], p_out[4];
    for (tw = -mes_var.n_twist; tw < mes_var.n_twist + 1; tw++)
    {
      restore_untwisted_links(u_gauge_old_f);
      double twist = (double)tw * 0.5; //(double)tw/(double)(mes_var.n_twist+1);
      lprintf("TEST", 0, "<p> before twist %1.6f\n", avr_plaquette());
      twist_XYZ_bc(twist * mes_var.px_in, twist * mes_var.py_in, twist * mes_var.pz_in);
      lprintf("TEST", 0, "<p> after twist %1.6f\n", avr_plaquette());

      p_in[0] = mom_in[0] * l;
      p_in[1] = mom_in[1] * l;
      p_in[2] = mom_in[2] * l;
      p_in[3] = mom_in[3] * l;

      for (k = 0; k < NF; ++k)
      {
        create_gauge_fixed_momentum_source(source, p_in[0], p_in[1], p_in[2], p_in[3], k);
        calc_propagator(prop_in + 4 * k, source, 4); //4 for spin components
      }

      if (mes_var.ne)
      {
        restore_untwisted_links(u_gauge_old_f);
        twist_XYZ_bc(twist * mes_var.px_out, twist * mes_var.py_out, twist * mes_var.pz_out);

        p_out[0] = mom_out[0] * l;
        p_out[1] = mom_out[1] * l;
        p_out[2] = mom_out[2] * l;
        p_out[3] = mom_out[3] * l;
        for (k = 0; k < NF; ++k)
        {
          create_gauge_fixed_momentum_source(source, p_out[0], p_out[1], p_out[2], p_out[3], k);
          calc_propagator(prop_out + 4 * k, source, 4); //4 for spin components
        }
      }
      else
      {
        p_out[0] = p_in[0];
        p_out[1] = p_in[1];
        p_out[2] = p_in[2];
        p_out[3] = p_in[3];
        for (j = 0; j < 4 * NF; j++)
          copy_spinor_field(&prop_out[j], &prop_in[j]);
      }
      lprintf("LOOK", 10, "%g%g%g%g %g%g%g%g twist %g", p_in[0], p_in[1], p_in[2], p_in[3], p_out[0], p_out[1], p_out[2], p_out[3], twist);
      char label[256];

      /* Tree-level self-check (see definitions above main()). Only valid
       * for exceptional kinematics (p_in == p_out) and needs twist == 0
       * so p != 0 without the twist term's cancellation risk. Runs once
       * per mom_idx at twist 0; skipped for mes_var.ne (non-exceptional)
       * input files -- in that case we still need the regular
       * (DONTSTORE) measure_renormalization() call below for the normal
       * print_renormalization() output. */
      // if (!mes_var.ne && tw == 0)
      if (!mes_var.ne)
      {
        data_storage_array *out_corr = NULL;
        suNf_propagator Sin_acc, Sout_acc, id_acc, g5_acc;
        int zcheck;
       
        /* Momentum seen by the Dirac operator: twist_XYZ_bc multiplies U_mu by
         * exp(i*pi*theta_mu/L_mu), theta_mu = twist * p_mu^in (T untwisted),
         * i.e. a shift of theta_mu/2 in units of 2*pi/L_mu. */
        double p_eff[4];
        p_eff[0] = p_in[0];
        p_eff[1] = p_in[1] + 0.5 * twist * mes_var.px_in;
        p_eff[2] = p_in[2] + 0.5 * twist * mes_var.py_in;
        p_eff[3] = p_in[3] + 0.5 * twist * mes_var.pz_in;
        measure_renormalization(prop_in, prop_out, nm, p_in[0], p_in[1], p_in[2], p_in[3], p_out[0], p_out[1], p_out[2], p_out[3],
                                STORE, &out_corr);

        propagator_from_storage(out_corr, idx_Sin, 0, &Sin_acc);
        propagator_from_storage(out_corr, idx_Sout, 0, &Sout_acc);
        propagator_from_storage(out_corr, idx_id, 0, &id_acc);
        propagator_from_storage(out_corr, idx_g5, 0, &g5_acc);


        zcheck = check_tree_level_z(&Sin_acc, &Sout_acc, &id_acc, &g5_acc,out_corr,p_eff, 1e-6);
        
        free_data_storage(out_corr);
        lprintf("CHECK_Z", 0, "mom_idx %d twist %d: tree-level Z check %s\n", num, tw,
                zcheck ? "PASSED" : "FAILED");
        if (!zcheck) return_value = 1;
      }
      else
      {
        measure_renormalization(prop_in, prop_out, nm, p_in[0], p_in[1], p_in[2], p_in[3], p_out[0], p_out[1], p_out[2], p_out[3],
                                DONTSTORE, NULL);
        if (mes_var.ne && l == 1 && tw == -mes_var.n_twist)
        {
          lprintf("CHECK_Z", 0, "mes:ne=1 (non-exceptional kinematics): tree-level self-check "
                              "skipped for this input file -- use a dedicated mes:ne=0 input for CI\n");
        }
      }

      sprintf(label, "NPR mom_idx %d twist %d ", num, tw);
      print_renormalization(0, nm, m, label, p_in[0], p_in[1], p_in[2], p_in[3], p_out[0], p_out[1], p_out[2], p_out[3]);

      num++;
    }
  }
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("MAIN", 0, "Unit configuration analysed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

  // The analytical (tree-level) check is now done inline above, per momentum
  // point, via check_tree_level_z() -- reading the 3-point functions
  // measure_renormalization(..., STORE, &out_corr) fills directly.
 {
    double p_gt[4] = { mom_in[0], mom_in[1], mom_in[2], mom_in[3] };
    if (!check_gauge_covariance(source, prop_in, prop_out, nm,m, p_gt, 1e-8)) return_value = 1;
  }


  global_sum_int(&return_value,1);

  if (return_value == 0)
    lprintf("MAIN", 0, "check_RIMOM: PASS (all tree-level Z factors = 1 to tolerance)\n");
  else
    lprintf("MAIN", 0, "check_RIMOM: FAIL (see CHECK_Z/CHECK GT lines above for which quantity and momentum)\n");

  free_spinor_field(source);
  free_spinor_field(prop_in);
  free_spinor_field(prop_out);

  free_propagator_eo();

  free_BCs();

  free_suNg_field(u_gauge);
  free_suNg_field(u_gauge_old);
  free_suNf_field(u_gauge_old_f);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_suNf_field(u_gauge_f);
#endif

  finalize_process();
  lprintf("MAIN", 0, "return_value= %d\n ",  return_value);

  return return_value;
}
