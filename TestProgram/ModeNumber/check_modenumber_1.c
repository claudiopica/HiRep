/*******************************************************************************
*
* NOCOMPILE= !WITH_CLOVER && !WITH_EXPCLOVER
* NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
* NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
*
* Check that with csw = 0 the clover (or exponential clover) mode number
* reproduces the Wilson one:
*
*  1. operator: g5Cphi_sq(m) = g5Dphi_sq(m) at csw = 0 (to rounding);
*  2. mode number: ModeNumber() of ModeNumber/modenumber.c (clover operator)
*     at csw = 0 equals the same code compiled with the Wilson operator
*     (modenumber_wilson.c), on identical Z2 sources;
*  3. sanity: at csw = 1 the two operators differ, i.e. the clover path is
*     really the one being tested.
*
* check_modenumber_1.approx.in was generated with
*   ModeNumber/approx_for_modenumber 0.05 0.05
* (minmax polynomial of degree 5, M/M* != 1). A loose approximation keeps the
* test cheap; the comparison does not depend on its accuracy.
*
*******************************************************************************/

#include "libhr.h"
#include "../../ModeNumber/modenumber.h"

void init_modenumber_wilson(double m, double inv, int nh, char *approxfile);
void free_modenumber_wilson(void);
double ModeNumber_wilson(double M2);

static double hmass = 0.1;
static char approx_file[] = "check_modenumber_1.approx.in";

/* set csw and force the clover term to be recomputed: set_csw() alone does not
   mark it stale, only represent_gauge_field() does */
static void use_csw(double csw) {
    set_csw(&csw);
    represent_gauge_field();
}

/* |g5Cphi_sq s - g5Dphi_sq s| / |g5Dphi_sq s| for a gaussian s */
static double operator_reldiff(double csw) {
    spinor_field *s = alloc_spinor_field(3, &glattice);
    use_csw(csw);
    gaussian_spinor_field(&s[0]);
    g5Dphi_sq(hmass, &s[1], &s[0]);
    g5Cphi_sq(hmass, &s[2], &s[0]);
    sub_assign_spinor_field(&s[2], &s[1]);
    double r = sqrt(sqnorm_spinor_field(&s[2]) / sqnorm_spinor_field(&s[1]));
    free_spinor_field(s);
    return r;
}

/* must_be_small = 1: pass if value < tol; 0: pass if value > tol */
static int check(int must_be_small, const char *what, double value, double tol) {
    int ok = must_be_small ? (value < tol) : (value > tol);
    lprintf("RESULT", 0, "%s: %e (required %s %e) %s\n", what, value, must_be_small ? "<" : ">", tol, ok ? "OK" : "FAILED");
    return ok ? 0 : 1;
}

int main(int argc, char *argv[]) {
    int return_value = 0;

    setup_process(&argc, &argv);
    setup_gauge_fields();

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);
    lprintf("MAIN", 0, "done.\n");

#ifdef WITH_EXPCLOVER
    lprintf("MAIN", 0, "Testing the exponential clover operator\n");
#else
    lprintf("MAIN", 0, "Testing the clover operator\n");
#endif

    /* 1. operator level, csw = 0 */
    double d0 = operator_reldiff(0.0);
    return_value += check(1, "csw=0: |g5Cphi_sq - g5Dphi_sq| / |g5Dphi_sq|", d0, 1.e-14);

    /* 3. sanity, csw = 1: the clover term must make a visible difference */
    double d1 = operator_reldiff(1.0);
    return_value += check(0, "csw=1: |g5Cphi_sq - g5Dphi_sq| / |g5Dphi_sq|", d1, 1.e-3);

    /* 2. mode number, csw = 0, same Z2 sources for both */
    use_csw(0.0);
    const double inverr2 = 1.e-24;
    const int nhits = 2;
    init_modenumber(hmass, inverr2, nhits, approx_file);
    init_modenumber_wilson(hmass, inverr2, nhits, approx_file);

    int nstate = rlxd_size();
    int *state = malloc(sizeof(int) * nstate);
    const double M[] = { 0.3, 1.0 };
    for (int k = 0; k < (int)(sizeof(M) / sizeof(M[0])); k++) {
        rlxd_get(state);
        double nu_c = ModeNumber(M[k] * M[k]);
        rlxd_reset(state);
        double nu_w = ModeNumber_wilson(M[k] * M[k]);
        double rel = fabs(nu_c - nu_w) / fabs(nu_w);
        lprintf("MAIN", 0, "M = %g: nu(clover, csw=0) = %.12e  nu(Wilson) = %.12e\n", M[k], nu_c, nu_w);
        char what[128];
        snprintf(what, sizeof(what), "csw=0: mode number M=%g, relative difference", M[k]);
        return_value += check(1, what, rel, 1.e-10);
    }

    /* negative control: at csw = 1 ModeNumber() must differ from Wilson, which
       shows that modenumber.c really applies the clover operator */
    use_csw(1.0);
    rlxd_get(state);
    double nu_c1 = ModeNumber(M[1] * M[1]);
    rlxd_reset(state);
    double nu_w1 = ModeNumber_wilson(M[1] * M[1]);
    lprintf("MAIN", 0, "M = %g: nu(clover, csw=1) = %.12e  nu(Wilson) = %.12e\n", M[1], nu_c1, nu_w1);
    return_value += check(0, "csw=1: mode number relative difference", fabs(nu_c1 - nu_w1) / fabs(nu_w1), 1.e-3);
    free(state);

    free_modenumber();
    free_modenumber_wilson();

    finalize_process();
    return return_value;
}
