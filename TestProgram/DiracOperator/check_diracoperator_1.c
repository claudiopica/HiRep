/*******************************************************************************
*
* Gauge covariance of the Dirac operator
*
*******************************************************************************/

#include "libhr.h"

static double hmass = 0.1;
static gtransf *g;

static void loc_D(spinor_field *out, spinor_field *in) {
    Dphi(hmass, out, in);
}

static void random_g(void) {
#ifdef WITH_FUSE_MASTER_FOR
    _FUSE_MASTER_FOR(&glattice, ix) {
        _FUSE_IDX(&glattice, ix);
#else
    _MASTER_FOR(&glattice, ix) {
#endif
        random_suNg(_FIELD_AT(g, ix));
    }
#ifdef WITH_GPU
    copy_to_gpu_gtransf(g);
#endif
}

static void transform_u(void) {
#ifdef WITH_GPU
    copy_from_gpu_suNg_field(u_gauge);
    copy_from_gpu_gtransf(g);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);
    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);
#endif
#ifdef WITH_FUSE_MASTER_FOR
    _FUSE_MASTER_FOR(&glattice, ix) {
        _FUSE_IDX(&glattice, ix);
#else
    _MASTER_FOR(&glattice, ix) {
#endif
        suNg v;
        for (int mu = 0; mu < 4; mu++) {
            int iy = iup(ix, mu);
            suNg *u = pu_gauge(ix, mu);
            _suNg_times_suNg_dagger(v, *u, *_FIELD_AT(g, iy));
            _suNg_times_suNg(*u, *_FIELD_AT(g, ix), v);
        }
    }

#ifdef WITH_GPU
    copy_to_gpu_suNg_field(u_gauge);
    copy_to_gpu_gtransf(g);
#endif

    start_sendrecv_suNg_field(u_gauge);
    represent_gauge_field();
    smear_gauge_field();
}

static void transform_s(spinor_field *out, spinor_field *in) {
#ifdef WITH_GPU
    copy_from_gpu_spinor_field(in);
#endif
#ifdef WITH_FUSE_MASTER_FOR
    _FUSE_MASTER_FOR(&glattice, ix) {
        _FUSE_IDX(&glattice, ix);
#else
    _MASTER_FOR(&glattice, ix) {
#endif
        suNf_spinor *s = _FIELD_AT(in, ix);
        suNf_spinor *r = _FIELD_AT(out, ix);
        suNf gfx;

        _group_represent2(&gfx, _FIELD_AT(g, ix));

        _suNf_multiply(r->c[0], gfx, s->c[0]);
        _suNf_multiply(r->c[1], gfx, s->c[1]);
        _suNf_multiply(r->c[2], gfx, s->c[2]);
        _suNf_multiply(r->c[3], gfx, s->c[3]);
    }

#ifdef WITH_GPU
    copy_to_gpu_spinor_field(out);
#endif
}

int main(int argc, char *argv[]) {
    int return_value = 1;
    double sig, tau;
    spinor_field *s0, *s1, *s2, *s3;
#ifdef WITH_GPU
    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary
#endif

    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    setup_gauge_fields();

    /* allocate additional memory */
    g = alloc_gtransf(&glattice);
    s0 = alloc_spinor_field(4, &glattice);
    s1 = s0 + 1;
    s2 = s1 + 1;
    s3 = s2 + 1;

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    fflush(stdout);
    random_u(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n");

    zero_spinor_field(s0);
    gaussian_spinor_field(&(s0[0]));
    tau = 1. / sqrt(sqnorm_spinor_field(s0));
    mul_spinor_field(s0, tau, s0);
    sig = sqnorm_spinor_field(s0);

    lprintf("MAIN", 0, "Normalized norm = %.2e\n", sig);

    lprintf("MAIN", 0, "Generating a random gauge transf... ");
    fflush(stdout);
    random_g();
    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);
    lprintf("MAIN", 0, "done.\n");

    lprintf("MAIN", 0, "Gauge covariance of the Dirac operator:\n");

    loc_D(s1, s0);
    transform_s(s2, s1);
    transform_s(s3, s0);
    transform_u();
    zero_spinor_field(s1);
    loc_D(s1, s3);

    mul_add_assign_spinor_field(s1, -1.0, s2);
    sig = sqnorm_spinor_field(s1);

    lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", sqrt(sig));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");

    if (sqrt(sig) < 10.e-14) { return_value = 0; }

    free_spinor_field(s0);

    free_gtransf(g);

    finalize_process();
    return return_value;
}
