/*******************************************************************************
  *
  * Test gauge fixing. 
  * Original : R. Arthur ?
  * 
  *******************************************************************************/

#include "libhr.h"

static double calc_plaq_diff(suNg_field *V, suNg_field *W) {
#ifdef WITH_GPU
    copy_from_gpu(V);
    V->comm_type = CPU_COMM;
    start_sendrecv_suNg_field(V);
    complete_sendrecv_suNg_field(V);
    copy_from_gpu(W);
    W->comm_type = CPU_COMM;
    start_sendrecv_suNg_field(W);
    complete_sendrecv_suNg_field(W);
#endif
    int mu, nu;
    int iy, iz;
    suNg *v1, *v2, *v3, *v4, w1, w2, w3;
    double pl, E = 0;

    int t, x, y, z, ix;
    for (t = 0; t < T; t++) {
        for (x = 0; x < X; x++) {
            for (y = 0; y < Y; y++) {
                for (z = 0; z < Z; z++) {
                    ix = ipt(t, x, y, z);
                    for (mu = 0; mu < 4; mu++) {
                        for (nu = mu + 1; nu < 4; nu++) {
                            iy = iup(ix, mu);
                            iz = iup(ix, nu);

                            v1 = _4FIELD_AT(V, ix, mu);
                            v2 = _4FIELD_AT(V, iy, nu);
                            v3 = _4FIELD_AT(V, iz, mu);
                            v4 = _4FIELD_AT(V, ix, nu);

                            _suNg_times_suNg(w1, (*v1), (*v2));
                            _suNg_times_suNg(w2, (*v4), (*v3));
                            _suNg_times_suNg_dagger(w3, w1, w2);

                            _suNg_trace_re(pl, w3);
                            E += pl;

                            v1 = _4FIELD_AT(W, ix, mu);
                            v2 = _4FIELD_AT(W, iy, nu);
                            v3 = _4FIELD_AT(W, iz, mu);
                            v4 = _4FIELD_AT(W, ix, nu);

                            _suNg_times_suNg(w1, (*v1), (*v2));
                            _suNg_times_suNg(w2, (*v4), (*v3));
                            _suNg_times_suNg_dagger(w3, w1, w2);

                            _suNg_trace_re(pl, w3);
                            E -= pl;

#ifdef PLAQ_WEIGHTS
                            if (plaq_weight != NULL) { pl *= plaq_weight[ix * 16 + mu * 4 + nu]; }
#endif
                        }
                    }
                }
            }
        }
    }

#ifdef WITH_GPU
    V->comm_type = GPU_COMM;
    W->comm_type = GPU_COMM;
#endif

    global_sum(&E, 1);
    return E / (6. * NG) / GLB_VOLUME;
}

static void random_g(gtransf *g) {
    _MASTER_FOR(g->type, ix) {
        random_suNg(_FIELD_AT(g, ix));
    }
#ifdef WITH_GPU
    copy_to_gpu(g);
#endif
}

static void transform_u(suNg_field *out, suNg_field *in, gtransf *g) {
    int iy, mu;
    suNg v;

#ifdef WITH_GPU
    copy_from_gpu(g);
    g->comm_type = CPU_COMM;
    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);
    copy_from_gpu(in);
    in->comm_type = CPU_COMM;
    start_sendrecv_suNg_field(in);
    complete_sendrecv_suNg_field(in);
    copy_from_gpu(out);
    out->comm_type = CPU_COMM;
    start_sendrecv_suNg_field(out);
    complete_sendrecv_suNg_field(out);
#endif

    _MASTER_FOR(&glattice, ix) {
        for (mu = 0; mu < 4; mu++) {
            iy = iup(ix, mu);
            _suNg_times_suNg_dagger(v, *_4FIELD_AT(in, ix, mu), *_FIELD_AT(g, iy));
            _suNg_times_suNg(*_4FIELD_AT(out, ix, mu), *_FIELD_AT(g, ix), v);
        }
    }

#ifdef WITH_GPU
    copy_to_gpu(g);
    g->comm_type = GPU_COMM;
    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);
    copy_to_gpu(in);
    in->comm_type = GPU_COMM;
    start_sendrecv_suNg_field(in);
    complete_sendrecv_suNg_field(in);
    copy_to_gpu(out);
    out->comm_type = GPU_COMM;
    start_sendrecv_suNg_field(out);
    complete_sendrecv_suNg_field(out);
#endif
}

int main(int argc, char *argv[]) {
    int return_value = 0;
    double p1, p2;
    double pdiff;
    double act;
    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary

    /* setup process communications */
    setup_process(&argc, &argv);
    setup_gauge_fields();
    gtransf *g = alloc_gtransf(&glattice);
    suNg_field *fixed_gauge = alloc_suNg_field(&glattice);

    // initialise random gauge field
    lprintf("TEST", 0, "Perform test gauge invariance of the gauge fixing with a random gauge field\n");
    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);

    // the gauge transformation
    random_g(g);
    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);

    p1 = avr_plaquette();
    lprintf("TEST", 0, "original gauge plaq %1.14f\n", p1);

    transform_u(fixed_gauge, u_gauge, g);
    start_sendrecv_suNg_field(fixed_gauge);
    complete_sendrecv_suNg_field(fixed_gauge);

    p2 = calc_plaq(fixed_gauge);
    lprintf("TEST", 0, "plaq after random gauge transformation %1.14f\n", p2);
    if (fabs(p1 - p2) > 1e-14) { return_value += 1; }

    pdiff = calc_plaq_diff(u_gauge, fixed_gauge);
    if (fabs(pdiff) > 1e-14) { return_value += 1; }
    lprintf("TEST", 0, "sum of difference of plaquettes original/gauge transformed field %1.14f\n", pdiff);

    act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb guage else Landau
                   1.8, //overrelax
                   10000, //maxit
                   1e-10, //tolerance
                   fixed_gauge //gauge
    );
    lprintf("TEST", 0, "action  %1.14f\n", act);

    p2 = calc_plaq(fixed_gauge);
    lprintf("TEST", 0, "plaq after random gauge transformation and gauge fixing %1.14f\n", p2);
    if (fabs(p1 - p2) > 1e-14) { return_value += 1; }

    pdiff = calc_plaq_diff(u_gauge, fixed_gauge);
    if (fabs(pdiff) > 1e-14) { return_value += 1; }
    lprintf("TEST", 0, "sum of difference of plaquettes original/gauge transformed and gauge fixed field %1.14f\n", pdiff);

    lprintf("TEST", 0, "Perform test gauge invariance of the gauge fixing with a unit gauge field\n");
    // initialise unit gauge field
    unit_gauge(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);

    // the gauge transformation
    random_g(g);
    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);

    p1 = calc_plaq(u_gauge);
    lprintf("TEST", 0, "original gauge plaq %1.14f\n", p1);

    transform_u(fixed_gauge, u_gauge, g);
    start_sendrecv_suNg_field(fixed_gauge);
    complete_sendrecv_suNg_field(fixed_gauge);

    p2 = calc_plaq(fixed_gauge);
    lprintf("TEST", 0, "plaq after random gauge transformation %1.14f\n", p2);
    if (fabs(p1 - p2) > 1e-14) { return_value += 1; }

    pdiff = calc_plaq_diff(u_gauge, fixed_gauge);
    if (fabs(pdiff) > 1e-14) { return_value += 1; }
    lprintf("TEST", 0, "sum of difference of plaquettes original/gauge transformed field %1.14f\n", pdiff);

    act = gaugefix(0, //= 0, 1, 2, 3 for Coulomb gauge else Landau
                   1.8, //overrelax
                   10000, //maxit
                   1e-10, //tolerance
                   fixed_gauge //gauge
    );
    lprintf("TEST", 0, "action  %1.14f\n", act);

    p2 = calc_plaq(fixed_gauge);
    lprintf("TEST", 0, "plaq after random gauge transformation and gauge fixing %1.14f\n", p2);
    if (fabs(p1 - p2) > 1e-14) { return_value += 1; }

    pdiff = calc_plaq_diff(u_gauge, fixed_gauge);
    if (fabs(pdiff) > 1e-14) { return_value += 1; }
    lprintf("TEST", 0, "sum of difference of plaquettes original/gauge transformed and gauge fixed field %1.14f\n", pdiff);

    global_sum_int(&return_value, 1);
    lprintf("MAIN", 0, "return_value= %d\n ", return_value);

    finalize_process();

    return return_value;
}
