/*
 * NOCOMPILE= BC_T_SF
 * NOCOMPILE= BC_T_SF_ROTATED
 * NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
 */
/*******************************************************************************
 *
 * Checks of propagator, spinmatrix and the sequential sources
 *
 *******************************************************************************/

#include "libhr.h"

#if defined(BC_T_SF_ROTATED) && defined(BC_T_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

/* Mesons parameters */
typedef struct input_mesons {
    char mstring[256];
    double precision;
    int ti;
    int tf;
    int ff_fixed_point;
    int dt;
    /* for the reading function */
    input_record_t read[8];
} input_mesons;

input_mesons mes_var;

// Check: g5 D^dag(x,0) g5 = D(0,x)
static int check_g5herm(spinor_field *prop1, int t1, spinor_field *prop2) {
    int beta, a, ix1, ix2;

    suNf_propagator sp1, sp2, spdag;
    hr_complex tr = 0.;
#ifdef WITH_MPI
    MPI_Status st;
    MPI_Request req;
#endif

    int sender_rank = 0;
    if (t1 >= zerocoord[0] && t1 < zerocoord[0] + T && COORD[1] == 0 && COORD[2] == 0 && COORD[3] == 0) {
        sender_rank = PID;
        ix1 = ipt(t1 - zerocoord[0], 0, 0, 0);
        for (a = 0; a < NF; ++a) {
            for (beta = 0; beta < 4; beta++) {
                _propagator_assign(sp1, *_FIELD_AT(&prop1[a * 4 + beta], ix1), a, beta); // S( (ti,0,0,0), (0,0,0,0) )
            }
        }

#ifdef WITH_MPI
        MPI_Isend(&sp1, NF * 4 * 4 * sizeof(suNf_vector) / sizeof(double), MPI_DOUBLE, 0, 999, GLB_COMM, &req);
#endif
    }

    global_sum_int(&sender_rank, 1);

    if (COORD[0] == 0 && COORD[1] == 0 && COORD[2] == 0 && COORD[3] == 0) {
#ifdef WITH_MPI
        MPI_Irecv(&sp1, NF * 4 * 4 * sizeof(suNf_vector) / sizeof(double), MPI_DOUBLE, sender_rank, 999, GLB_COMM, &req);
        MPI_Wait(&req, &st);
#endif
        ix2 = ipt(0, 0, 0, 0);

        for (a = 0; a < NF; ++a) {
            for (beta = 0; beta < 4; beta++) {
                _propagator_assign(sp2, *_FIELD_AT(&prop2[a * 4 + beta], ix2), a, beta); // S( (0,0,0,0), (ti,0,0) )
            }
        }
        _propagator_dagger(spdag, sp2);
        _g5_propagator(sp2, spdag);
        _propagator_g5(spdag, sp2);
        lprintf("CK_G5HERM", 0, "Propagator1\n");
        // print_prop(sp1);
        lprintf("CK_G5HERM", 0, "g5 Propagator2^dagger g5 \n");
        // print_prop(spdag);
        _propagator_sub(sp2, sp1, spdag);
        lprintf("CK_G5HERM", 0, "Propagator1 - g5 Propagator2^dagger g5 \n");
        // print_prop(sp2);
        _propagator_trace(tr, sp2);
        lprintf("CK_G5HERM", 0, "Tr[ g5 Propagator1^dag g5 - Propagator2 ] = %g + I%g\n", creal(tr), cimag(tr));
    }

    global_sum((double *)&tr, 2);

    if (cabs(tr) > 1.e-14) {
        return 1;
    } else {
        return 0;
    }
}

// source = g5 prop( x, 0 ) delta( x, (tf,0,0,0) )
void create_sequential_source_point(spinor_field *source, int tf, spinor_field *prop) {
    int beta, a, ix;

    suNf_propagator sp0, sp1;

    for (beta = 0; beta < 4 * NF; ++beta) {
#ifdef WITH_GPU
        zero_spinor_field_cpu(&source[beta]);
#endif
        zero_spinor_field(&source[beta]);
    }

    ix = ipt(tf - zerocoord[0], 0, 0, 0);
    for (a = 0; a < NF; ++a) {
        for (beta = 0; beta < 4; beta++) {
            _propagator_assign(sp0, *_FIELD_AT(&prop[a * 4 + beta], ix), a, beta);
        }
    }
    _g5_propagator(sp1, sp0); // g5 Prop
    _propagator_transpose(sp0, sp1);

    for (a = 0; a < NF; ++a) {
        for (beta = 0; beta < 4; beta++) {
            *_FIELD_AT(&source[a * 4 + beta], ix) = sp0.c[a].c[beta];
        }
    }

    for (beta = 0; beta < 4 * NF; ++beta) {
        start_sendrecv_spinor_field(source + beta);
        complete_sendrecv_spinor_field(source + beta);
    }
}

int main(int argc, char *argv[]) {
    int k, tmp;
    int return_value = 0;
    double m[256];
    struct timeval start, end, etime;
    /* setup process id and communications */

    logger_map("DEBUG", "debug");

    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary
    setup_process(&argc, &argv);

    setup_gauge_fields();

    spinor_field *source;
    spinor_field *prop_1;
    spinor_field *prop_2;
    spinor_field *source_seq;
    spinor_field *prop_seq;

    source = alloc_spinor_field(4, &glattice);
    source_seq = alloc_spinor_field(4 * NF, &glattice);
    prop_1 = alloc_spinor_field(4 * NF, &glattice);
    prop_2 = alloc_spinor_field(4 * NF, &glattice);
    prop_seq = alloc_spinor_field(4 * NF, &glattice);

    for (k = 0; k < 4 * NF; k++) {
        zero_spinor_field(prop_1 + k);
        zero_spinor_field(prop_2 + k);
        zero_spinor_field(prop_seq + k);
    }

    mes_var.precision = 1e-24;
    mes_var.ti = GLB_T / 4;
    mes_var.tf = 3 * GLB_T / 4;
    m[0] = 20.12;

    lprintf("MAIN", 0, "ti = %d\n", mes_var.ti);
    lprintf("MAIN", 0, "tf = %d\n", mes_var.tf);
    lprintf("MAIN", 0, "m = %f\n", m[0]);
    mes_var.ff_fixed_point = 1;
    lprintf("MAIN", 0, "Inverter precision = %e\n", mes_var.precision);

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);
    represent_gauge_field();

    gettimeofday(&start, 0);

    init_propagator_eo(1, m, mes_var.precision); // 1 for number of masses
    for (k = 0; k < NF; ++k) {
        create_point_source(source, 0, k);
        calc_propagator(prop_1 + 4 * k, source, 4); // 4 for spin components
        create_point_source(source, mes_var.ti, k);
        calc_propagator(prop_2 + 4 * k, source, 4); // 4 for spin components

#ifdef WITH_GPU
        for (int beta = 0; beta < 4; beta++) {
            copy_from_gpu(prop_1 + beta);
            copy_from_gpu(prop_2 + beta);
        }
#endif
    }
    tmp = check_g5herm(prop_1, mes_var.ti, prop_2);
    return_value += tmp;
    //  create_sequential_source_point(source_seq, mes_var.ti, prop_1);
    //  calc_propagator(prop_seq, source_seq, 4 * NF);
    //  tmp = check_sequential_point(prop_1, prop_2, prop_seq, mes_var.ti);
    //  return_value += tmp;
    //  create_sequential_source(source_seq, mes_var.tf, prop_1);
    //  calc_propagator(prop_seq, source_seq, 4 * NF);
    //  tmp = check_sequential(prop_seq, prop_1, 0, mes_var.tf);
    //  return_value += tmp;
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MAIN", 0, "Random configuration analysed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

    free_propagator_eo();
    free_spinor_field(source);
    free_spinor_field(source_seq);
    free_spinor_field(prop_1);
    free_spinor_field(prop_2);
    free_spinor_field(prop_seq);

    lprintf("MAIN", 0, "return_value = %d\n", return_value);
    finalize_process();

    return return_value;
}
