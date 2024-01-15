/*******************************************************************************
* 
* Check GPU plaquette against CPU plaquette
*
* NOCOMPILE= !WITH_GPU
*
*******************************************************************************/

#include "libhr.h"

static int errors = 0;

int main(int argc, char *argv[]) {
    std_comm_t = ALL_COMMS;
    logger_map("DEBUG", "debug");

    setup_process(&argc, &argv);
    setup_random_gauge_fields();

    /* Average plaquette */
    start_sendrecv_suNg_field(u_gauge);
    double plaq_cpu = avr_plaquette_cpu();
    start_sendrecv_suNg_field(u_gauge);
    double plaq_gpu = avr_plaquette_gpu();
    compare_diff(errors, plaq_cpu, plaq_gpu, "CHECK PLAQUETTE", EPSILON_TEST);

    /* full plaquette */
    // Compare output by eye
    full_plaquette_gpu();
    full_plaquette_cpu();

    double plaqt_cpu[GLB_T];
    double plaqt_gpu[GLB_T];
    double plaqs_cpu[GLB_T];
    double plaqs_gpu[GLB_T];

    scalar_field *sgpu = alloc_scalar_field(1, &glattice);
    scalar_field *scpu = alloc_scalar_field(1, &glattice);

    local_plaquette(sgpu);
    _MASTER_FOR(&glattice, ix) {
        *_FIELD_AT(scpu, ix) = local_plaq(ix);
    }

    copy_to_gpu_scalar_field(scpu);
    lprintf("SANITY CHECK", 0, "L2 diff: %0.15e, %0.15e\n", sqnorm(scpu), sqnorm(sgpu));
    sub_assign(scpu, sgpu);
    lprintf("LOCAL PLAQUETTE", 0, "L2 diff: %0.15e\n", sqnorm(scpu));

    avr_plaquette_time_cpu(plaqt_cpu, plaqs_cpu);
    avr_plaquette_time_gpu(plaqt_gpu, plaqs_gpu);

    for (int nt = 0; nt < GLB_T; nt++) {
        compare_diff(errors, plaqt_cpu[nt], plaqt_gpu[nt], "CHECK PLAQUETTE", EPSILON_TEST);
        compare_diff(errors, plaqs_cpu[nt], plaqs_gpu[nt], "CHECK_PLAQUETTE", EPSILON_TEST);
    }

    finalize_process();
    return errors;
}
