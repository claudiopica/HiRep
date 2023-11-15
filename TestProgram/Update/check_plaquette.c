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
    start_sendrecv_gfield(u_gauge);
    double plaq_cpu = avr_plaquette_cpu();
    start_sendrecv_gfield(u_gauge);
    double plaq_gpu = avr_plaquette_gpu();
    compare_diff(errors, plaq_cpu, plaq_gpu);

    /* full plaquette */
    // Compare output by eye
    full_plaquette_gpu();
    full_plaquette_cpu();

    double plaqt_cpu[GLB_T];
    double plaqt_gpu[GLB_T];
    double plaqs_cpu[GLB_T];
    double plaqs_gpu[GLB_T];

    avr_plaquette_time_cpu(plaqt_cpu, plaqs_cpu);
    avr_plaquette_time_gpu(plaqt_gpu, plaqs_gpu);

    for (int nt = 0; nt < GLB_T; nt++) {
        compare_diff(errors, plaqt_cpu[nt], plaqt_gpu[nt]);
        compare_diff(errors, plaqs_cpu[nt], plaqs_gpu[nt]);
    }

    finalize_process();
    return errors;
}
