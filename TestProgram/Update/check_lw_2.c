/******************************************************************************
 *
 * NOCOMPILE= !WITH_GPU
 *
******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary

    int return_val = 0;
    setup_process(&argc, &argv);
    setup_gauge_fields();

    random_u(u_gauge);
    represent_gauge_field();
    copy_from_gpu_suNg_field(u_gauge);
    copy_from_gpu_suNf_field(u_gauge_f);

    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);
    start_sendrecv_suNf_field(u_gauge_f);
    complete_sendrecv_suNf_field(u_gauge_f);

    double beta = 5.8;
    double c0 = 1.6667;
    double c1 = (1. - c0) / 8.;

    double res_gpu = lw_action_gpu(beta, c0, c1);
    double res_cpu = lw_action_cpu(beta, c0, c1);

    lprintf("RESULT", 0, "Result on GPU: %0.15e\n", res_gpu);
    lprintf("RESULT", 0, "Result on CPU: %0.15e\n", res_cpu);

    return_val += check_diff_norm(res_gpu - res_cpu, 1e-12);

    finalize_process();

    return return_val;
}