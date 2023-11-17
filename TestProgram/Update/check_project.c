/******************************************************************************
 *
 * NOCOMPILE= !WITH_GPU
 *
******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    std_comm_t = ALL_COMMS;
    int return_val = 0;
    setup_process(&argc, &argv);
    setup_gauge_fields();

    random_u(u_gauge);
    represent_gauge_field();
    copy_from_gpu_suNg_field(u_gauge);
    copy_from_gpu_suNf_field(u_gauge_f);

    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);

    u_gauge_f->comm_type = ALL_COMMS;
    start_sendrecv_suNf_field(u_gauge_f);
    complete_sendrecv_suNf_field(u_gauge_f);

    exec_project();
    _MASTER_FOR(&glattice, ix) {
        project_to_suNg(pu_gauge(ix, 0));
        project_to_suNg(pu_gauge(ix, 1));
        project_to_suNg(pu_gauge(ix, 2));
        project_to_suNg(pu_gauge(ix, 3));
    }

    suNg_field *tmp = alloc_suNg_field(&glattice);
    cudaMemcpy(tmp->gpu_ptr, u_gauge->gpu_ptr, 4 * sizeof(suNg) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNg_field(tmp);

    sub_assign_suNg_field_cpu(u_gauge, tmp);

    double sqnorm = sqrt(sqnorm_suNg_field_cpu(u_gauge));
    return_val += check_diff_norm(sqnorm, 1e-13);

    finalize_process();
    return return_val;
}