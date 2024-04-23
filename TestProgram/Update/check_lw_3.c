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

    double beta = 5.8;
    double c0 = 1.6667;
    double c1 = (1. - c0) / 8.;

    scalar_field *loc_action = alloc_scalar_field(1, &glattice);
    scalar_field *tmp = alloc_scalar_field(1, &glattice);

    suNg_field *gfield_tmp = alloc_suNg_field(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
    suNf_field *gfield_f_tmp = alloc_suNf_field(&glattice);
#endif

    random_scalar_field_cpu(loc_action);
    copy_to_gpu_scalar_field(loc_action);

    lw_local_action_gpu(loc_action, beta, c0, c1);
    lw_local_action_cpu(loc_action, beta, c0, c1);

    cudaMemcpy(tmp->gpu_ptr, loc_action->gpu_ptr, sizeof(double) * glattice.gsize_spinor, cudaMemcpyDeviceToDevice);
    copy_from_gpu_scalar_field(tmp);

    lprintf("SANITY", 0, "cpu loc action sqnorm: %0.2e\n", sqnorm_cpu(loc_action));
    lprintf("SANITY", 0, "gpu loc action sqnorm: %0.2e\n", sqnorm_cpu(tmp));

    lprintf("TEST", 0, "Checking local action\n");
    sub_assign_cpu(loc_action, tmp);
    double sqnorm = sqrt(sqnorm_cpu(loc_action));
    return_val += check_diff_norm(sqnorm, 1e-12);

    lprintf("TEST", 0, "Testing gauge field\n");
    cudaMemcpy(gfield_tmp->gpu_ptr, u_gauge->gpu_ptr, 4 * sizeof(suNg) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNg_field(gfield_tmp);
    sub_assign_cpu(u_gauge, gfield_tmp);
    sqnorm = sqnorm_cpu(u_gauge);
    return_val += check_diff_norm_zero(sqnorm);

#ifdef ALLOCATE_REPR_GAUGE_FIELD
    lprintf("TEST", 0, "Testing represented gauge field\n");
    cudaMemcpy(gfield_f_tmp->gpu_ptr, u_gauge_f->gpu_ptr, 4 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNf_field(gfield_f_tmp);
    sub_assign_cpu(u_gauge_f, gfield_f_tmp);
    sqnorm = sqnorm_cpu(u_gauge_f);
    return_val += check_diff_norm_zero(sqnorm);
#endif

    finalize_process();

    return return_val;
}