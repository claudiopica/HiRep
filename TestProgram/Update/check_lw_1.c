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

    double dt = 0.1;

    suNg_av_field *force = alloc(force, 1, &glattice);
    random_suNg_av_field_cpu(force);
    copy_to_gpu_suNg_av_field(force);

    suNg_av_field *force_tmp = alloc(force_tmp, 1, &glattice);
    suNg_field *gfield_tmp = alloc(gfield_tmp, 1, &glattice);
    suNf_field *gfield_f_tmp = alloc(gfield_f_tmp, 1, &glattice);

    force_gauge_par *par = (force_gauge_par *)malloc(sizeof(force_gauge_par));
    par->momenta = &force;
    par->beta = 5.8;
    par->c0 = 1.6667;
    par->c1 = (1. - par->c0) / 8.;

    lw_force_gpu(dt, par);
    lw_force_cpu(dt, par);

    cudaMemcpy(force_tmp->gpu_ptr, force->gpu_ptr, 4 * sizeof(suNg_algebra_vector) * glattice.gsize_gauge,
               cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNg_av_field(force_tmp);

    cudaMemcpy(gfield_tmp->gpu_ptr, u_gauge->gpu_ptr, 4 * sizeof(suNg) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNg_field(gfield_tmp);

    cudaMemcpy(gfield_f_tmp->gpu_ptr, u_gauge_f->gpu_ptr, 4 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNf_field(gfield_f_tmp);

    lprintf("SANITY", 0, "sqnorm cpu: %0.6e\n", sqnorm(force));
    lprintf("SANITY", 0, "sqnorm cpu: %0.6e\n", sqnorm(force_tmp));

    lprintf("TEST", 0, "Checking force field\n");
    sub_assign_cpu(force, force_tmp);
    double sqnorm = sqrt(sqnorm_cpu(force));
    return_val += check_diff_norm(sqnorm, 1e-13);

    lprintf("TEST", 0, "Checking gauge field\n");
    sub_assign_cpu(u_gauge, gfield_tmp);
    sqnorm = sqrt(sqnorm_cpu(u_gauge));
    return_val += check_diff_norm_zero(sqnorm);

#ifdef ALLOCATE_REPR_GAUGE_FIELD
    lprintf("TEST", 0, "Checking represented gauge field\n");
    sub_assign_cpu(u_gauge_f, gfield_f_tmp);
    sqnorm = sqrt(sqnorm_cpu(u_gauge_f));
    return_val += check_diff_norm_zero(sqnorm);
#endif

    finalize_process();

    return return_val;
}