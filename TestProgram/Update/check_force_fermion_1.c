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
    double sqnorm = 0;

    double dt = 0.1;
    double residue = 0.1;

    setup_random_gauge_fields();
    suNg_field *gfield_tmp = alloc_gfield(&glattice);
    suNf_field *gfield_f_tmp = alloc_gfield_f(&glattice);

    spinor_field *X = alloc_spinor_field_f(1, &glattice);
    spinor_field *X_tmp = alloc_spinor_field_f(1, &glattice);
    spinor_field *Y = alloc_spinor_field_f(1, &glattice);
    spinor_field *Y_tmp = alloc_spinor_field_f(1, &glattice);

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    suNfc_field *cl_term_tmp = alloc_clover_term(&glattice);
    suNf_field *cl_force_tmp = alloc_clover_force(&glattice);
    ldl_field *cl_ldl_tmp = alloc_clover_ldl(&glattice);
#endif

    gaussian_spinor_field(X);
    copy_from_gpu_spinor_field_f(X);
    gaussian_spinor_field(Y);
    copy_from_gpu_spinor_field_f(Y);

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    random_clover_force_cpu(cl_force);
    copy_to_gpu_clover_force(cl_force);
    copy_from_gpu_clover_term(cl_term);

    cudaMemcpy(cl_term_tmp->gpu_ptr, cl_term->gpu_ptr, 4 * sizeof(suNfc) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_term(cl_term_tmp);

    cudaMemcpy(cl_force_tmp->gpu_ptr, cl_force->gpu_ptr, 6 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_force(cl_force_tmp);

    cudaMemcpy(cl_ldl_tmp->gpu_ptr, cl_ldl->gpu_ptr, sizeof(ldl_t) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_ldl(cl_ldl_tmp);

#endif

    suNg_av_field *force = alloc_avfield(&glattice);
    suNg_av_field *force_tmp = alloc_avfield(&glattice);
    random_avfield_cpu(force);
    copy_to_gpu_avfield(force);

    start_sendrecv_avfield(force);
    complete_sendrecv_avfield(force);

    fermion_force_begin_gpu();
    fermion_force_begin_cpu();

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    cudaMemcpy(cl_term_tmp->gpu_ptr, cl_term->gpu_ptr, 4 * sizeof(suNfc) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_term(cl_term_tmp);

    cudaMemcpy(cl_force_tmp->gpu_ptr, cl_force->gpu_ptr, 6 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_force(cl_force_tmp);

    cudaMemcpy(cl_ldl_tmp->gpu_ptr, cl_ldl->gpu_ptr, sizeof(ldl_t) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_ldl(cl_ldl_tmp);
#endif

    force_fermion_core_gpu(X, Y, 1, dt, residue);
    force_fermion_core_cpu(X, Y, 1, dt, residue);

    cudaMemcpy(X_tmp->gpu_ptr, X->gpu_ptr, sizeof(suNf_spinor) * glattice.gsize_spinor, cudaMemcpyDeviceToDevice);
    cudaMemcpy(Y_tmp->gpu_ptr, Y->gpu_ptr, sizeof(suNf_spinor) * glattice.gsize_spinor, cudaMemcpyDeviceToDevice);
    copy_from_gpu_spinor_field_f(X_tmp);
    copy_from_gpu_spinor_field_f(Y_tmp);

    cudaMemcpy(gfield_tmp->gpu_ptr, u_gauge->gpu_ptr, 4 * sizeof(suNg) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_gfield(gfield_tmp);

#ifdef ALLOCATE_REPR_GAUGE_FIELD
    cudaMemcpy(gfield_f_tmp->gpu_ptr, u_gauge_f->gpu_ptr, 4 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_gfield_f(gfield_f_tmp);
#endif

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    cudaMemcpy(cl_term_tmp->gpu_ptr, cl_term->gpu_ptr, 4 * sizeof(suNfc) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_term(cl_term_tmp);

    cudaMemcpy(cl_force_tmp->gpu_ptr, cl_force->gpu_ptr, 6 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_force(cl_force_tmp);

    cudaMemcpy(cl_ldl_tmp->gpu_ptr, cl_ldl->gpu_ptr, sizeof(ldl_t) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_ldl(cl_ldl_tmp);

    lprintf("TEST", 0, "Checking clover term\n");
    lprintf("SANITY", 0, "Clover term CPU: %0.2e\n", sqrt(sqnorm_clover_term_cpu(cl_term)));
    lprintf("SANITY", 0, "Clover term GPU: %0.2e\n", sqrt(sqnorm_clover_term_cpu(cl_term_tmp)));
    sub_assign_clover_term_cpu(cl_term, cl_term_tmp);
    sqnorm = sqrt(sqnorm_clover_term_cpu(cl_term));
    return_val += check_diff_norm(sqnorm, 1e-12);

    lprintf("TEST", 0, "Checking clover force\n");
    lprintf("SANITY", 0, "Clover force CPU: %0.2e\n", sqrt(sqnorm_clover_force_cpu(cl_force)));
    lprintf("SANITY", 0, "Clover force GPU: %0.2e\n", sqrt(sqnorm_clover_force_cpu(cl_force_tmp)));
    sub_assign_clover_force_cpu(cl_force, cl_force_tmp);
    sqnorm = sqrt(sqnorm_clover_force_cpu(cl_force));
    return_val += check_diff_norm(sqnorm, 1e-12);

    lprintf("TEST", 0, "Checking clover ldl\n");
    lprintf("SANITY", 0, "Clover ldl CPU: %0.2e\n", sqrt(sqnorm_clover_ldl_cpu(cl_ldl)));
    lprintf("SANITY", 0, "Clover ldl GPU: %0.2e\n", sqrt(sqnorm_clover_ldl_cpu(cl_ldl_tmp)));
    sub_assign_clover_ldl_cpu(cl_ldl, cl_ldl_tmp);
    sqnorm = sqrt(sqnorm_clover_ldl_cpu(cl_ldl));
    return_val += check_diff_norm(sqnorm, 1e-12);
#endif

    lprintf("TEST", 0, "Checking X\n");
    spinor_field_sub_assign_f_cpu(X, X_tmp);
    sqnorm = sqrt(spinor_field_sqnorm_f_cpu(X));
    return_val += check_diff_norm(sqnorm, 1e-12);

    lprintf("TEST", 0, "Checking Y\n");
    spinor_field_sub_assign_f_cpu(Y, Y_tmp);
    sqnorm = sqrt(spinor_field_sqnorm_f_cpu(Y));
    return_val += check_diff_norm(sqnorm, 1e-12);

    lprintf("TEST", 0, "Checking gauge field\n");
    sub_assign_gfield_cpu(u_gauge, gfield_tmp);
    sqnorm = sqrt(sqnorm_gfield_cpu(u_gauge));
    return_val += check_diff_norm_zero(sqnorm);

#ifdef ALLOCATE_REPR_GAUGE_FIELD
    lprintf("TEST", 0, "Checking represented gauge field\n");
    sub_assign_gfield_f_cpu(u_gauge_f, gfield_f_tmp);
    sqnorm = sqrt(sqnorm_gfield_f_cpu(u_gauge_f));
    return_val += check_diff_norm_zero(sqnorm);
#endif

    free_gfield(gfield_tmp);
    free_gfield_f(gfield_f_tmp);
    free_spinor_field_f(X);
    free_spinor_field_f(Y);
    free_spinor_field_f(X_tmp);
    free_spinor_field_f(Y_tmp);
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    free_clover_term(cl_term_tmp);
    free_clover_force(cl_force_tmp);
    free_clover_ldl(cl_ldl_tmp);
#endif
    finalize_process();

    return return_val;
}