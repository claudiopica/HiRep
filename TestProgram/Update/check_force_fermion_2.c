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

    double dt = 0.1;
    double residue = 0.1;

    setup_random_gauge_fields();
    suNg_field *gfield_tmp = alloc_suNg_field(&glattice);

#ifdef ALLOCATE_REPR_GAUGE_FIELD
    suNf_field *gfield_f_tmp = alloc_suNf_field(&glattice);
#endif

    spinor_field *X = alloc_spinor_field(1, &glattice);
    spinor_field *X_tmp = alloc_spinor_field(1, &glattice);
    spinor_field *Y = alloc_spinor_field(1, &glattice);
    spinor_field *Y_tmp = alloc_spinor_field(1, &glattice);

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    clover_term *cl_term_tmp = alloc_clover_term(&glattice);
    clover_force *cl_force_tmp = alloc_clover_force(&glattice);
    ldl_field *cl_ldl_tmp = alloc_ldl_field(&glattice);
    random_ldl_field_cpu(cl_ldl);
    copy_to_gpu_ldl_field(cl_ldl);
#endif

    gaussian_spinor_field(X);
    copy_from_gpu_spinor_field(X);
    gaussian_spinor_field(Y);
    copy_from_gpu_spinor_field(Y);

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    random_clover_force_cpu(cl_force);
    copy_to_gpu_clover_force(cl_force);
    copy_from_gpu_clover_term(cl_term);

    cudaMemcpy(cl_term_tmp->gpu_ptr, cl_term->gpu_ptr, 4 * sizeof(suNfc) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_term(cl_term_tmp);

    cudaMemcpy(cl_force_tmp->gpu_ptr, cl_force->gpu_ptr, 6 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_force(cl_force_tmp);

    cudaMemcpy(cl_ldl_tmp->gpu_ptr, cl_ldl->gpu_ptr, sizeof(ldl_t) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_ldl_field(cl_ldl_tmp);

    lprintf("SANITY", 0, "Clover force before begin CPU: %0.2e\n", sqrt(sqnorm_clover_force_cpu(cl_force)));
    lprintf("SANITY", 0, "Clover force before begin GPU: %0.2e\n", sqrt(sqnorm_clover_force_cpu(cl_force_tmp)));

    lprintf("SANITY", 0, "Clover term before begin CPU: %0.2e\n", sqrt(sqnorm_clover_term_cpu(cl_term)));
    lprintf("SANITY", 0, "Clover term before begin GPU: %0.2e\n", sqrt(sqnorm_clover_term_cpu(cl_term_tmp)));

#if 0
    lprintf("SANITY", 0, "Clover ldl before begin CPU: %0.2e\n", sqrt(sqnorm_ldl_field_cpu(cl_ldl)));
    lprintf("SANITY", 0, "Clover ldl before begin GPU: %0.2e\n", sqrt(sqnorm_ldl_field_cpu(cl_ldl_tmp)));
#endif
#endif

    suNg_av_field *force = alloc_suNg_av_field(&glattice);
    suNg_av_field *force_tmp = alloc_suNg_av_field(&glattice);
    random_suNg_av_field_cpu(force);
    copy_to_gpu_suNg_av_field(force);

    start_sendrecv_suNg_av_field(force);
    complete_sendrecv_suNg_av_field(force);

    fermion_force_begin_gpu();
    fermion_force_begin_cpu();

    lprintf("SANITY", 0, "before X gpu norm: %0.2e\n", sqnorm_spinor_field(X));
    lprintf("SANITY", 0, "before X cpu norm: %0.2e\n", sqnorm_spinor_field_cpu(X));

    lprintf("SANITY", 0, "before Y gpu norm: %0.2e\n", sqnorm_spinor_field(Y));
    lprintf("SANITY", 0, "before Y cpu norm: %0.2e\n", sqnorm_spinor_field_cpu(Y));

    force_fermion_core_gpu(X, Y, 1, dt, residue);
    force_fermion_core_cpu(X, Y, 1, dt, residue);

    fermion_force_end_gpu(dt, force);
    fermion_force_end_cpu(dt, force);

    lprintf("SANITY", 0, "X gpu norm: %0.2e\n", sqnorm_spinor_field(X));
    lprintf("SANITY", 0, "X cpu norm: %0.2e\n", sqnorm_spinor_field_cpu(X));

    lprintf("SANITY", 0, "Y gpu norm: %0.2e\n", sqnorm_spinor_field(Y));
    lprintf("SANITY", 0, "Y cpu norm: %0.2e\n", sqnorm_spinor_field_cpu(Y));

    cudaMemcpy(force_tmp->gpu_ptr, force->gpu_ptr, 4 * sizeof(suNg_algebra_vector) * glattice.gsize_gauge,
               cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNg_av_field(force_tmp);

    lprintf("SANITY", 0, "Force GPU: %0.2e\n", sqnorm_suNg_av_field_cpu(force_tmp));
    lprintf("SANITY", 0, "Force CPU: %0.2e\n", sqnorm_suNg_av_field_cpu(force));

    cudaMemcpy(X_tmp->gpu_ptr, X->gpu_ptr, sizeof(suNf_spinor) * glattice.gsize_spinor, cudaMemcpyDeviceToDevice);
    cudaMemcpy(Y_tmp->gpu_ptr, Y->gpu_ptr, sizeof(suNf_spinor) * glattice.gsize_spinor, cudaMemcpyDeviceToDevice);
    copy_from_gpu_spinor_field(X_tmp);
    copy_from_gpu_spinor_field(Y_tmp);

    cudaMemcpy(gfield_tmp->gpu_ptr, u_gauge->gpu_ptr, 4 * sizeof(suNg) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNg_field(gfield_tmp);

#ifdef ALLOCATE_REPR_GAUGE_FIELD
    cudaMemcpy(gfield_f_tmp->gpu_ptr, u_gauge_f->gpu_ptr, 4 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_suNf_field(gfield_f_tmp);
#endif

    double sqnorm = 0;

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    cudaMemcpy(cl_term_tmp->gpu_ptr, cl_term->gpu_ptr, 4 * sizeof(suNfc) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_term(cl_term_tmp);

    cudaMemcpy(cl_force_tmp->gpu_ptr, cl_force->gpu_ptr, 6 * sizeof(suNf) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_clover_force(cl_force_tmp);

    cudaMemcpy(cl_ldl_tmp->gpu_ptr, cl_ldl->gpu_ptr, sizeof(ldl_t) * glattice.gsize_gauge, cudaMemcpyDeviceToDevice);
    copy_from_gpu_ldl_field(cl_ldl_tmp);

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

#if 0
    lprintf("TEST", 0, "Checking clover ldl\n");
    lprintf("SANITY", 0, "Clover ldl CPU: %0.2e\n", sqrt(sqnorm_ldl_field_cpu(cl_ldl)));
    lprintf("SANITY", 0, "Clover ldl GPU: %0.2e\n", sqrt(sqnorm_ldl_field_cpu(cl_ldl_tmp)));
    sub_assign_ldl_field_cpu(cl_ldl, cl_ldl_tmp);
    sqnorm = sqrt(sqnorm_ldl_field_cpu(cl_ldl));
    return_val += check_diff_norm(sqnorm, 1e-12);
#endif
#endif

    lprintf("TEST", 0, "Checking force\n");
    sub_assign_suNg_av_field_cpu(force, force_tmp);
    sqnorm = sqrt(sqnorm_suNg_av_field_cpu(force));
    return_val += check_diff_norm(sqnorm, 1e-12);

    lprintf("TEST", 0, "Checking X\n");
    sub_assign_spinor_field(X, X_tmp);
    sqnorm = sqrt(sqnorm_spinor_field_cpu(X));
    return_val += check_diff_norm(sqnorm, 1e-12);

    lprintf("TEST", 0, "Checking Y\n");
    sub_assign_spinor_field(Y, Y_tmp);
    sqnorm = sqrt(sqnorm_spinor_field_cpu(Y));
    return_val += check_diff_norm(sqnorm, 1e-12);

    lprintf("TEST", 0, "Checking gauge field\n");
    sub_assign_suNg_field_cpu(u_gauge, gfield_tmp);
    sqnorm = sqrt(sqnorm_suNg_field_cpu(u_gauge));
    return_val += check_diff_norm_zero(sqnorm);

#ifdef ALLOCATE_REPR_GAUGE_FIELD
    lprintf("TEST", 0, "Checking represented gauge field\n");
    sub_assign_suNf_field_cpu(u_gauge_f, gfield_f_tmp);
    sqnorm = sqrt(sqnorm_suNf_field_cpu(u_gauge_f));
    return_val += check_diff_norm_zero(sqnorm);
#endif

    finalize_process();

    return return_val;
}