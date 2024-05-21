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
    setup_random_gauge_fields();

    scalar_field *sfield, *tmp;
    sfield = alloc(sfield, 2, &glattice);
    tmp = sfield + 1;
    random_scalar_field_cpu(sfield);
    copy_to_gpu(sfield);

    suNg_av_field *momenta, *momenta_tmp;
    momenta = alloc(momenta, 2, &glattice);
    momenta_tmp = momenta + 1;
    random_suNg_av_field_cpu(momenta);
    copy_to_gpu(momenta);

    local_hmc_action_cpu(NEW, sfield, momenta, NULL);
    local_hmc_action_gpu(NEW, sfield, momenta, NULL);

    lprintf("TEST", 0, "cpu: %0.15e, gpu: %0.15e\n", sqnorm(sfield), sqnorm_cpu(sfield));

    copy(tmp, sfield);
    copy(momenta_tmp, momenta);

    copy_from_gpu(tmp);
    copy_from_gpu(momenta_tmp);

    sub_assign_cpu(sfield, tmp);
    sub_assign_cpu(momenta, momenta_tmp);

    double sdiff = max_cpu(sfield);
    double mdiff = max_cpu(momenta);

    lprintf("TEST", 0, "Checking scalar field\n");
    return_val += check_diff_norm(sdiff, 1e-15);

    lprintf("TEST", 0, "Checking momenta\n");
    return_val += check_diff_norm(mdiff, 1e-15);

    // CHECK DELTA
    random_scalar_field_cpu(sfield);
    copy_to_gpu(sfield);
    random_suNg_av_field_cpu(momenta);
    copy_to_gpu(momenta);

    local_hmc_action_cpu(DELTA, sfield, momenta, NULL);
    local_hmc_action_gpu(DELTA, sfield, momenta, NULL);

    lprintf("TEST", 0, "cpu: %0.15e, gpu: %0.15e\n", sqnorm(sfield), sqnorm_cpu(sfield));

    copy(tmp, sfield);
    copy(momenta_tmp, momenta);

    copy_from_gpu(tmp);
    copy_from_gpu(momenta_tmp);

    sub_assign_cpu(sfield, tmp);
    sub_assign_cpu(momenta, momenta_tmp);

    sdiff = max_cpu(sfield);
    mdiff = max_cpu(momenta);

    lprintf("TEST", 0, "Checking scalar field\n");
    return_val += check_diff_norm(sdiff, 1e-15);

    lprintf("TEST", 0, "Checking momenta\n");
    return_val += check_diff_norm(mdiff, 1e-15);

    // CHECK PF
    random_scalar_field_cpu(sfield);
    copy_to_gpu(sfield);

    spinor_field *pf, *pf_tmp;
    pf = alloc(pf, 2, &glattice);
    pf_tmp = pf + 1;
    gaussian_spinor_field(pf);

    pf_local_action_cpu(sfield, pf);
    pf_local_action_gpu(sfield, pf);

    lprintf("TEST", 0, "cpu: %0.15e, gpu: %0.15e\n", sqnorm(sfield), sqnorm_cpu(sfield));

    copy(tmp, sfield);
    copy(pf_tmp, pf);

    copy_from_gpu(tmp);
    copy_from_gpu(pf_tmp);

    sub_assign_cpu(sfield, tmp);
    sub_assign_cpu(pf, pf_tmp);

    sdiff = max_cpu(sfield);
    mdiff = max_cpu(pf);

    lprintf("TEST", 0, "Checking scalar field\n");
    return_val += check_diff_norm(sdiff, 1e-14);

    lprintf("TEST", 0, "Checking pf\n");
    return_val += check_diff_norm(mdiff, 1e-15);

    finalize_process();

    return return_val;
}