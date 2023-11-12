/******************************************************************************
 *
 * NOCOMPILE= !WITH_GPU
 *
******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int return_val = 0;
    setup_process(&argc, &argv);

    // Random gauge
    suNg_field *in_gpu, *in_cpu;
    in_gpu = alloc_gfield(&glattice);
    in_cpu = alloc_gfield(&glattice);

    random_gfield_cpu(in_gpu);
    copy_gfield_cpu(in_cpu, in_gpu);
    copy_to_gpu_gfield(in_gpu);

    // Random force
    suNg_av_field *force_gpu = alloc_avfield(&glattice);
    suNg_av_field *force_cpu = alloc_avfield(&glattice);

    random_avfield_cpu(force_gpu);
    copy_avfield_cpu(force_cpu, force_gpu);
    copy_to_gpu_avfield(force_gpu);

    double dt = 0.1;

    // GPU operation
    exec_field_update(in_gpu, force_gpu, dt);

    // CPU operation
    _MASTER_FOR(&glattice, ix) {
        ExpX(dt, _4FIELD_AT(force_cpu, ix, 0), _4FIELD_AT(in_cpu, ix, 0));
        ExpX(dt, _4FIELD_AT(force_cpu, ix, 1), _4FIELD_AT(in_cpu, ix, 1));
        ExpX(dt, _4FIELD_AT(force_cpu, ix, 2), _4FIELD_AT(in_cpu, ix, 2));
        ExpX(dt, _4FIELD_AT(force_cpu, ix, 3), _4FIELD_AT(in_cpu, ix, 3));
    }

    // Check that both force and gfield are identical.
    copy_from_gpu_gfield(in_gpu);
    copy_from_gpu_avfield(force_gpu);

    lprintf("SANITY CHECK", 0, "in field not zero: %0.15e\n", sqnorm_gfield_cpu(in_gpu));
    lprintf("SANITY CHECK", 0, "cpu not zero: %0.15e\n", sqnorm_gfield_cpu(in_cpu));

    lprintf("SANITY CHECK", 0, "avfield not zero: %0.15e\n", sqnorm_avfield_cpu(force_gpu));
    lprintf("SANITY CHECK", 0, "cpu not zero: %0.15e\n", sqnorm_avfield_cpu(force_cpu));

    sub_assign_gfield_cpu(in_cpu, in_gpu);
    sub_assign_avfield_cpu(force_cpu, force_gpu);

    lprintf("TEST", 0, "Checking Gauge Field\n");
    double diff_norm = sqrt(sqnorm_gfield_cpu(in_cpu));
    return_val += check_diff_norm(diff_norm, 1e-13);

    lprintf("TEST", 0, "Checking force field\n");
    diff_norm = sqrt(sqnorm_avfield_cpu(force_cpu));
    return_val += check_diff_norm_zero(diff_norm);

    finalize_process();

    return return_val;
}