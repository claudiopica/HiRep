/******************************************************************************
 *
 * NOCOMPILE= !WITH_GPU
 *
******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int return_val = 0;
    setup_process(&argc, &argv);
    setup_gauge_fields();

    // Random gauge
    suNg_field *in_cpu;
    in_cpu = alloc_suNg_field(&glattice);

    random_suNg_field_cpu(u_gauge);
    copy_to_gpu_suNg_field(u_gauge);

    // Random force
    suNg_av_field *force_gpu = alloc_suNg_av_field(&glattice);
    suNg_av_field *force_cpu = alloc_suNg_av_field(&glattice);

    random_suNg_av_field_cpu(force_gpu);
    copy_suNg_av_field_cpu(force_cpu, force_gpu);
    copy_to_gpu_suNg_av_field(force_gpu);

    double coeff = 0.2;

    lprintf("SANITY CHECK", 0, "suNg_av_field not zero: %0.15e\n", sqnorm_suNg_av_field_cpu(force_gpu));
    lprintf("SANITY CHECK", 0, "cpu not zero: %0.15e\n", sqnorm_suNg_av_field_cpu(force_cpu));

    // GPU operation
    force0_kernel_gpu(force_gpu, coeff);

    // CPU operation
    _MASTER_FOR(&glattice, ix) {
        suNg s1, s2;
        suNg_algebra_vector f;
        for (int mu = 0; mu < 4; ++mu) {
            staples(ix, mu, &s1);
            _suNg_times_suNg_dagger(s2, *_4FIELD_AT(u_gauge, ix, mu), s1);

            /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
            _fund_algebra_project(f, s2);
            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_cpu, ix, mu), coeff, f);
        }
    }

    // Check that both force and gfield are identical.
    copy_suNg_field_cpu(in_cpu, u_gauge);
    copy_from_gpu_suNg_field(u_gauge);
    copy_from_gpu_suNg_av_field(force_gpu);

    lprintf("SANITY CHECK", 0, "in field not zero: %0.15e\n", sqnorm_suNg_field_cpu(u_gauge));
    lprintf("SANITY CHECK", 0, "cpu not zero: %0.15e\n", sqnorm_suNg_field_cpu(in_cpu));

    lprintf("SANITY CHECK", 0, "suNg_av_field not zero: %0.15e\n", sqnorm_suNg_av_field_cpu(force_gpu));
    lprintf("SANITY CHECK", 0, "cpu not zero: %0.15e\n", sqnorm_suNg_av_field_cpu(force_cpu));

    sub_assign_suNg_field_cpu(in_cpu, u_gauge);
    sub_assign_suNg_av_field_cpu(force_cpu, force_gpu);

    lprintf("TEST", 0, "Checking Gauge Field\n");
    double diff_norm = sqrt(sqnorm_suNg_field_cpu(in_cpu));
    return_val += check_diff_norm(diff_norm, 1e-13);

    lprintf("TEST", 0, "Checking force field\n");
    diff_norm = sqrt(sqnorm_suNg_av_field_cpu(force_cpu));
    return_val += check_diff_norm_zero(diff_norm);

    finalize_process();

    return return_val;
}