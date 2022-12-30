/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* coherence of the dirac float with the dirac op
*
*******************************************************************************/

#include "libhr.h"

static double hmass = 0.1;

static void loc_D(spinor_field *out, spinor_field *in)
{
    Dphi(hmass, out, in);
}

static void loc_D_flt(spinor_field_flt *out, spinor_field_flt *in)
{
    Dphi_flt(hmass, out, in);
}

static void loc_D_cpu(spinor_field *out, spinor_field *in) 
{
    Dphi_cpu(hmass, out, in);
}

static void loc_D_flt_cpu(spinor_field_flt *out, spinor_field_flt *in) 
{
    Dphi_flt_cpu(hmass, out, in);
}

int test_coherence_gpu() 
{
    lprintf("INFO", 0, "Testing for GPU\n");
    int return_value = 0;

    spinor_field *s0, *s1;
    spinor_field_flt *f0, *f1;
    double sig, tau;

    s0 = alloc_spinor_field_f(1, &glattice);
    s1 = alloc_spinor_field_f(1, &glattice);
    f0 = alloc_spinor_field_f_flt(1, &glattice);
    f1 = alloc_spinor_field_f_flt(1, &glattice);

    random_u(u_gauge);
    represent_gauge_field();
    gaussian_spinor_field(s0);
    tau = 1. / sqrt(spinor_field_sqnorm_f_cpu(s0));
    spinor_field_mul_f_cpu(s0, tau, s0);
    assign_sd2s(f0, s0);
    assign_ud2u_f();

    copy_to_gpu_spinor_field_f(s0);
    copy_to_gpu_spinor_field_f(s1);
    copy_to_gpu_spinor_field_f_flt(f0);
    copy_to_gpu_spinor_field_f_flt(f1);
    copy_to_gpu_gfield_f(u_gauge_f);
    copy_to_gpu_gfield_f_flt(u_gauge_f_flt);

    #ifdef WITH_MPI
        start_sendrecv_gpu_gfield_f(u_gauge_f);
        complete_sendrecv_gpu_gfield_f(u_gauge_f);

        start_sendrecv_gpu_gfield_f_flt(u_gauge_f_flt);
        complete_sendrecv_gpu_gfield_f_flt(u_gauge_f_flt);
    #endif

    loc_D(s1, s0);
    loc_D_flt(f1, f0);

    #ifdef WITH_GPU
        copy_from_gpu_spinor_field_f(s0);
        copy_from_gpu_spinor_field_f(s1);
        copy_from_gpu_spinor_field_f_flt(f0);
        copy_from_gpu_spinor_field_f_flt(f1);
    #endif

    lprintf("INFO", 0, "Spinor sqnorm result double precision: %0.2e\n", spinor_field_sqnorm_f_cpu(s1));
    lprintf("INFO", 0, "Spinor sqnorm result single precision: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(f1));

    assign_sd2s(f0, s1);

    spinor_field_mul_add_assign_f_flt_cpu(f0, -1.0, f1);
    sig = spinor_field_sqnorm_f_flt_cpu(f0);

    lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", sqrt(sig));
    lprintf("MAIN", 0, "(should be around 1*10^(-8) or so)\n\n");

    if (sqrt(sig) > 10.e-7)
        return_value = 1;

    free_spinor_field_f(s0);
    free_spinor_field_f_flt(f0);
    free_spinor_field_f(s1);
    free_spinor_field_f_flt(f1);
    return return_value;
}

int test_coherence_cpu() 
{
    lprintf("INFO", 0, "Testing for CPU\n");
    int return_value = 0;

    spinor_field *s0, *s1;
    spinor_field_flt *f0, *f1;
    double sig, tau;

    s0 = alloc_spinor_field_f(1, &glattice);
    s1 = alloc_spinor_field_f(1, &glattice);
    f0 = alloc_spinor_field_f_flt(1, &glattice);
    f1 = alloc_spinor_field_f_flt(1, &glattice);

    random_u(u_gauge);
    represent_gauge_field();
    gaussian_spinor_field(s0);
    tau = 1. / sqrt(spinor_field_sqnorm_f_cpu(s0));
    spinor_field_mul_f_cpu(s0, tau, s0);
    assign_sd2s(f0, s0);
    assign_ud2u_f();

    loc_D_cpu(s1, s0);
    loc_D_flt_cpu(f1, f0);

    lprintf("INFO", 0, "Spinor sqnorm result double precision: %0.2e\n", spinor_field_sqnorm_f_cpu(s1));
    lprintf("INFO", 0, "Spinor sqnorm result single precision: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(f1));

    assign_sd2s(f0, s1);

    spinor_field_mul_add_assign_f_flt_cpu(f0, -1.0, f1);
    sig = spinor_field_sqnorm_f_flt_cpu(f0);

    lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", sqrt(sig));
    lprintf("MAIN", 0, "(should be around 1*10^(-8) or so)\n\n");

    if (sqrt(sig) > 10.e-7)
        return_value = 1;

    free_spinor_field_f(s0);
    free_spinor_field_f_flt(f0);
    free_spinor_field_f(s1);
    free_spinor_field_f_flt(f1);
    return return_value;
}

int main(int argc, char *argv[]) 
{
    int return_value = 0;

    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    setup_gauge_fields();
    u_gauge_f_flt = alloc_gfield_f_flt(&glattice); // TODO: This should be done in the geom setup (SAM)

    return_value += test_coherence_gpu();
    return_value += test_coherence_cpu();

    finalize_process();
    return return_value;
}
