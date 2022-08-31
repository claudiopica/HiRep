#define MAIN_PROGRAM

#include <stdbool.h>
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "global.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "logger.h"
#include "setup.h"
#include "hr_complex.h"
#include "test_hermitian.h"
#include "helper_functions.h"

// TODO: Add explanation/comments

int main(int argc, char *argv[])
{
    // Init
    int pass;
    init_test(argc, argv);

    // Test Block
    //run_test(test_hermiticity(&I_operator, &I_operator_cpu, "Unit operator"), pass);
    //cudaDeviceReset();

    run_test(test_hermiticity(&Q_operator, &Q_operator_cpu, "Q = g5Dphi"), pass);

    finalize_process();
    return pass;
}

bool test_hermiticity(spinor_operator S, spinor_operator S_cpu, char *name)
{
    lprintf("INFO", 0, "[Testing %s]\n", name);

    spinor_field *s1, *s2, *S_s1, *S_s2, *S_s1_cpu, *S_s2_cpu;
    int return_val = 0;
    
    hr_complex tau, tau_cpu;
    s1 = setup_infield();
    s2 = setup_infield();
    S_s1 = alloc_spinor_field_f(1, &glattice);
    S_s2 = alloc_spinor_field_f(1, &glattice);
    //S_s1_cpu = alloc_spinor_field_f(1, &glattice);
    //S_s2_cpu = alloc_spinor_field_f(1, &glattice);
    lprintf("INFO", 0, "[Infield norm before: %0.15lf]\n", spinor_field_sqnorm_f(s1));
    lprintf("INFO", 0, "[Infield2 norm before: %0.15lf]\n", spinor_field_sqnorm_f(s2));
    lprintf("INFO", 0, "[Outfield norm before: %0.15lf]\n", spinor_field_sqnorm_f(S_s1));
    lprintf("INFO", 0, "[Outfield2 norm before: %0.15lf]\n", spinor_field_sqnorm_f(S_s2));
    bool pass_gpu = is_hermitian_on_GPU(s1, s2, S_s1, S_s2, S);
    //bool pass_cpu = is_hermitian_on_CPU(s1, s2, S_s1_cpu, S_s2_cpu, S_cpu);
    lprintf("INFO", 0, "[Infield norm: %0.15lf]\n", spinor_field_sqnorm_f(s1));
    lprintf("INFO", 0, "[Infield2 norm: %0.15lf]\n", spinor_field_sqnorm_f(s2));
    lprintf("INFO", 0, "[Outfield norm: %0.15lf]\n", spinor_field_sqnorm_f(S_s1));
    lprintf("INFO", 0, "[Outfield2 norm: %0.15lf]\n", spinor_field_sqnorm_f(S_s2));
    
    bool pass_sanity_check = result_spinor_fields_not_identically_zero_gpu(S_s1, S_s2);
    
    //bool pass_sanity_check_cpu = result_spinor_fields_not_identically_zero_cpu(S_s1_cpu, S_s2_cpu);
    //bool are_copies_identical = gpu_and_cpu_copies_identical(S_s1, S_s2);

    /*spinor_field_copy_from_gpu_f(S_s1);
    spinor_field_copy_from_gpu_f(S_s2);

    spinor_field *diff_1, *diff_2;
    diff_1 = alloc_spinor_field_f(1, &glattice);
    diff_2 = alloc_spinor_field_f(1, &glattice);
    spinor_field_sub_f_cpu(diff_1, S_s1, S_s1_cpu);
    spinor_field_sub_f_cpu(diff_2, S_s2, S_s2_cpu);

    lprintf("INFO", 0, "[Diff norm gpu-cpu: %0.15lf]\n", spinor_field_sqnorm_f_cpu(diff_1));
    lprintf("INFO", 0, "[Diff norm gpu-cpu: %0.15lf]\n", spinor_field_sqnorm_f_cpu(diff_2));*/

    free_spinors(&s1, &s2, &S_s1, &S_s2);
    return pass_gpu && pass_sanity_check;
}

bool is_hermitian_on_GPU(spinor_field *s1, spinor_field *s2, 
			       spinor_field *S_s1, spinor_field *S_s2, 
			       spinor_operator S) 
{
    hr_complex tau, N;
    S(S_s1, s1);
    S(S_s2, s2);
    N = sqrt(spinor_field_sqnorm_f(s1)*spinor_field_sqnorm_f(s2));
    tau = (spinor_field_prod_f(S_s2, s1) - spinor_field_prod_f(s2, S_s1))/N;

    bool pass = fabs(_complex_re(tau)) < 1.e-14 && fabs(_complex_im(tau)) < 1.e-14;
    if (!pass) 
    {
        lprintf("FAILED", 0, "The operator is not hermitian on the GPU.\n");
    }
    lprintf("RESULT", 0, "[diff gpu = %0.20lf + i%0.20lf]\n", _complex_re(tau), _complex_im(tau));
    return pass;
}

bool is_hermitian_on_CPU(spinor_field *s1, spinor_field *s2, 
		 	       spinor_field *S_s1, spinor_field *S_s2, 
			       spinor_operator S) 
{
    hr_complex tau_cpu, N;
    S(S_s1, s1);
    S(S_s2, s2);
    N = sqrt(spinor_field_sqnorm_f_cpu(s1)*spinor_field_sqnorm_f_cpu(s2));
    tau_cpu = (spinor_field_prod_f_cpu(S_s2, s1) - spinor_field_prod_f_cpu(s2, S_s1))/N;

    bool pass = fabs(_complex_re(tau_cpu)) < 1.e-14 && fabs(_complex_im(tau_cpu)) < 1.e-14;
    if (!pass) 
    {
        lprintf("FAILED", 0, "The operator is not hermitian on the CPU.\n");
    }
    lprintf("RESULT", 0, "[diff cpu = %0.20lf + i%0.20lf]\n", 
		                         _complex_re(tau_cpu), _complex_im(tau_cpu));
    return pass;
}

bool result_spinor_fields_not_identically_zero_gpu(spinor_field *S_s1, spinor_field *S_s2) 
{
    bool pass_gpu = spinor_field_sqnorm_f(S_s1)!=0 && spinor_field_sqnorm_f(S_s2)!=0;
    if (!pass_gpu) 
    {
        lprintf("FAILED", 0, "Result spinor fields are identically zero on GPU.\n");
    }
    return pass_gpu;
}

bool result_spinor_fields_not_identically_zero_cpu(spinor_field *S_s1_cpu, spinor_field *S_s2_cpu) 
{
    bool pass_cpu = spinor_field_sqnorm_f_cpu(S_s1_cpu)!=0 && spinor_field_sqnorm_f_cpu(S_s2_cpu)!=0;
    if (!pass_cpu) 
    {
        lprintf("FAILED", 0, "Result spinor fields are identically zero on CPU.\n");
    }
    return pass_cpu;
}

bool gpu_and_cpu_copies_identical(spinor_field *S_s1, spinor_field *S_s2) 
{
    double diff_1, diff_2;
    
    diff_1 = spinor_field_sqnorm_f(S_s1) - spinor_field_sqnorm_f_cpu(S_s1);
    diff_2 = spinor_field_sqnorm_f(S_s2) - spinor_field_sqnorm_f_cpu(S_s2);

    bool pass = fabs(diff_1) < 1.e-11 && fabs(diff_2) < 1.e-11;
    if (!pass) {
        lprintf("FAILED", 0, "Operations on CPU and GPU do not yield the same result.\n");
    }
    lprintf("RESULT", 0, "[diff gpu-cpu s1: %0.20lf]\n", diff_1);
    lprintf("RESULT", 0, "[diff gpu-cpu s2: %0.20lf]\n", diff_2);
    return pass;
}

/*
 * Frees used spinor fields
 *
 * in: (**)spinor_field       All the fields to be freed.
 */
void free_spinors(spinor_field **s1,
                  spinor_field **s2,
                  spinor_field **S_s1,
                  spinor_field **S_s2)
{
    free_spinor_field_f(*s1);
    free_spinor_field_f(*s2);
    free_spinor_field_f(*S_s1);
    free_spinor_field_f(*S_s2);
}

/* ============== OPERATOR DEFINITIONS ==========================*/

static double hmass = 0.1;

void Q_operator(spinor_field *out, spinor_field *in)
{
    /*
     * Implement this for the actual GPU operator.
     * Issues here are
     *
     * 1. The operator does not work yet, because of an old implementation that does not
     *    work for newer versions of CUDA yet
     * 2. There is some spinor allocation issue inside this operator. Using g5Dphi on the
     *    GPU will also affect the CPU copy (for some reason) and as a result both test segments
     *    will fail. This is, however, not an issue of the CPU code.
     * 3. There is a temporary replacement of g5Dphi of spinor multiplication which will pass.
     *
     * */
    //spinor_field_mul_f(out, 2.5, in);
    g5Dphi(-hmass, out, in);
}

void Q_operator_cpu(spinor_field *out, spinor_field *in) 
{
    g5Dphi_cpu(-hmass, out, in);
}

void I_operator(spinor_field *out, spinor_field *in) 
{
    spinor_field_mul_f(out, 1, in);
}

void I_operator_cpu(spinor_field *out, spinor_field *in) 
{
    spinor_field_mul_f_cpu(out, 1, in);
}
