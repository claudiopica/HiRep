/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that CPU and GPU copies are identical
*
********************************************************************************/

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
#include "random.h"
#include "representation.h"

int test_hermiticity(spinor_operator, spinor_operator, char*);
void Q_operator(spinor_field*, spinor_field*);
void Q_operator_cpu(spinor_field*, spinor_field*);
void I_operator(spinor_field*, spinor_field*);
void I_operator_cpu(spinor_field*, spinor_field*);

int main(int argc, char *argv[])
{
    // Init
    int return_val = 0;

    // Setup process and communication
    setup_process(&argc, &argv);

    // Setup gauge field
    setup_gauge_fields();
    random_u(u_gauge);
    represent_gauge_field();
    copy_to_gpu_gfield_f(u_gauge_f);

    // Test Block
    return_val += test_hermiticity(&I_operator, &I_operator_cpu, "Unit operator");
    return_val += test_hermiticity(&Q_operator, &Q_operator_cpu, "Q = g5Dphi");

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_hermiticity(spinor_operator S, spinor_operator S_cpu, char *name)
{
    lprintf("INFO", 0, "[Testing %s]\n", name);

    spinor_field *s, *S_s, *S_s_cpu, *diff;
    int return_val = 0;
    
    hr_complex tau, tau_cpu;

    s = alloc_spinor_field_f(1, &glattice);
    S_s = alloc_spinor_field_f(1, &glattice);
    S_s_cpu = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(s);
    copy_to_gpu_spinor_field_f(s);

    lprintf("INFO", 0, "Input spinor field norm CPU: %0.2e\n", spinor_field_sqnorm_f_cpu(s));
    lprintf("INFO", 0, "Input spinor field norm GPU: %0.2e\n", spinor_field_sqnorm_f(s));

    S(S_s, s);
    S_cpu(S_s_cpu, s);

    copy_from_gpu_spinor_field_f(S_s);
    
    // Sanity checks: Norms are not identically zero
    lprintf("INFO", 0, "Output spinor field norm GPU: %0.2e\n", spinor_field_sqnorm_f_cpu(S_s));
    lprintf("INFO", 0, "Output spinor field norm CPU: %0.2e\n", spinor_field_sqnorm_f_cpu(S_s_cpu));

    spinor_field_sub_assign_f_cpu(S_s_cpu, S_s);

    double diff_norm = spinor_field_sqnorm_f_cpu(S_s_cpu);
    if (fabs(diff_norm) > 1e-14 || !isfinite(diff_norm)) 
    {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    
    lprintf("RESULT", 0, "[Diff norm gpu-cpu %0.2e]\n", diff_norm);

    free_spinor_field_f(s);
    free_spinor_field_f(S_s);
    free_spinor_field_f(S_s_cpu);
    return return_val;
}

/* ============== OPERATOR DEFINITIONS ==========================*/

static double hmass = 0.1;

void Q_operator(spinor_field *out, spinor_field *in)
{
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
