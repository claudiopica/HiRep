/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that CPU and GPU copies are identical
*
********************************************************************************/

#include "libhr.h"

// TODO: test more operators

int test_identical(spinor_operator, spinor_operator, char*);
int test_identical_flt(spinor_operator_flt, spinor_operator_flt, char*);
int test_identical_massless(spinor_operator, spinor_operator, char*, geometry_descriptor*, geometry_descriptor*);
void Q_operator(spinor_field*, spinor_field*);
void Q_operator_cpu(spinor_field*, spinor_field*);
void Q_operator_flt(spinor_field_flt*, spinor_field_flt*);
void Q_operator_flt_cpu(spinor_field_flt*, spinor_field_flt*);
void D_operator(spinor_field*, spinor_field*);
void D_operator_cpu(spinor_field*, spinor_field*);
void D_operator_flt(spinor_field_flt*, spinor_field_flt*);
void D_operator_flt_cpu(spinor_field_flt*, spinor_field_flt*);
void I_operator(spinor_field*, spinor_field*);
void I_operator_cpu(spinor_field*, spinor_field*);
void D_massless(spinor_field*, spinor_field*);
void D_massless_cpu(spinor_field*, spinor_field*);

int main(int argc, char *argv[])
{
    // Init
    int return_val = 0;

    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary

    // Setup process and communication
    setup_process(&argc, &argv);

    // Setup gauge field
    setup_gauge_fields();

    random_u(u_gauge);
    represent_gauge_field();

    sync_single_precision_gauge_field();

    copy_from_gpu_gfield_f(u_gauge_f);
    copy_from_gpu_gfield_f_flt(u_gauge_f_flt);

    // Test Block
    return_val += test_identical(&I_operator, &I_operator_cpu, "Unit operator Full Lattice");
    return_val += test_identical(&Q_operator, &Q_operator_cpu, "Q = g5Dphi Full Lattice");
    return_val += test_identical_flt(&Q_operator_flt, &Q_operator_flt_cpu, "Q = g5Dphi_flt Full Lattice");
    return_val += test_identical(&D_operator, &D_operator_cpu, "D = Dphi Full Lattice");
    return_val += test_identical_flt(&D_operator_flt, &D_operator_flt_cpu, "D = Dphi_flt Full Lattice");

    return_val += test_identical_massless(&D_massless, &D_massless_cpu, "D = Dphi Massless Even Lattice", &glat_even, &glat_odd);
    return_val += test_identical_massless(&D_massless, &D_massless_cpu, "D = Dphi Massless Odd Lattice", &glat_odd, &glat_even);
    
    // Finalize and return
    finalize_process();
    return return_val;
}

int test_identical(spinor_operator S, spinor_operator S_cpu, char *name)
{
    lprintf("INFO", 0, "[Testing %s]\n", name);

    spinor_field *s, *S_s, *S_s_cpu;
    int return_val = 0;
    
    s = alloc_spinor_field_f(1, &glattice);
    S_s = alloc_spinor_field_f(1, &glattice);
    S_s_cpu = alloc_spinor_field_f(1, &glattice);

    gaussian_spinor_field(s);
    copy_to_gpu_spinor_field_f(s);

    lprintf("INFO", 0, "Input spinor field norm CPU: %0.2e\n", spinor_field_sqnorm_f_cpu(s));
    lprintf("INFO", 0, "Input spinor field norm GPU: %0.2e\n", spinor_field_sqnorm_f(s));

    #ifdef WITH_MPI
        start_sendrecv_gfield_f(u_gauge_f);
        complete_sendrecv_gfield_f(u_gauge_f);
    #endif

    S(S_s, s);
    S_cpu(S_s_cpu, s);

    //copy_from_gpu_spinor_field_f(S_s);
    copy_to_gpu_spinor_field_f(S_s_cpu);
    
    // Sanity checks: Norms are not identically zero
    lprintf("INFO", 0, "Output spinor field norm GPU: %0.2e\n", spinor_field_sqnorm_f(S_s));
    lprintf("INFO", 0, "Output spinor field norm CPU: %0.2e\n", spinor_field_sqnorm_f(S_s_cpu));

    spinor_field_sub_assign_f(S_s_cpu, S_s);

    double diff_norm = sqrt(spinor_field_sqnorm_f(S_s_cpu));
    return_val += check_finiteness(diff_norm);
    return_val += check_diff_norm(diff_norm, 1e-13);

    free_spinor_field_f(s);
    free_spinor_field_f(S_s);
    free_spinor_field_f(S_s_cpu);
    return return_val;
}

int test_identical_massless(spinor_operator S, spinor_operator S_cpu, char *name, geometry_descriptor* gd1, geometry_descriptor* gd2) 
{
    lprintf("INFO", 0, "[Testing %s]\n", name);

    spinor_field *s, *S_s, *S_s_cpu;
    int return_val = 0;

    s = alloc_spinor_field_f(1, gd1);
    S_s = alloc_spinor_field_f(1, gd2);
    S_s_cpu = alloc_spinor_field_f(1, gd2);

    gaussian_spinor_field(s);
    spinor_field_zero_f(S_s);
    copy_to_gpu_spinor_field_f(s);

    lprintf("INFO", 0, "Input spinor field norm CPU: %0.2e\n", spinor_field_sqnorm_f_cpu(s));
    lprintf("INFO", 0, "Input spinor field norm GPU: %0.2e\n", spinor_field_sqnorm_f(s));

    #ifdef WITH_MPI
        start_sendrecv_gfield_f(u_gauge_f);
        complete_sendrecv_gfield_f(u_gauge_f);
    #endif

    S(S_s, s);
    S_cpu(S_s_cpu, s);

    copy_from_gpu_spinor_field_f(S_s);
    
    // Sanity checks: Norms are not identically zero
    lprintf("INFO", 0, "Output spinor field norm GPU: %0.2e\n", spinor_field_sqnorm_f_cpu(S_s));
    lprintf("INFO", 0, "Output spinor field norm CPU: %0.2e\n", spinor_field_sqnorm_f_cpu(S_s_cpu));

    spinor_field_sub_assign_f_cpu(S_s_cpu, S_s);

    double diff_norm = sqrt(spinor_field_sqnorm_f_cpu(S_s_cpu));
    return_val += check_finiteness(diff_norm);
    return_val += check_diff_norm(diff_norm, 1e-13);

    free_spinor_field_f(s);
    free_spinor_field_f(S_s);
    free_spinor_field_f(S_s_cpu);
    return return_val;
}

int test_identical_flt(spinor_operator_flt S, spinor_operator_flt S_cpu, char *name)
{
    lprintf("INFO", 0, "[Testing %s]\n", name);

    spinor_field_flt *s, *S_s, *S_s_cpu;
    int return_val = 0;
    
    s = alloc_spinor_field_f_flt(1, &glattice);
    S_s = alloc_spinor_field_f_flt(1, &glattice);
    S_s_cpu = alloc_spinor_field_f_flt(1, &glattice);

    gaussian_spinor_field_flt(s);
    copy_to_gpu_spinor_field_f_flt(s);

    lprintf("INFO", 0, "Input spinor field norm CPU: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(s));
    lprintf("INFO", 0, "Input spinor field norm GPU: %0.2e\n", spinor_field_sqnorm_f_flt(s));

    #ifdef WITH_MPI
        start_sendrecv_gfield_f_flt(u_gauge_f_flt);
        complete_sendrecv_gfield_f_flt(u_gauge_f_flt);
    #endif

    S(S_s, s);
    S_cpu(S_s_cpu, s);

    copy_from_gpu_spinor_field_f_flt(S_s);
    
    // Sanity checks: Norms are not identically zero
    lprintf("INFO", 0, "Output spinor field norm GPU: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(S_s));
    lprintf("INFO", 0, "Output spinor field norm CPU: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(S_s_cpu));

    spinor_field_sub_assign_f_flt_cpu(S_s_cpu, S_s);

    double diff_norm = sqrt(spinor_field_sqnorm_f_flt_cpu(S_s_cpu));
    return_val += check_finiteness(diff_norm);
    return_val += check_diff_norm(diff_norm, 1e-4);

    free_spinor_field_f_flt(s);
    free_spinor_field_f_flt(S_s);
    free_spinor_field_f_flt(S_s_cpu);
    return return_val;
}

/* ============== OPERATOR DEFINITIONS ==========================*/

static double hmass = 0.1;

void D_massless(spinor_field *out, spinor_field *in) 
{
    Dphi_(out, in);
}

void D_massless_cpu(spinor_field *out, spinor_field *in) 
{
    Dphi_cpu_(out, in);
}

void D_operator(spinor_field *out, spinor_field *in) 
{
    Dphi(-hmass, out, in);
}

void D_operator_cpu(spinor_field *out, spinor_field *in) 
{
    Dphi_cpu(-hmass, out, in);
}

void D_operator_flt(spinor_field_flt *out, spinor_field_flt *in) 
{
    Dphi_flt(-hmass, out, in);
}

void D_operator_flt_cpu(spinor_field_flt *out, spinor_field_flt *in) 
{
    Dphi_flt_cpu(-hmass, out, in);
}

void Q_operator(spinor_field *out, spinor_field *in)
{
    g5Dphi(-hmass, out, in);
}

void Q_operator_cpu(spinor_field *out, spinor_field *in) 
{
    g5Dphi_cpu(-hmass, out, in);
}

void Q_operator_flt(spinor_field_flt *out, spinor_field_flt* in) 
{
    g5Dphi_flt(-hmass, out, in);
}

void Q_operator_flt_cpu(spinor_field_flt *out, spinor_field_flt* in) 
{
    g5Dphi_flt_cpu(-hmass, out, in);
}

void I_operator(spinor_field *out, spinor_field *in) 
{
    spinor_field_mul_f(out, 1, in);
}

void I_operator_cpu(spinor_field *out, spinor_field *in) 
{
    spinor_field_mul_f_cpu(out, 1, in);
}