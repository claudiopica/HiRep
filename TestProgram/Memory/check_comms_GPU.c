/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
* NOCOMPILE= !WITH_MPI
*
* Check that copy sync and buffer comms execute without errors
* and don't change the sqnorm of the field
*
*******************************************************************************/

// TODO: Odd sfield comms & single precision spinor not working

#include "libhr.h"
#include <string.h>

int test_comms_spinor_field_f(geometry_descriptor*);
int test_comms_spinor_field_f_flt(geometry_descriptor*);
int test_comms_sfield(geometry_descriptor*);
int test_comms_gfield();
int test_comms_gfield_flt();
int test_comms_gfield_f();
int test_comms_gfield_f_flt();
int test_comms_suNg_scalar_field();
int test_comms_avfield();
int test_comms_gtransf();
int test_comms_clover_ldl();
int test_comms_clover_term();
int test_comms_clover_force();

int test_sync_identical_to_cpu_spinor_field_f(geometry_descriptor*);
int test_sync_identical_to_cpu_spinor_field_f_flt(geometry_descriptor *gd);
int test_sync_identical_to_cpu_sfield(geometry_descriptor *gd);
int test_sync_identical_to_cpu_gfield();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Run tests
      /* Double precision */
    lprintf("INFO", 0, "\n\nFull lattice tests\n\n");
    return_val += test_comms_spinor_field_f(&glattice);
    return_val += test_comms_sfield(&glattice);
    return_val += test_comms_gfield();
    return_val += test_comms_gfield_f();
    return_val += test_comms_suNg_scalar_field();
    return_val += test_comms_avfield();
    return_val += test_comms_gtransf();
    return_val += test_comms_clover_ldl();
    return_val += test_comms_clover_term();
    return_val += test_comms_clover_force();

      /* Single precision */
    return_val += test_comms_spinor_field_f_flt(&glattice);
    return_val += test_comms_gfield_flt();
    return_val += test_comms_gfield_f_flt();

    lprintf("INFO", 0, "\n\nSpinor tests on even lattice\n\n");
    return_val += test_comms_spinor_field_f(&glat_even);
    return_val += test_comms_spinor_field_f_flt(&glat_even);
    return_val += test_comms_sfield(&glat_even);

    lprintf("INFO", 0, "\n\nSpinor tests on odd lattice\n\n");
    return_val += test_comms_spinor_field_f(&glat_odd);
    return_val += test_comms_spinor_field_f_flt(&glat_odd);
    return_val += test_comms_sfield(&glat_odd);

    /* Test sync */
    return_val += test_sync_identical_to_cpu_spinor_field_f(&glattice);
    return_val += test_sync_identical_to_cpu_spinor_field_f(&glat_even);
    return_val += test_sync_identical_to_cpu_spinor_field_f(&glat_odd);

    return_val += test_sync_identical_to_cpu_spinor_field_f_flt(&glattice);
    return_val += test_sync_identical_to_cpu_spinor_field_f_flt(&glat_even);
    return_val += test_sync_identical_to_cpu_spinor_field_f_flt(&glat_odd);

    return_val += test_sync_identical_to_cpu_sfield(&glattice);
    return_val += test_sync_identical_to_cpu_sfield(&glat_even);
    return_val += test_sync_identical_to_cpu_sfield(&glat_odd);

    return_val += test_sync_identical_to_cpu_gfield();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_comms_spinor_field_f(geometry_descriptor *gd) 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    lprintf("INFO", 0, "Sqnorm testing\n");

    // Setup fields on GPU
    int return_val = 0;
    spinor_field *f = alloc_spinor_field_f(1, gd);
    gaussian_spinor_field(f);
    copy_to_gpu_spinor_field_f(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = spinor_field_sqnorm_f(f);
    lprintf("SANITY CHECK", 0, "[In field GPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    
    // Execute communications
    start_sendrecv_gpu_spinor_field_f(f);
    complete_sendrecv_gpu_spinor_field_f(f);

    // Evaluate sqnorm after comms
    double sqnorm_end = spinor_field_sqnorm_f(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);// TODO: Check diff not diff norm (SAM)
    
    free_spinor_field_f(f);

    return return_val;
}

int test_sync_identical_to_cpu_spinor_field_f(geometry_descriptor *gd) 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    int return_val = 0;
    spinor_field *f = alloc_spinor_field_f(1, gd);
    gaussian_spinor_field(f);
    copy_to_gpu_spinor_field_f(f);

    //gd = &glattice;
    int sendbuf_len = 0;
    box_t *L = geometryBoxes->next;
    int m = 0;
    while (L && m < glattice.nbuffers_spinor) 
    {
        sendbuf_len += boxVolume(L->sendBox);
        L=L->next; m++;
    }

    // Sync to buffer on CPU and save the sendbuffer in an array
    sync_field(gd, sizeof(*(f->ptr)), 1, f->ptr, f->sendbuf_ptr);
    suNf_spinor* sendbuf_cpu = (suNf_spinor*)malloc(sendbuf_len*sizeof(suNf_spinor));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, sendbuf_len*sizeof(suNf_spinor));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_spinor_field_f(f);
    suNf_spinor* sendbuf_gpu = (suNf_spinor*)malloc(sendbuf_len*sizeof(suNf_spinor));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, sendbuf_len*sizeof(suNf_spinor), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNf_spinor *spinor_gpu = (suNf_spinor*)malloc(sizeof(suNf_spinor));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            suNf_spinor *spinor_cpu = sendbuf_cpu + j + gd->sbuf_start[i];
            suNf_spinor *in_block = sendbuf_gpu + gd->sbuf_start[i];
            int stride = gd->sbuf_len[i];
            read_gpu_suNf_spinor(stride, *(spinor_gpu), in_block, j, 0);

            for (int k = 0; k < 4; ++k) 
            {
                for (int comp = 0; comp < NF; ++comp) 
                {
                    diff += creal((*spinor_cpu).c[k].c[comp])-creal((*spinor_gpu).c[k].c[comp]);
                    diff += cimag((*spinor_cpu).c[k].c[comp])-cimag((*spinor_gpu).c[k].c[comp]);
                }
            }
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_spinor_field_f(f);
    return return_val;
}

int test_comms_spinor_field_f_flt(geometry_descriptor *gd) 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_f_flt(1, gd);
    gaussian_spinor_field_flt(f);
    copy_to_gpu_spinor_field_f_flt(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = spinor_field_sqnorm_f_flt(f);
    lprintf("SANITY CHECK", 0, "[In field GPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    
    // Execute communications
    start_sendrecv_gpu_spinor_field_f_flt(f);
    complete_sendrecv_gpu_spinor_field_f_flt(f);

    // Evaluate sqnorm after comms
    double sqnorm_end = spinor_field_sqnorm_f_flt(f);
    lprintf("SANITY CHECK", 0, "[Out field GPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);// TODO: Check diff not diff norm (SAM)

    free_spinor_field_f_flt(f);

    return return_val;
}

int test_sync_identical_to_cpu_spinor_field_f_flt(geometry_descriptor *gd) 
{
    lprintf("INFO", 0, " ======= TEST SINGLE PRECISION SPINOR FIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_f_flt(1, gd);
    gaussian_spinor_field_flt(f);
    copy_to_gpu_spinor_field_f_flt(f);

    //gd = &glattice;
    int sendbuf_len = 0;
    box_t *L = geometryBoxes->next;
    int m = 0;
    while (L && m < glattice.nbuffers_spinor) 
    {
        sendbuf_len += boxVolume(L->sendBox);
        L=L->next; m++;
    }

    // Sync to buffer on CPU and save the sendbuffer in an array
    sync_field(gd, sizeof(*(f->ptr)), 1, f->ptr, f->sendbuf_ptr);
    suNf_spinor_flt* sendbuf_cpu = (suNf_spinor_flt*)malloc(sendbuf_len*sizeof(suNf_spinor_flt));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, sendbuf_len*sizeof(suNf_spinor_flt));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_spinor_field_f_flt(f);
    suNf_spinor_flt* sendbuf_gpu = (suNf_spinor_flt*)malloc(sendbuf_len*sizeof(suNf_spinor_flt));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, sendbuf_len*sizeof(suNf_spinor_flt), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNf_spinor_flt *spinor_gpu = (suNf_spinor_flt*)malloc(sizeof(suNf_spinor_flt));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            suNf_spinor_flt *spinor_cpu = sendbuf_cpu + j + gd->sbuf_start[i];
            suNf_spinor_flt *in_block = sendbuf_gpu + gd->sbuf_start[i];
            int stride = gd->sbuf_len[i];
            read_gpu_suNf_spinor_flt(stride, *(spinor_gpu), in_block, j, 0);

            for (int k = 0; k < 4; ++k) 
            {
                for (int comp = 0; comp < NF; ++comp) 
                {
                    diff += creal((*spinor_cpu).c[k].c[comp])-creal((*spinor_gpu).c[k].c[comp]);
                    diff += cimag((*spinor_cpu).c[k].c[comp])-cimag((*spinor_gpu).c[k].c[comp]);
                }
            }
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_spinor_field_f_flt(f);
    return return_val;
}

int test_comms_sfield(geometry_descriptor *gd) 
{
    lprintf("INFO", 0, " ======= TEST SFIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    scalar_field *f = alloc_sfield(1, gd);
    random_sfield_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_sfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_sfield(f);
    
    // Execute communications
    start_sendrecv_gpu_sfield(f);
    complete_sendrecv_gpu_sfield(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_sfield(f);
    double sqnorm_end = sqnorm_sfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_sfield(f);

    return return_val;
}

int test_sync_identical_to_cpu_sfield(geometry_descriptor *gd) 
{
    lprintf("INFO", 0, " ======= TEST SFIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    int return_val = 0;
    scalar_field *f = alloc_sfield(1, gd);
    random_sfield_cpu(f);
    copy_to_gpu_sfield(f);

    //gd = &glattice;
    int sendbuf_len = 0;
    box_t *L = geometryBoxes->next;
    int m = 0;
    while (L && m < glattice.nbuffers_spinor) 
    {
        sendbuf_len += boxVolume(L->sendBox);
        L=L->next; m++;
    }

    // Sync to buffer on CPU and save the sendbuffer in an array
    sync_field(gd, sizeof(*(f->ptr)), 1, f->ptr, f->sendbuf_ptr);
    double* sendbuf_cpu = (double*)malloc(sendbuf_len*sizeof(double));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, sendbuf_len*sizeof(double));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_sfield(f);
    double* sendbuf_gpu = (double*)malloc(sendbuf_len*sizeof(double));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, sendbuf_len*sizeof(double), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    double *spinor_gpu = (double*)malloc(sizeof(double));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            double *spinor_cpu = sendbuf_cpu + j + gd->sbuf_start[i];
            double *in_block = sendbuf_gpu + gd->sbuf_start[i];
            read_gpu_double(gd->sbuf_len[i], *(spinor_gpu), in_block, j, 0);
            diff = (*spinor_cpu) - (*spinor_gpu);
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_sfield(f);
    return return_val;
}

int test_comms_gfield() 
{
    lprintf("INFO", 0, " ======= TEST GFIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_field *f = alloc_gfield(&glattice);
    random_gfield_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gfield(f);
    
    // Execute communications
    start_sendrecv_gpu_gfield(f);
    complete_sendrecv_gpu_gfield(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gfield(f);
    double sqnorm_end = sqnorm_gfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gfield(f);

    return return_val;
}

int test_sync_identical_to_cpu_gfield() 
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNg_field *f = alloc_gfield(gd);
    random_gfield_cpu(f);
    copy_to_gpu_gfield(f);

    //gd = &glattice;
    int sendbuf_len = 0;
    box_t *L = geometryBoxes->next;
    int m = 0;
    while (L && m < glattice.nbuffers_spinor) 
    {
        sendbuf_len += boxVolume(L->sendBox);
        L=L->next; m++;
    }

    // Sync to buffer on CPU and save the sendbuffer in an array
    sync_field(gd, sizeof(*(f->ptr)), 1, f->ptr, f->sendbuf_ptr);
    suNg* sendbuf_cpu = (suNg*)malloc(sendbuf_len*sizeof(suNg));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, sendbuf_len*sizeof(suNg));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_gfield(f);
    suNg* sendbuf_gpu = (suNg*)malloc(sendbuf_len*sizeof(suNg));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, sendbuf_len*sizeof(suNg), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNg *mat_gpu = (suNg*)malloc(sizeof(suNg));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            suNg *in_block = sendbuf_gpu + 4*gd->sbuf_start[i];
            int stride = gd->sbuf_len[i];

            for (int comp = 0; comp < NG; ++comp) 
            {
                int base_index = j + gd->sbuf_start[i];
                suNg *mat_cpu = sendbuf_cpu + coord_to_index(base_index, comp);
                read_gpu_suNg(stride, *(mat_gpu), in_block, j, comp);
                printf("ix: %d, cpu: %0.5e + i%0.5e, gpu: %0.5e + i%0.5e\n", j + gd->sbuf_start[i], 
                    creal((*mat_cpu).c[comp]), cimag((*mat_cpu).c[comp]), 
                    creal((*mat_gpu).c[comp]), cimag((*mat_gpu).c[comp]));
                diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
                printf("diff: %0.5e\n", diff);
            }
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_gfield(f);
    return return_val;
}

int test_comms_gfield_flt() 
{
    lprintf("INFO", 0, " ======= TEST GFIELD SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_field_flt *f = alloc_gfield_flt(&glattice);
    random_gfield_flt_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gfield_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gfield_flt(f);
    
    // Execute communications
    start_sendrecv_gpu_gfield_flt(f);
    complete_sendrecv_gpu_gfield_flt(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gfield_flt(f);
    double sqnorm_end = sqnorm_gfield_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gfield_flt(f);

    return return_val;
}

int test_comms_gfield_f() 
{
    lprintf("INFO", 0, " ======= TEST GFIELD FERMION REP ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNf_field *f = alloc_gfield_f(&glattice);
    random_gfield_f_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gfield_f_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gfield_f(f);
    
    // Execute communications
    start_sendrecv_gpu_gfield_f(f);
    complete_sendrecv_gpu_gfield_f(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gfield_f(f);
    double sqnorm_end = sqnorm_gfield_f_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gfield_f(f);

    return return_val;
}

int test_comms_gfield_f_flt() 
{
    lprintf("INFO", 0, " ======= TEST GFIELD FERMION REP SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNf_field_flt *f = alloc_gfield_f_flt(&glattice);
    random_gfield_f_flt_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gfield_f_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gfield_f_flt(f);
    
    // Execute communications
    start_sendrecv_gpu_gfield_f_flt(f);
    complete_sendrecv_gpu_gfield_f_flt(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gfield_f_flt(f);
    double sqnorm_end = sqnorm_gfield_f_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gfield_f_flt(f);

    return return_val;
}

int test_comms_suNg_scalar_field() 
{
    lprintf("INFO", 0, " ======= TEST SU(NG) SCALAR FIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_scalar_field *f = alloc_suNg_scalar_field(&glattice);
    random_suNg_scalar_field_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_suNg_scalar_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_suNg_scalar_field(f);
    
    // Execute communications
    start_sendrecv_gpu_suNg_scalar_field(f);
    complete_sendrecv_gpu_suNg_scalar_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_suNg_scalar_field(f);
    double sqnorm_end = sqnorm_suNg_scalar_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_suNg_scalar_field(f);

    return return_val;
}

int test_comms_avfield() 
{
    lprintf("INFO", 0, " ======= TEST ALGEBRA VECTOR FIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_av_field *f = alloc_avfield(&glattice);
    random_avfield_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_avfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_avfield(f);
    
    // Execute communications
    start_sendrecv_gpu_avfield(f);
    complete_sendrecv_gpu_avfield(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_avfield(f);
    double sqnorm_end = sqnorm_avfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_avfield(f);

    return return_val;
}

int test_comms_gtransf() 
{
    lprintf("INFO", 0, " ======= TEST GAUGE TRANSFORMATION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_field *f = alloc_gtransf(&glattice);
    random_gtransf_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gtransf_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gtransf(f);
    
    // Execute communications
    start_sendrecv_gpu_gtransf(f);
    complete_sendrecv_gpu_gtransf(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gtransf(f);
    double sqnorm_end = sqnorm_gtransf_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gtransf(f);

    return return_val;
}

int test_comms_clover_ldl() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER LDL ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    ldl_field *f = alloc_clover_ldl(&glattice);
    random_clover_ldl_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_clover_ldl_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_clover_ldl(f);
    
    // Execute communications
    start_sendrecv_gpu_clover_ldl(f);
    complete_sendrecv_gpu_clover_ldl(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_clover_ldl(f);
    double sqnorm_end = sqnorm_clover_ldl_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_clover_ldl(f);

    return return_val;
}

int test_comms_clover_term() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER TERM ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNfc_field *f = alloc_clover_term(&glattice);
    random_clover_term_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_clover_term_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_clover_term(f);
    
    // Execute communications
    start_sendrecv_gpu_clover_term(f);
    complete_sendrecv_gpu_clover_term(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_clover_term(f);
    double sqnorm_end = sqnorm_clover_term_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_clover_term(f);

    return return_val;
}

int test_comms_clover_force() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER FORCE ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNf_field *f = alloc_clover_force(&glattice);
    random_clover_force_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_clover_force_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_clover_force(f);
    
    // Execute communications
    start_sendrecv_gpu_clover_force(f);
    complete_sendrecv_gpu_clover_force(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_clover_force(f);
    double sqnorm_end = sqnorm_clover_force_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_clover_force(f);

    return return_val;
}
