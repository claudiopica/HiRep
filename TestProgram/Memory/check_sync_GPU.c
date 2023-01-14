/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
* NOCOMPILE= !WITH_MPI
*
* Check that syncs yields the same sendbuffers up to geometry conversion
*
*******************************************************************************/

// TODO: Test does not work for gfield

#include "libhr.h"
#include <string.h>

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
    /*return_val += test_sync_identical_to_cpu_spinor_field_f(&glattice);
    return_val += test_sync_identical_to_cpu_spinor_field_f(&glat_even);
    return_val += test_sync_identical_to_cpu_spinor_field_f(&glat_odd);

    return_val += test_sync_identical_to_cpu_spinor_field_f_flt(&glattice);
    return_val += test_sync_identical_to_cpu_spinor_field_f_flt(&glat_even);
    return_val += test_sync_identical_to_cpu_spinor_field_f_flt(&glat_odd);

    return_val += test_sync_identical_to_cpu_sfield(&glattice);
    return_val += test_sync_identical_to_cpu_sfield(&glat_even);
    return_val += test_sync_identical_to_cpu_sfield(&glat_odd);*/

    return_val += test_sync_identical_to_cpu_gfield();

    // Finalize and return
    finalize_process();
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
    while (L && m < glattice.nbuffers_gauge) 
    {
        sendbuf_len += boxVolume(L->sendBox);
        L=L->next; m++;
    }

    // Sync to buffer on CPU and save the sendbuffer in an array
    sync_field(gd, sizeof(*(f->ptr)), 1, f->ptr, f->sendbuf_ptr);
    suNg* sendbuf_cpu = (suNg*)malloc(4*sendbuf_len*sizeof(suNg));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, 4*sendbuf_len*sizeof(suNg));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_gfield(f);
    suNg* sendbuf_gpu = (suNg*)malloc(4*sendbuf_len*sizeof(suNg));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, 4*sendbuf_len*sizeof(suNg), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNg *mat_gpu = (suNg*)malloc(sizeof(suNg));
    suNg *mat_cpu = (suNg*)malloc(sizeof(suNg));
    suNg *in_block;
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            int stride = gd->sbuf_len[i];
            int index = gd->sbuf_start[i];
            mat_cpu = _FIELD_AT_PTR(sendbuf_cpu, index, 0);
            in_block = _FIELD_AT_PTR(sendbuf_gpu, gd->sbuf_start[i], 0);
            for (int comp = 0; comp < NG*NG; ++comp) {
                read_gpu_suNg(stride, (*mat_gpu), in_block, j, comp);
                if (comp==0 && j==0) printf("ix: %d, comp: %d, cpu: %0.5e + i%0.5e, gpu: %0.5e + i%0.5e\n", 
                    j + gd->sbuf_start[i], comp, 
                    creal((*mat_cpu).c[comp]), cimag((*mat_cpu).c[comp]), 
                    creal((*mat_gpu).c[comp]), cimag((*mat_gpu).c[comp]));
                diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
            }            
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_gfield(f);
    return return_val;
}