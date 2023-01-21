/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
* NOCOMPILE= !WITH_MPI
*
* Check that syncs yields the same sendbuffers up to geometry conversion
*
*******************************************************************************/

// TODO: Differences over elementary site types

#include "libhr.h"
#include <string.h>

int test_sync_identical_to_cpu_spinor_field_f(geometry_descriptor*);
int test_sync_identical_to_cpu_spinor_field_f_flt(geometry_descriptor *gd);
int test_sync_identical_to_cpu_sfield(geometry_descriptor *gd);
int test_sync_identical_to_cpu_gfield();
int test_sync_identical_to_cpu_gfield_f();
int test_sync_identical_to_cpu_gfield_flt();
int test_sync_identical_to_cpu_gfield_f_flt();
int test_sync_identical_to_cpu_suNg_scalar_field();
int test_sync_identical_to_cpu_avfield();
int test_sync_identical_to_cpu_gtransf();
int test_sync_identical_to_cpu_clover_ldl();
int test_sync_identical_to_cpu_clover_term();
int test_sync_identical_to_cpu_clover_force();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Run tests
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
    return_val += test_sync_identical_to_cpu_gfield_f();
    return_val += test_sync_identical_to_cpu_gfield_flt();
    return_val += test_sync_identical_to_cpu_gfield_f_flt();
    return_val += test_sync_identical_to_cpu_suNg_scalar_field();
    return_val += test_sync_identical_to_cpu_avfield();
    return_val += test_sync_identical_to_cpu_gtransf();
    return_val += test_sync_identical_to_cpu_clover_ldl();
    return_val += test_sync_identical_to_cpu_clover_term();
    return_val += test_sync_identical_to_cpu_clover_force();

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
            read_gpu_suNf_spinor(stride, *(spinor_gpu), in_block, j, 0, 1);

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
            read_gpu_suNf_spinor_flt(stride, *(spinor_gpu), in_block, j, 0, 1);

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
            read_gpu_double(gd->sbuf_len[i], *(spinor_gpu), in_block, j, 0, 1);
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
    sync_field(f->type, 4*sizeof(*f->ptr), 0, f->ptr, f->sendbuf_ptr);
    suNg* sendbuf_cpu = (suNg*)malloc(4*sendbuf_len*sizeof(suNg));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, 4*sendbuf_len*sizeof(suNg));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_gfield(f);
    suNg* sendbuf_gpu = (suNg*)malloc(4*sendbuf_len*sizeof(suNg));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, 4*sendbuf_len*sizeof(suNg), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNg *mat_gpu = (suNg*)malloc(sizeof(suNg));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            for (int mu = 0; mu < 4; ++mu) {
                int idx = gd->sbuf_start[i] + j;
                int offset = gd->sbuf_start[i];
                suNg *mat_cpu = _4FIELD_AT_PTR(sendbuf_cpu, idx, mu, 0);
                suNg *in_block = _4FIELD_AT_PTR(sendbuf_gpu, offset, 0, 0);
                int stride = gd->sbuf_len[i];
                read_gpu_suNg(stride, (*mat_gpu), in_block, j, mu, 4);

                for (int comp = 0; comp < NG*NG; ++comp) {
                    diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                    diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
                }
            }     
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_gfield(f);
    return return_val;
}

int test_sync_identical_to_cpu_gfield_f() 
{
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNf_field *f = alloc_gfield_f(gd);
    random_gfield_f_cpu(f);
    copy_to_gpu_gfield_f(f);

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
    sync_field(f->type, 4*sizeof(*f->ptr), 0, f->ptr, f->sendbuf_ptr);
    suNf* sendbuf_cpu = (suNf*)malloc(4*sendbuf_len*sizeof(suNf));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, 4*sendbuf_len*sizeof(suNf));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_gfield_f(f);
    suNf* sendbuf_gpu = (suNf*)malloc(4*sendbuf_len*sizeof(suNf));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, 4*sendbuf_len*sizeof(suNf), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNf *mat_gpu = (suNf*)malloc(sizeof(suNf));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            for (int mu = 0; mu < 4; ++mu) {
                int idx = gd->sbuf_start[i] + j;
                int offset = gd->sbuf_start[i];
                suNf *mat_cpu = _4FIELD_AT_PTR(sendbuf_cpu, idx, mu, 0);
                suNf *in_block = _4FIELD_AT_PTR(sendbuf_gpu, offset, 0, 0);
                int stride = gd->sbuf_len[i];
                read_gpu_suNf(stride, (*mat_gpu), in_block, j, mu, 4);

                for (int comp = 0; comp < NF*NF; ++comp) {
                    diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                    diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
                }
            }     
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_gfield_f(f);
    return return_val;
}

int test_sync_identical_to_cpu_gfield_flt() 
{
    lprintf("INFO", 0, " ======= TEST SINGLE PRECISION GAUGE FIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNg_field_flt *f = alloc_gfield_flt(gd);
    random_gfield_flt_cpu(f);
    copy_to_gpu_gfield_flt(f);

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
    sync_field(f->type, 4*sizeof(*f->ptr), 0, f->ptr, f->sendbuf_ptr);
    suNg_flt* sendbuf_cpu = (suNg_flt*)malloc(4*sendbuf_len*sizeof(suNg_flt));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, 4*sendbuf_len*sizeof(suNg_flt));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_gfield_flt(f);
    suNg_flt* sendbuf_gpu = (suNg_flt*)malloc(4*sendbuf_len*sizeof(suNg_flt));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, 4*sendbuf_len*sizeof(suNg_flt), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNg_flt *mat_gpu = (suNg_flt*)malloc(sizeof(suNg_flt));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            for (int mu = 0; mu < 4; ++mu) {
                int idx = gd->sbuf_start[i] + j;
                int offset = gd->sbuf_start[i];
                suNg_flt *mat_cpu = _4FIELD_AT_PTR(sendbuf_cpu, idx, mu, 0);
                suNg_flt *in_block = _4FIELD_AT_PTR(sendbuf_gpu, offset, 0, 0);
                int stride = gd->sbuf_len[i];
                read_gpu_suNg_flt(stride, (*mat_gpu), in_block, j, mu, 4);

                for (int comp = 0; comp < NG*NG; ++comp) {
                    diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                    diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
                }
            }     
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_gfield_flt(f);
    return return_val;
}

int test_sync_identical_to_cpu_gfield_f_flt() 
{
    lprintf("INFO", 0, " ======= TEST REPRESENTED SINGLE PRECISION GAUGE FIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNf_field_flt *f = alloc_gfield_f_flt(gd);
    random_gfield_f_flt_cpu(f);
    copy_to_gpu_gfield_f_flt(f);

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
    sync_field(f->type, 4*sizeof(*f->ptr), 0, f->ptr, f->sendbuf_ptr);
    suNf_flt* sendbuf_cpu = (suNf_flt*)malloc(4*sendbuf_len*sizeof(suNf_flt));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, 4*sendbuf_len*sizeof(suNf_flt));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_gfield_f_flt(f);
    suNf_flt* sendbuf_gpu = (suNf_flt*)malloc(4*sendbuf_len*sizeof(suNf_flt));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, 4*sendbuf_len*sizeof(suNf_flt), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNf_flt *mat_gpu = (suNf_flt*)malloc(sizeof(suNf_flt));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            for (int mu = 0; mu < 4; ++mu) {
                int idx = gd->sbuf_start[i] + j;
                int offset = gd->sbuf_start[i];
                suNf_flt *mat_cpu = _4FIELD_AT_PTR(sendbuf_cpu, idx, mu, 0);
                suNf_flt *in_block = _4FIELD_AT_PTR(sendbuf_gpu, offset, 0, 0);
                int stride = gd->sbuf_len[i];
                read_gpu_suNf_flt(stride, (*mat_gpu), in_block, j, mu, 4);

                for (int comp = 0; comp < NF*NF; ++comp) {
                    diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                    diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
                }
            }     
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_gfield_f_flt(f);
    return return_val;
}

int test_sync_identical_to_cpu_suNg_scalar_field() 
{
    lprintf("INFO", 0, " ======= TEST SU(NG) SCALAR FIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNg_scalar_field *f = alloc_suNg_scalar_field(gd);
    random_suNg_scalar_field_cpu(f);
    copy_to_gpu_suNg_scalar_field(f);

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
    sync_field(f->type, sizeof(*(f->ptr)), 0, f->ptr, f->sendbuf_ptr);
    suNg_vector* sendbuf_cpu = (suNg_vector*)malloc(sendbuf_len*sizeof(suNg_vector));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, sendbuf_len*sizeof(suNg_vector));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_suNg_scalar_field(f);
    suNg_vector* sendbuf_gpu = (suNg_vector*)malloc(sendbuf_len*sizeof(suNg_vector));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, sendbuf_len*sizeof(suNg_vector), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNg_vector *mat_gpu = (suNg_vector*)malloc(sizeof(suNg_vector));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            int idx = gd->sbuf_start[i] + j;
            int offset = gd->sbuf_start[i];
            suNg_vector *mat_cpu = _FIELD_AT_PTR(sendbuf_cpu, idx, 0);
            suNg_vector *in_block = _FIELD_AT_PTR(sendbuf_gpu, offset, 0);
            int stride = gd->sbuf_len[i];
            read_gpu_suNg_vector(stride, (*mat_gpu), in_block, j, 0, 1);

            for (int comp = 0; comp < NG; ++comp) {
                diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
            } 
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_suNg_scalar_field(f);
    return return_val;
}

int test_sync_identical_to_cpu_avfield() 
{
    lprintf("INFO", 0, " ======= TEST AVFIELD ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNg_av_field *f = alloc_avfield(gd);
    random_avfield_cpu(f);
    copy_to_gpu_avfield(f);

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
    sync_field(f->type, 4*sizeof(*(f->ptr)), 0, f->ptr, f->sendbuf_ptr);
    suNg_algebra_vector* sendbuf_cpu = (suNg_algebra_vector*)malloc(4*sendbuf_len*sizeof(suNg_algebra_vector));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, 4*sendbuf_len*sizeof(suNg_algebra_vector));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_avfield(f);
    suNg_algebra_vector* sendbuf_gpu = (suNg_algebra_vector*)malloc(4*sendbuf_len*sizeof(suNg_algebra_vector));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, 4*sendbuf_len*sizeof(suNg_algebra_vector), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNg_algebra_vector *mat_gpu = (suNg_algebra_vector*)malloc(sizeof(suNg_algebra_vector));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            int idx = gd->sbuf_start[i] + j;
            int offset = gd->sbuf_start[i];
            suNg_algebra_vector *mat_cpu = _4FIELD_AT_PTR(sendbuf_cpu, idx, 0, 0);
            suNg_algebra_vector *in_block = _4FIELD_AT_PTR(sendbuf_gpu, offset, 0, 0);
            int stride = gd->sbuf_len[i];
            read_gpu_suNg_algebra_vector(stride, (*mat_gpu), in_block, j, 0, 4);

            for (int comp = 0; comp < NG*NG-1; ++comp) {
                diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
            } 
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_avfield(f);
    return return_val;
}

int test_sync_identical_to_cpu_gtransf() 
{
    lprintf("INFO", 0, " ======= TEST GTRANSF ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNg_field *f = alloc_gtransf(gd);
    random_gtransf_cpu(f);
    copy_to_gpu_gtransf(f);

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
    sync_field(f->type, sizeof(*(f->ptr)), 0, f->ptr, f->sendbuf_ptr);
    suNg* sendbuf_cpu = (suNg*)malloc(sendbuf_len*sizeof(suNg));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, sendbuf_len*sizeof(suNg));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_gtransf(f);
    suNg* sendbuf_gpu = (suNg*)malloc(sendbuf_len*sizeof(suNg));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, sendbuf_len*sizeof(suNg), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNg *mat_gpu = (suNg*)malloc(sizeof(suNg));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            int idx = gd->sbuf_start[i] + j;
            int offset = gd->sbuf_start[i];
            suNg *mat_cpu = _FIELD_AT_PTR(sendbuf_cpu, idx, 0);
            suNg *in_block = _FIELD_AT_PTR(sendbuf_gpu, offset, 0);
            int stride = gd->sbuf_len[i];
            read_gpu_suNg(stride, (*mat_gpu), in_block, j, 0, 1);

            for (int comp = 0; comp < NG*NG; ++comp) {
                diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
            } 
        }
    }
    return_val += check_diff_norm_zero(diff);
    free_gtransf(f);
    return return_val;
}

int test_sync_identical_to_cpu_clover_ldl() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER LDL ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    ldl_field *f = alloc_clover_ldl(gd);
    random_clover_ldl_cpu(f);
    copy_to_gpu_clover_ldl(f);

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
    sync_field(f->type, sizeof(*(f->ptr)), 0, f->ptr, f->sendbuf_ptr);
    ldl_t* sendbuf_cpu = (ldl_t*)malloc(sendbuf_len*sizeof(ldl_t));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, sendbuf_len*sizeof(ldl_t));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_clover_ldl(f);
    ldl_t* sendbuf_gpu = (ldl_t*)malloc(sendbuf_len*sizeof(ldl_t));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, sendbuf_len*sizeof(ldl_t), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    ldl_t *mat_gpu = (ldl_t*)malloc(sizeof(ldl_t));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            int idx = gd->sbuf_start[i] + j;
            int offset = gd->sbuf_start[i];
            ldl_t *mat_cpu = _FIELD_AT_PTR(sendbuf_cpu, idx, 0);
            ldl_t *in_block = _FIELD_AT_PTR(sendbuf_gpu, offset, 0);
            int stride = gd->sbuf_len[i];
            read_gpu_ldl_t(stride, (*mat_gpu), in_block, j, 0, 1);

            for (int comp = 0; comp < NF*(2*NF-1); ++comp) {
                diff += creal((*mat_cpu).up[comp])-creal((*mat_gpu).up[comp]);
                diff += cimag((*mat_cpu).up[comp])-cimag((*mat_gpu).up[comp]);
                diff += creal((*mat_cpu).dn[comp])-creal((*mat_gpu).dn[comp]);
                diff += cimag((*mat_cpu).dn[comp])-cimag((*mat_gpu).dn[comp]);
            } 
        }
    }
    return_val += check_diff_norm_zero(diff);
    free_clover_ldl(f);
    return return_val;
}

int test_sync_identical_to_cpu_clover_term() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER TERM ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNfc_field *f = alloc_clover_term(gd);
    random_clover_term_cpu(f);
    copy_to_gpu_clover_term(f);

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
    sync_field(f->type, 4*sizeof(*f->ptr), 0, f->ptr, f->sendbuf_ptr);
    suNfc* sendbuf_cpu = (suNfc*)malloc(4*sendbuf_len*sizeof(suNfc));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, 4*sendbuf_len*sizeof(suNfc));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_clover_term(f);
    suNfc* sendbuf_gpu = (suNfc*)malloc(4*sendbuf_len*sizeof(suNfc));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, 4*sendbuf_len*sizeof(suNfc), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNfc *mat_gpu = (suNfc*)malloc(sizeof(suNfc));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            for (int mu = 0; mu < 4; ++mu) {
                int idx = gd->sbuf_start[i] + j;
                int offset = gd->sbuf_start[i];
                suNfc *mat_cpu = _4FIELD_AT_PTR(sendbuf_cpu, idx, mu, 0);
                suNfc *in_block = _4FIELD_AT_PTR(sendbuf_gpu, offset, 0, 0);
                int stride = gd->sbuf_len[i];
                read_gpu_suNfc(stride, (*mat_gpu), in_block, j, mu, 4);

                for (int comp = 0; comp < NG*NG; ++comp) {
                    diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                    diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
                }
            }     
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_clover_term(f);
    return return_val;
}

int test_sync_identical_to_cpu_clover_force() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER FORCE ======= \n");
    lprintf("INFO", 0, "Check sendbuffers filled by sync are identical on CPU and GPU\n");

    // Setup fields on GPU
    geometry_descriptor* gd = &glattice;
    int return_val = 0;
    suNf_field *f = alloc_clover_force(gd);
    random_clover_force_cpu(f);
    copy_to_gpu_clover_force(f);

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
    sync_field(f->type, 6*sizeof(*f->ptr), 0, f->ptr, f->sendbuf_ptr);
    suNf* sendbuf_cpu = (suNf*)malloc(6*sendbuf_len*sizeof(suNf));
    memcpy(sendbuf_cpu, f->sendbuf_ptr, 6*sendbuf_len*sizeof(suNf));

    // Sync to buffer on GPU and save the sendbuffer in another array
    sync_gpu_clover_force(f);
    suNf* sendbuf_gpu = (suNf*)malloc(6*sendbuf_len*sizeof(suNf));
    cudaMemcpy(sendbuf_gpu, f->sendbuf_gpu_ptr, 6*sendbuf_len*sizeof(suNf), cudaMemcpyDeviceToHost);

    // Iterate over the arrays and check
    double diff = 0.0;
    suNf *mat_gpu = (suNf*)malloc(sizeof(suNf));
    for (int i = 0; i < gd->nbuffers_spinor; ++i) 
    {
        for (int j = 0; j < gd->sbuf_len[i]; ++j) 
        {
            for (int mu = 0; mu < 4; ++mu) {
                int idx = gd->sbuf_start[i] + j;
                int offset = gd->sbuf_start[i];
                suNf *mat_cpu = _6FIELD_AT_PTR(sendbuf_cpu, idx, mu, 0);
                suNf *in_block = _6FIELD_AT_PTR(sendbuf_gpu, offset, 0, 0);
                int stride = gd->sbuf_len[i];
                read_gpu_suNf(stride, (*mat_gpu), in_block, j, mu, 6);

                for (int comp = 0; comp < NG*NG; ++comp) {
                    diff += creal((*mat_cpu).c[comp])-creal((*mat_gpu).c[comp]);
                    diff += cimag((*mat_cpu).c[comp])-cimag((*mat_gpu).c[comp]);
                }
            }     
        }
    }

    return_val += check_diff_norm_zero(diff);
    free_clover_force(f);
    return return_val;
}