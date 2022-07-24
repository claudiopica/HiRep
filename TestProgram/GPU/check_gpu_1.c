/*******************************************************************************
*
* NOCOMPILE= WITH_MPI
*
* Test of modules
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "global.h"
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "setup.h"
#include "logger.h"

static double EPSILON=1.e-12;
void random_spinor_field_cpu(spinor_field *s)
{
    geometry_descriptor *type = s->type;
    _PIECE_FOR(type, ixp){
        int start = type->master_start[ixp];
        int N = type->master_end[ixp] - type->master_start[ixp]+1;
        gauss((double*)(_FIELD_AT(s, start)), N * sizeof(suNf_spinor) / sizeof(double));
    }
}


void unit_array(double *a, int len)
{
    for(int i=0;i<len;i++){
        a[i]=1.;
    }
}


void unit_spinor_field_cpu(spinor_field *s)
{
    geometry_descriptor *type = s->type;
    _PIECE_FOR(type, ixp){
        int start = type->master_start[ixp];
        int N = type->master_end[ixp] - type->master_start[ixp]+1;
        unit_array((double*)(_FIELD_AT(s, start)), N * sizeof(suNf_spinor) / sizeof(double));
    }
}


int main(int argc, char *argv[])
{
    spinor_field *sf1, *sf2, *sf3;
    int sfsize = 1;
    double norm_cpu, norm_cpu2, norm_cpu3;

    /* setup process id and communications */
    //logger_setlevel(0,10000);
    logger_map("DEBUG", "debug");

    setup_process(&argc, &argv);

    // Allocate memory for CPU and GPU spinor fields
    sf1=alloc_spinor_field_f(sfsize, &glattice);
    sf2=alloc_spinor_field_f(sfsize, &glattice);
    sf3=alloc_spinor_field_f(sfsize, &glattice);

    // Assign random field to CPU
    for (int i=0;i<sfsize;i++){
        random_spinor_field_cpu(&sf1[i]);
        //unit_spinor_field_cpu(&sf1[i]);
    }
    // Copy spinor field to GPU and back to CPU
    for (int i=0;i<sfsize;i++){
        spinor_field_copy_f_cpu(&sf2[i],&sf1[i]);
               //print_spinor_field_cpu(&sf2[i]);
        spinor_field_copy_to_gpu_f(&sf2[i]);
        spinor_field_zero_f_cpu(&sf2[i]);
               //print_spinor_field_cpu(&sf2[i]);
        spinor_field_copy_f(&sf3[i],&sf2[i]);
        spinor_field_copy_from_gpu_f(&sf2[i]);
               //print_spinor_field_cpu(&sf2[i]);
        spinor_field_copy_from_gpu_f(&sf3[i]);
    }

    // Calculate norm on CPU
    norm_cpu = spinor_field_sqnorm_f_cpu(&sf1[0]);
    lprintf("GPU TEST",1,"Norm CPU 1: %.2e\n",norm_cpu);
    norm_cpu2 = spinor_field_sqnorm_f_cpu(&sf2[0]);
    lprintf("GPU TEST",1,"Norm CPU 2: %.2e\n",norm_cpu2);
    norm_cpu3 = spinor_field_sqnorm_f_cpu(&sf3[0]);
    lprintf("GPU TEST",1,"Norm CPU 3: %.2e\n",norm_cpu3);
    double d1=norm_cpu-norm_cpu2;
    lprintf("GPU TEST",1,"diff 1: %.2e\n",d1);
    double d2=norm_cpu-norm_cpu3;
    lprintf("GPU TEST",1,"diff 2: %.2e\n",d2);
}
