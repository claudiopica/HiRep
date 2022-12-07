#include <math.h>
#include "logger.h"
#include "ranlux.h"

void test_setup() 
{
    // TODO: other settings
    rlxd_init(1, 205);
    rlxs_init(2, 208);
}


int check_diff_norm(double diff_norm, double tol) 
{
    int return_val = 0;
    if (fabs(diff_norm) > tol) {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    } else 
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);
    return return_val;
}

int check_diff_norm_zero(double diff_norm) 
{
    int return_val = 0;
    if (fabs(diff_norm) != 0) 
    {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    } else 
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);
    return return_val;
}