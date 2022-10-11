#include <math.h>
#include "logger.h"

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