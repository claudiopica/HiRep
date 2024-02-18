/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Test of C++ complex number implementation
*
*******************************************************************************/

#include "libhr.h"


int main(int argc, char *argv[]) {
    int return_value = 0;
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    return_value += test_overload_plus_rhs_double();
    return_value += test_overload_plus_lhs_double();
    return_value += test_overload_prod_rhs_double();
    return_value += test_overload_prod_lhs_double();
    return_value += test_overload_div_rhs_double();
    return_value += test_overload_div_lhs_double();
    return_value += test_overload_plus_hr_complex();
    return_value += test_overload_minus_hr_complex();
    return_value += test_overload_prod_hr_complex();
    return_value += test_overload_div_hr_complex();
    return_value += test_negate();
    return_value += test_cast_double();
    return_value += test_overload_plus_rhs_integer();
    return_value += test_overload_plus_lhs_integer();
    return_value += test_overload_div_rhs_integer();
    return_value += test_overload_div_lhs_integer();
    return_value += test_I_add();
    return_value += test_I_prod();

    return_value += check_plus();
    return_value += check_minus();
    return_value += check_mul();
    return_value += check_div();

    if(return_value == 0) lprintf("RESULT", 0, "PASS");
    else lprintf("RESULT", 0, "FAIL");
    finalize_process();
    return return_value;
}