/*
* NOCOMPILE = !WITH_GPU
*/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int test1;
    int gpu_test1;
    int test2;
    int test3;
    int test4;
    int test5;
    int test6;
    int test7;
    int test8;
    int test9;
    int gpu_test9;
    int test10;
    int test11;
    int test12;
    int test13;
    int test14;
    int test15;
    int test16;
    int test17;
    int test18;

    test1 = test_overload_plus_rhs_double();
    gpu_test1 = test_overload_plus_rhs_double_gpu();
    test2 = test_overload_plus_lhs_double();
    test3 = test_overload_prod_rhs_double();
    test4 = test_overload_prod_lhs_double();
    test5 = test_overload_div_rhs_double();
    test6 = test_overload_div_lhs_double();
    test7 = test_overload_plus_hr_complex();
    test8 = test_overload_minus_hr_complex();
    test9 = test_overload_prod_hr_complex();
    gpu_test9 = test_overload_prod_hr_complex_gpu();
    test10 = test_overload_div_hr_complex();
    test11 = test_negate();
    test12 = test_cast_double();
    test13 = test_overload_plus_rhs_integer();
    test14 = test_overload_plus_lhs_integer();
    test15 = test_overload_div_rhs_integer();
    test16 = test_overload_div_lhs_integer();
    test17 = test_I_add();
    test18 = test_I_prod();

    printf("test1: %d\n", test1);
    printf("gpu_test1: %d\n", gpu_test1);
    printf("test2: %d\n", test2);
    printf("test3: %d\n", test3);
    printf("test4: %d\n", test4);
    printf("test5: %d\n", test5);
    printf("test6: %d\n", test6);
    printf("test7: %d\n", test7);
    printf("test8: %d\n", test8);
    printf("test9: %d\n", test9);
    printf("gpu_test9: %d\n", gpu_test9);
    printf("test10: %d\n", test10);
    printf("test11: %d\n", test11);
    printf("test12: %d\n", test12);
    printf("test13: %d\n", test13);
    printf("test14: %d\n", test14);
    printf("test15: %d\n", test15);
    printf("test16: %d\n", test16);
    printf("test17: %d\n", test17);
    printf("test18: %d\n", test18);
}
