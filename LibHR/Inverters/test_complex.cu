#ifdef WITH_GPU
#include "inverters.h"

int test_overload_plus_rhs_double() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = c + b;
    if (abs(4.4 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(2.2 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b) > pow(10, -6)) { result += int(pow(10, 4)); }
    return result;
}

__global__ void test_gpu1(double b, hr_complex c, hr_complex *a) {
    *a = c + b;
}

int test_overload_plus_rhs_double_gpu() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    hr_complex *a_d;
    double b = 3.3;
    cudaMalloc((void **)&a_d, sizeof(hr_complex));

    test_gpu1<<<1, 1, 0, 0>>>(b, c, a_d);
    cudaDeviceSynchronize();
    cudaMemcpy(&a, a_d, sizeof(hr_complex), cudaMemcpyDeviceToHost);
    cudaFree(a_d);
    if (abs(4.4 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(2.2 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b) > pow(10, -6)) { result += int(pow(10, 4)); }
    return result;
}

int test_overload_plus_lhs_double() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = b + c;
    if (abs(4.4 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(2.2 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b) > pow(10, -6)) { result += int(pow(10, 4)); }
    return result;
}

int test_overload_prod_rhs_double() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = c * b;
    if (abs(3.63 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(7.26 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b) > pow(10, -6)) { result += int(pow(10, 4)); }
    return result;
}

int test_overload_prod_lhs_double() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = b * c;
    if (abs(3.63 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(7.26 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b) > pow(10, -6)) { result += int(pow(10, 4)); }
    return result;
}

int test_overload_div_rhs_double() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = c / b;
    if (abs(1.0 / 3.0 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(2.0 / 3.0 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b) > pow(10, -6)) { result += int(pow(10, 4)); }
    return result;
}

int test_overload_div_lhs_double() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = b / c;
    if (abs(0.6 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(-1.2 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b) > pow(10, -6)) { result += int(pow(10, 4)); }
    return result;
}

int test_overload_plus_hr_complex() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    a = b + c;
    if (abs(4.4 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(6.6 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4.4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

int test_overload_minus_hr_complex() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    a = b - c;
    if (abs(2.2 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(2.2 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4.4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

int test_overload_prod_hr_complex() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    a = b * c;
    if (abs(-6.05 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(12.1 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4.4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

__global__ void test_gpu9(hr_complex b, hr_complex c, hr_complex *a) {
    *a = b * c;
}

int test_overload_prod_hr_complex_gpu() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    hr_complex *a_d;
    cudaMalloc((void **)&a_d, sizeof(hr_complex));

    test_gpu9<<<1, 1, 0, 0>>>(b, c, a_d);
    cudaDeviceSynchronize();
    cudaMemcpy(&a, a_d, sizeof(hr_complex), cudaMemcpyDeviceToHost);
    cudaFree(a_d);
    if (abs(-6.05 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(12.1 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4.4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

int test_overload_div_hr_complex() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    a = b / c;
    if (abs(2.2 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(-0.4 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3.3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4.4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

int test_negate() {
    int result = 0;
    hr_complex b = hr_complex(1.1, 2.2);
    hr_complex a;
    a = -b;
    if (abs(-1.1 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(-2.2 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - b.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - b.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    return result;
}

int test_cast_double() {
    int result = 0;
    hr_complex a;
    double b = 1.1;
    a = (hr_complex)b;
    if (abs(1.1 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(0.0 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - b) > pow(10, -6)) { result += int(pow(10, 2)); }
    return result;
}

int test_overload_plus_rhs_integer() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    hr_complex_int b = hr_complex_int(3, 4);
    a = c + b;
    if (abs(4.1 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(6.2 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

int test_overload_plus_lhs_integer() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    hr_complex_int b = hr_complex_int(3, 4);
    a = b + c;
    if (abs(4.1 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(6.2 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

int test_overload_div_rhs_integer() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex_int b = hr_complex_int(3, 4);
    hr_complex a;
    a = c / b;
    if (abs(0.484 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(0.088 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

int test_overload_div_lhs_integer() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex_int b = hr_complex_int(3, 4);
    hr_complex a;
    a = b / c;
    if (abs(2.0 - a.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(-4.0 / 11.0 - a.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    if (abs(3 - b.re) > pow(10, -6)) { result += int(pow(10, 4)); }
    if (abs(4 - b.im) > pow(10, -6)) { result += int(pow(10, 5)); }
    return result;
}

int test_I_add() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b;
    b = c + I;
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(1.1 - b.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(3.2 - b.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    return result;
}

int test_I_prod() {
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b;
    b = c * I;
    if (abs(1.1 - c.re) > pow(10, -6)) { result += int(pow(10, 0)); }
    if (abs(2.2 - c.im) > pow(10, -6)) { result += int(pow(10, 1)); }
    if (abs(-2.2 - b.re) > pow(10, -6)) { result += int(pow(10, 2)); }
    if (abs(1.1 - b.im) > pow(10, -6)) { result += int(pow(10, 3)); }
    return result;
}

#endif //WITH_GPU
