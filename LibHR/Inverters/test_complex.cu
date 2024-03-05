#ifdef WITH_GPU
#ifndef WITH_MPI
#include "inverters.h"
#include "utils.h"
#include "io.h"
#include "random.h"

static int nreps = 10;

// CPU Consistency tests
int test_overload_plus_rhs_double() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = c + b;

    result += check_diff_norm(4.4 - a.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b, EPSILON_TEST);
    return result;
}

int test_overload_plus_lhs_double() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = b + c;

    result += check_diff_norm(4.4 - a.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b, EPSILON_TEST);
    return result;
}

int test_overload_prod_rhs_double() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = c * b;

    result += check_diff_norm(3.63 - a.re, EPSILON_TEST);
    result += check_diff_norm(7.26 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b, EPSILON_TEST);
    return result;
}

int test_overload_prod_lhs_double() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = b * c;

    result += check_diff_norm(3.63 - a.re, EPSILON_TEST);
    result += check_diff_norm(7.26 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b, EPSILON_TEST);
    return result;
}

int test_overload_div_rhs_double() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = c / b;

    result += check_diff_norm(1.0 / 3.0 - a.re, EPSILON_TEST);
    result += check_diff_norm(2.0 / 3.0 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b, EPSILON_TEST);
    return result;
}

int test_overload_div_lhs_double() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    double b = 3.3;
    a = b / c;

    result += check_diff_norm(0.6 - a.re, EPSILON_TEST);
    result += check_diff_norm(-1.2 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b, EPSILON_TEST);
    return result;
}

int test_overload_plus_hr_complex() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    a = b + c;

    result += check_diff_norm(4.4 - a.re, EPSILON_TEST);
    result += check_diff_norm(6.6 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b.re, EPSILON_TEST);
    result += check_diff_norm(4.4 - b.im, EPSILON_TEST);

    return result;
}

int test_overload_minus_hr_complex() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    a = b - c;

    result += check_diff_norm(2.2 - a.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b.re, EPSILON_TEST);
    result += check_diff_norm(4.4 - b.im, EPSILON_TEST);
    return result;
}

int test_overload_prod_hr_complex() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    a = b * c;

    result += check_diff_norm(-6.05 - a.re, EPSILON_TEST);
    result += check_diff_norm(12.1 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b.re, EPSILON_TEST);
    result += check_diff_norm(4.4 - b.im, EPSILON_TEST);
    return result;
}

int test_overload_div_hr_complex() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b = hr_complex(3.3, 4.4);
    hr_complex a;
    a = b / c;

    result += check_diff_norm(2.2 - a.re, EPSILON_TEST);
    result += check_diff_norm(-0.4 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3.3 - b.re, EPSILON_TEST);
    result += check_diff_norm(4.4 - b.im, EPSILON_TEST);
    return result;
}

int test_negate() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex b = hr_complex(1.1, 2.2);
    hr_complex a;
    a = -b;

    result += check_diff_norm(-1.1 - a.re, EPSILON_TEST);
    result += check_diff_norm(-2.2 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - b.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - b.im, EPSILON_TEST);
    return result;
}

int test_cast_double() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex a;
    double b = 1.1;
    a = (hr_complex)b;

    result += check_diff_norm(1.1 - a.re, EPSILON_TEST);
    result += check_diff_norm(0.0 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - b, EPSILON_TEST);
    return result;
}

int test_overload_plus_rhs_integer() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    hr_complex_int b = hr_complex_int(3, 4);
    a = c + b;

    result += check_diff_norm(4.1 - a.re, EPSILON_TEST);
    result += check_diff_norm(6.2 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3 - b.re, EPSILON_TEST);
    result += check_diff_norm(4 - b.im, EPSILON_TEST);
    return result;
}

int test_overload_plus_lhs_integer() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex a;
    hr_complex_int b = hr_complex_int(3, 4);
    a = b + c;

    result += check_diff_norm(4.1 - a.re, EPSILON_TEST);
    result += check_diff_norm(6.2 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3 - b.re, EPSILON_TEST);
    result += check_diff_norm(4 - b.im, EPSILON_TEST);
    return result;
}

int test_overload_div_rhs_integer() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex_int b = hr_complex_int(3, 4);
    hr_complex a;
    a = c / b;

    result += check_diff_norm(0.484 - a.re, EPSILON_TEST);
    result += check_diff_norm(0.088 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3 - b.re, EPSILON_TEST);
    result += check_diff_norm(4 - b.im, EPSILON_TEST);
    return result;
}

int test_overload_div_lhs_integer() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex_int b = hr_complex_int(3, 4);
    hr_complex a;
    a = b / c;

    result += check_diff_norm(2.0 - a.re, EPSILON_TEST);
    result += check_diff_norm(-4.0 / 11.0 - a.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(3 - b.re, EPSILON_TEST);
    result += check_diff_norm(4 - b.im, EPSILON_TEST);
    return result;
}

int test_I_add() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b;
    b = c + I;

    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(1.1 - b.re, EPSILON_TEST);
    result += check_diff_norm(3.2 - b.im, EPSILON_TEST);
    return result;
}

int test_I_prod() {
    lprintf("TEST", 0, __func__);
    int result = 0;
    hr_complex c = hr_complex(1.1, 2.2);
    hr_complex b;
    b = c * I;

    result += check_diff_norm(1.1 - c.re, EPSILON_TEST);
    result += check_diff_norm(2.2 - c.im, EPSILON_TEST);
    result += check_diff_norm(-2.2 - b.re, EPSILON_TEST);
    result += check_diff_norm(1.1 - b.im, EPSILON_TEST);
    return result;
}

#define MAXVAL 10

template <class L> void random_number(L &c) {
    double rand_no[1];
    rand_no[0] = 0;
    while (rand_no[0] == 0) {
        ranlxd(rand_no, 1);
    }
    c = (rand_no[0] * (L)MAXVAL - (L)(MAXVAL / 2)) * (L)MAXVAL;
}

template <class L> void random_number(hr_complex_t<L> &c) {
    double rand_no[2];
    rand_no[0] = 0.0;
    rand_no[1] = 0.0;
    while (rand_no[0] == 0 || rand_no[1] == 0) {
        ranlxd(rand_no, 2);
        rand_no[0] = (L)(rand_no[0] * (L)MAXVAL - (L)(MAXVAL / 2)) * (L)MAXVAL;
        rand_no[1] = (L)(rand_no[0] * (L)MAXVAL - (L)(MAXVAL / 2)) * (L)MAXVAL;
    }
    c = hr_complex_t<L>(rand_no[0], rand_no[1]);
}

#define CHECK_KERNEL(__operator__, __fname__)                                              \
    template <class L, class T, class M> __global__ void __fname__##_gpu(L b, T c, M *a) { \
        *a = c __operator__ b;                                                             \
    }

#define safe_diff(__res__, __a__, __b__)        \
    if (__a__ != 0) {                           \
        __res__ = abs((__a__ - __b__) / __a__); \
    } else {                                    \
        __res__ = abs(__a__ - __b__);           \
    }

#define CHECK_OPERATOR(__operator__, __fname__, __debugprint__)                       \
    template <typename CPX1, typename CPX2> int __fname__##_tmpl(double prec) {       \
        lprintf("TEST", 0, __func__);                                                 \
        lprintf("TEST", 0, "\n");                                                     \
        CPX1 a, a_gpu, *a_d, c, c_gpu;                                                \
        int result = 0;                                                               \
        CPX2 b, b_gpu;                                                                \
        cudaMalloc((void **)&a_d, sizeof(CPX1));                                      \
        double res_vec[12 * nreps];                                                   \
                                                                                      \
        for (int i = 0; i < nreps; i++) {                                             \
            random_number(c);                                                         \
            random_number(b);                                                         \
            b_gpu = b;                                                                \
            c_gpu = c;                                                                \
                                                                                      \
            a = c __operator__ b;                                                     \
            __fname__##_gpu<<<1, 1>>>(b_gpu, c_gpu, a_d);                             \
                                                                                      \
            cudaDeviceSynchronize();                                                  \
            cudaMemcpy(&a_gpu, a_d, sizeof(CPX1), cudaMemcpyDeviceToHost);            \
                                                                                      \
            auto b2 = hr_complex_t(b);                                                \
            auto b_gpu2 = hr_complex_t(b_gpu);                                        \
            safe_diff(res_vec[6 * i], a.re, a_gpu.re);                                \
            safe_diff(res_vec[6 * i + 1], a.im, a_gpu.im);                            \
            safe_diff(res_vec[6 * i + 2], c.re, c_gpu.re);                            \
            safe_diff(res_vec[6 * i + 3], c.im, c_gpu.im);                            \
            safe_diff(res_vec[6 * i + 4], b2.re, b_gpu2.re);                          \
            safe_diff(res_vec[6 * i + 5], b2.im, b_gpu2.im);                          \
        }                                                                             \
        for (int i = 0; i < nreps; i++) {                                             \
            random_number(c);                                                         \
            random_number(b);                                                         \
            b_gpu = b;                                                                \
            c_gpu = c;                                                                \
                                                                                      \
            a = b __operator__ c;                                                     \
            __fname__##_gpu<<<1, 1>>>(c_gpu, b_gpu, a_d);                             \
                                                                                      \
            cudaDeviceSynchronize();                                                  \
            cudaMemcpy(&a_gpu, a_d, sizeof(CPX1), cudaMemcpyDeviceToHost);            \
                                                                                      \
            auto b2 = hr_complex_t(b);                                                \
            auto b_gpu2 = hr_complex_t(b_gpu);                                        \
            safe_diff(res_vec[6 * nreps + 6 * i], a.re, a_gpu.re);                    \
            safe_diff(res_vec[6 * nreps + 6 * i + 1], a.im, a_gpu.im);                \
            safe_diff(res_vec[6 * nreps + 6 * i + 2], c.re, c_gpu.re);                \
            safe_diff(res_vec[6 * nreps + 6 * i + 3], c.im, c_gpu.im);                \
            safe_diff(res_vec[6 * nreps + 6 * i + 4], b2.re, b_gpu2.re);              \
            safe_diff(res_vec[6 * nreps + 6 * i + 5], b2.im, b_gpu2.im);              \
        }                                                                             \
                                                                                      \
        double max = res_vec[0];                                                      \
        for (int i = 1; i < 12 * nreps; i++) {                                        \
            max = res_vec[i] > max ? res_vec[i] : max;                                \
        }                                                                             \
        result = check_diff_norm(max, prec);                                          \
        cudaFree(a_d);                                                                \
        return result;                                                                \
    }                                                                                 \
                                                                                      \
    int __fname__() {                                                                 \
        int result = 0;                                                               \
        lprintf("TEST", 0, "Check hr_complex -- double\n");                           \
        result += __fname__##_tmpl<hr_complex, double>(EPSILON_TEST);                 \
        lprintf("TEST", 0, "Check hr_complex -- float\n");                            \
        result += __fname__##_tmpl<hr_complex, float>(EPSILON_TEST);                  \
        lprintf("TEST", 0, "Check hr_complex -- int\n");                              \
        result += __fname__##_tmpl<hr_complex, int>(EPSILON_TEST);                    \
        lprintf("TEST", 0, "Check hr_complex_flt -- double\n");                       \
        result += __fname__##_tmpl<hr_complex_flt, double>(EPSILON_FLT_TEST);         \
        lprintf("TEST", 0, "Check hr_complex_flt -- float\n");                        \
        result += __fname__##_tmpl<hr_complex_flt, float>(EPSILON_FLT_TEST);          \
        lprintf("TEST", 0, "Check hr_complex_flt -- int\n");                          \
        result += __fname__##_tmpl<hr_complex_flt, int>(EPSILON_FLT_TEST);            \
        lprintf("TEST", 0, "Check hr_complex_int -- double\n");                       \
        result += __fname__##_tmpl<hr_complex_int, double>(EPSILON_TEST);             \
        lprintf("TEST", 0, "Check hr_complex_int -- float\n");                        \
        result += __fname__##_tmpl<hr_complex_int, float>(EPSILON_TEST);              \
        lprintf("TEST", 0, "Check hr_complex_int -- int\n");                          \
        result += __fname__##_tmpl<hr_complex_int, int>(EPSILON_TEST);                \
        lprintf("TEST", 0, "Check hr_complex -- hr_complex\n");                       \
        result += __fname__##_tmpl<hr_complex, hr_complex>(EPSILON_TEST);             \
        lprintf("TEST", 0, "Check hr_complex -- hr_complex_flt\n");                   \
        result += __fname__##_tmpl<hr_complex, hr_complex_flt>(EPSILON_FLT_TEST);     \
        lprintf("TEST", 0, "Check hr_complex -- hr_complex_int\n");                   \
        result += __fname__##_tmpl<hr_complex, hr_complex_int>(EPSILON_TEST);         \
        lprintf("TEST", 0, "Check hr_complex_flt -- hr_complex\n");                   \
        result += __fname__##_tmpl<hr_complex_flt, hr_complex>(EPSILON_FLT_TEST);     \
        lprintf("TEST", 0, "Check hr_complex_flt -- hr_complex_flt\n");               \
        result += __fname__##_tmpl<hr_complex_flt, hr_complex_flt>(EPSILON_FLT_TEST); \
        lprintf("TEST", 0, "Check hr_complex_flt -- hr_complex_int\n");               \
        result += __fname__##_tmpl<hr_complex_flt, hr_complex_flt>(EPSILON_FLT_TEST); \
        lprintf("TEST", 0, "Check hr_complex_int -- hr_complex\n");                   \
        result += __fname__##_tmpl<hr_complex_int, hr_complex>(EPSILON_TEST);         \
        lprintf("TEST", 0, "Check hr_complex_int -- hr_complex_flt\n");               \
        result += __fname__##_tmpl<hr_complex_int, hr_complex_flt>(EPSILON_FLT_TEST); \
        lprintf("TEST", 0, "Check hr_complex_int -- hr_complex_int\n");               \
        result += __fname__##_tmpl<hr_complex_int, hr_complex_int>(0.0);              \
        return result;                                                                \
    }

CHECK_KERNEL(+, check_plus)
CHECK_OPERATOR(+, check_plus, 0)
CHECK_KERNEL(-, check_minus)
CHECK_OPERATOR(-, check_minus, 0)
CHECK_KERNEL(*, check_mul)
CHECK_OPERATOR(*, check_mul, 1)
CHECK_KERNEL(/, check_div)
CHECK_OPERATOR(/, check_div, 0)

#endif //WITH_MPI
#endif //WITH_GPU
