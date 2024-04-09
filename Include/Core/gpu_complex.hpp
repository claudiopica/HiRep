/***************************************************************************\
* Copyright (c) 2008-2024, Sofie Martins, Erik Kjellgren, Claudio Pica      *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file gpu_complex.h
 * @brief Type definitions and macros for complex numbers used in C++ and CUDA
 */

#ifndef GPU_COMPLEX_H
#define GPU_COMPLEX_H

#include "Core/gpu.h"
#define PI 3.141592653589793238462643383279502884197

template <class T> struct __align__(sizeof(T)) hr_complex_t {
    T re, im;

    /** CONSTRUCTORS **/
    constexpr visible hr_complex_t(const T a = (T)0, const T b = (T)0)
        : re(a)
        , im(b) {
    }

    template <class L>
    constexpr visible hr_complex_t(const hr_complex_t<L> a)
        : re(a.re)
        , im(a.im) {
    }

    /** ASSIGN **/
    template <class L> hr_complex_t visible inline __attribute__((always_inline)) &operator=(const L x) {
        re = x;
        im = 0;
        return *this;
    }

    template <class L> hr_complex_t visible inline __attribute__((always_inline)) &operator=(const hr_complex_t<L> x) {
        re = x.re;
        im = x.im;
        return *this;
    }

    hr_complex_t visible inline __attribute__((always_inline)) operator-(void) const {
        return hr_complex_t(-re, -im);
    }

    /** ADD ASSIGN **/
    template <class L> hr_complex_t visible inline __attribute__((always_inline)) &operator+=(const L x) {
        re += x;
        return *this;
    }

    template <class L> hr_complex_t visible inline __attribute__((always_inline)) &operator+=(const hr_complex_t<L> x) {
        re += x.re;
        im += x.im;
        return *this;
    }

    /** SUB ASSIGN **/
    template <class L> hr_complex_t visible inline __attribute__((always_inline)) &operator-=(const L x) {
        re -= x;
        return *this;
    }

    template <class L> auto visible inline __attribute__((always_inline)) &operator-=(const hr_complex_t<L> x) {
        re -= x.re;
        im -= x.im;
        return *this;
    }

    /** MUL ASSIGN **/
    template <class L> auto visible inline __attribute__((always_inline)) &operator*=(const L x) {
        re *= x;
        im *= x;
        return *this;
    }

    template <class L> auto visible inline __attribute__((always_inline)) &operator*=(const hr_complex_t<L> x) {
        T re_new = re * x.re - im * x.im;
        T im_new = re * x.im + im * x.re;
        re = re_new;
        im = im_new;
        return *this;
    }

    /** DIV ASSIGN **/
    template <class L> auto visible inline __attribute__((always_inline)) &operator/=(const L x) {
        re /= x;
        im /= x;
        return *this;
    }

    template <class L> auto visible inline __attribute__((always_inline)) &operator/=(const hr_complex_t<L> x) {
        T re_new = (re * x.re + im * x.im) / (x.re * x.re + x.im * x.im);
        T im_new = (im * x.re - re * x.im) / (x.re * x.re + x.im * x.im);
        re = re_new;
        im = im_new;
        return *this;
    }

    /** COMPLEX CONJUGATION **/
    hr_complex_t visible inline __attribute__((always_inline)) conj() {
        return hr_complex_t(re, -im);
    }
};

/** ADD **/
template <class L, class T> auto visible inline __attribute__((always_inline)) operator+(const hr_complex_t<L> &c, const T x) {
    return hr_complex_t(c.re + x, c.im + (T)0);
}

template <class L, class T> auto visible inline __attribute__((always_inline)) operator+(const L x, const hr_complex_t<T> &c) {
    return hr_complex_t(c.re + x, c.im + (L)0);
}

template <class L, class T>
auto visible inline __attribute__((always_inline)) operator+(const hr_complex_t<L> &x, const hr_complex_t<T> &c) {
    return hr_complex_t(c.re + x.re, c.im + x.im);
}

/** SUBTRACT **/
template <class T, class L> auto visible inline __attribute__((always_inline)) operator-(const hr_complex_t<T> &c, const L x) {
    return hr_complex_t(c.re - x, c.im - (L)0);
}

template <class T, class L> auto visible inline __attribute__((always_inline)) operator-(const T x, const hr_complex_t<L> &c) {
    return hr_complex_t(x - c.re, (T)0 - c.im);
}

template <class T, class L>
auto visible inline __attribute__((always_inline)) operator-(const hr_complex_t<L> &x, const hr_complex_t<T> &c) {
    return hr_complex_t(x.re - c.re, x.im - c.im);
}

/** MULTIPLY **/
template <class T, class L> auto visible inline __attribute__((always_inline)) operator*(const hr_complex_t<T> &c, const L x) {
    return hr_complex_t(c.re * x, c.im * x);
}

template <class T, class L> auto visible inline __attribute__((always_inline)) operator*(const T x, const hr_complex_t<L> &c) {
    return hr_complex_t(c.re * x, c.im * x);
}

template <class T, class L>
auto visible inline __attribute__((always_inline)) operator*(const hr_complex_t<T> &a, const hr_complex_t<L> &b) {
    return hr_complex_t(a.re * b.re - b.im * a.im, a.im * b.re + b.im * a.re);
}

/** DIVIDE **/
template <class T, class L> auto visible inline __attribute__((always_inline)) operator/(const hr_complex_t<T> &c, const L x) {
    return hr_complex_t(c.re / x, c.im / x);
}

template <class T, class L> auto visible inline __attribute__((always_inline)) operator/(const T x, const hr_complex_t<L> &c) {
    return hr_complex_t(x * c.re / (c.re * c.re + c.im * c.im), -x * c.im / (c.re * c.re + c.im * c.im));
}

template <class T, class L>
auto visible inline __attribute__((always_inline)) operator/(const hr_complex_t<T> x, const hr_complex_t<L> c) {
    return hr_complex_t((x.re * c.re + x.im * c.im) / (c.re * c.re + c.im * c.im),
                        (x.im * c.re - x.re * c.im) / (c.re * c.re + c.im * c.im));
}

typedef struct hr_complex_t<int> hr_complex_int;
typedef struct hr_complex_t<double> hr_complex;
typedef struct hr_complex_t<float> hr_complex_flt;

#define I (hr_complex_int(0, 1))
#define creal(a) ((a).re)
#define cimag(a) ((a).im)
#define conj(a) ((a).conj())

visible double carg(hr_complex c);

visible hr_complex cpow(hr_complex c, double pow);

#endif
