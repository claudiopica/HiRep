/* This is part of the ranlux 3.3 distribution by Martin Luesher */

/*******************************************************************************
 *
 * File ranlxd.c
 *
 * Copyright (C) 2005 Martin Luescher
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Random number generator "ranlxd". See the notes
 *
 *   "User's guide for ranlxs and ranlxd v3.2" (December 2005)
 *
 *   "Algorithms used in ranlux v3.0" (May 2001)
 *
 * for a detailed description
 *
 * The externally accessible functions are
 *
 *   void ranlxd(double r[],int n)
 *     Computes the next n double-precision random numbers and
 *     assigns them to the elements r[0],...,r[n-1] of the array r[]
 *
 *   void rlxd_init(int level,int seed)
 *     Initialization of the generator
 *
 *   int rlxd_size(void)
 *     Returns the number of integers required to save the state of
 *     the generator
 *
 *   void rlxd_get(int state[])
 *     Extracts the current state of the generator and stores the
 *     information in the array state[N] where N>=rlxd_size()
 *
 *   void rlxd_reset(int state[])
 *     Resets the generator to the state defined by the array state[N]
 *
 *******************************************************************************/

#define RANLXD_C

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ranlux.h"
#include "global.h"

static void local_error(int no)
{
    switch (no)
    {
    case 1:
        printf("Error in subroutine rlxd_init\n");
        printf("Bad choice of luxury level (should be 1 or 2)\n");
        break;
    case 2:
        printf("Error in subroutine rlxd_init\n");
        printf("Bad choice of seed (should be between 1 and 2^31-1)\n");
        break;
    case 3:
        printf("Error in rlxd_get\n");
        printf("Undefined state (ranlxd is not initialized\n");
        break;
    case 5:
        printf("Error in rlxd_reset\n");
        printf("Unexpected input data\n");
        break;
    case 6:
        printf("Unitialized random seed\n");
        break;
    }
    printf("Program aborted\n");
    exit(0);
}

#if (defined SSE)

typedef struct
{
    float c1, c2, c3, c4;
} vec_t __attribute__((aligned(16)));

typedef struct
{
    vec_t c1, c2;
} dble_vec_t __attribute__((aligned(16)));

static int init = 0, pr, prm, ir, jr, is, is_old, next[96];
static vec_t one, one_bit, carry;

static union
{
    dble_vec_t vec[12];
    float num[96];
} x __attribute__((aligned(16)));

#define STEP(pi, pj)                                       \
    __asm__ __volatile__("movaps %4, %%xmm4 \n\t"          \
                         "movaps %%xmm2, %%xmm3 \n\t"      \
                         "subps %2, %%xmm4 \n\t"           \
                         "movaps %%xmm1, %%xmm5 \n\t"      \
                         "cmpps $0x6, %%xmm4, %%xmm2 \n\t" \
                         "andps %%xmm2, %%xmm5 \n\t"       \
                         "subps %%xmm3, %%xmm4 \n\t"       \
                         "andps %%xmm0, %%xmm2 \n\t"       \
                         "addps %%xmm4, %%xmm5 \n\t"       \
                         "movaps %%xmm5, %0 \n\t"          \
                         "movaps %5, %%xmm6 \n\t"          \
                         "movaps %%xmm2, %%xmm3 \n\t"      \
                         "subps %3, %%xmm6 \n\t"           \
                         "movaps %%xmm1, %%xmm7 \n\t"      \
                         "cmpps $0x6, %%xmm6, %%xmm2 \n\t" \
                         "andps %%xmm2, %%xmm7 \n\t"       \
                         "subps %%xmm3, %%xmm6 \n\t"       \
                         "andps %%xmm0, %%xmm2 \n\t"       \
                         "addps %%xmm6, %%xmm7 \n\t"       \
                         "movaps %%xmm7, %1"               \
                         : "=m"((*pi).c1),                 \
                           "=m"((*pi).c2)                  \
                         : "m"((*pi).c1),                  \
                           "m"((*pi).c2),                  \
                           "m"((*pj).c1),                  \
                           "m"((*pj).c2)                   \
                         : "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7")

static void update(void)
{
    int k, kmax;
    dble_vec_t *pmin, *pmax, *pi, *pj;

    kmax = pr;
    pmin = &x.vec[0];
    pmax = pmin + 12;
    pi = &x.vec[ir];
    pj = &x.vec[jr];

    __asm__ __volatile__("movaps %0, %%xmm0 \n\t"
                         "movaps %1, %%xmm1 \n\t"
                         "movaps %2, %%xmm2"
                         :
                         : "m"(one_bit),
                           "m"(one),
                           "m"(carry)
                         : "xmm0", "xmm1", "xmm2");

    for (k = 0; k < kmax; k++)
    {
        STEP(pi, pj);
        pi += 1;
        pj += 1;
        if (pi == pmax)
            pi = pmin;
        if (pj == pmax)
            pj = pmin;
    }

    __asm__ __volatile__("movaps %%xmm2, %0"
                         : "=m"(carry));

    ir += prm;
    jr += prm;
    if (ir >= 12)
        ir -= 12;
    if (jr >= 12)
        jr -= 12;
    is = 8 * ir;
    is_old = is;
}

static void define_constants(void)
{
    int k;
    float b;

    one.c1 = 1.0f;
    one.c2 = 1.0f;
    one.c3 = 1.0f;
    one.c4 = 1.0f;

    b = (float)(ldexp(1.0, -24));
    one_bit.c1 = b;
    one_bit.c2 = b;
    one_bit.c3 = b;
    one_bit.c4 = b;

    for (k = 0; k < 96; k++)
    {
        next[k] = (k + 1) % 96;
        if ((k % 4) == 3)
            next[k] = (k + 5) % 96;
    }
}

void rlxd_init(int level, int seed)
{
    int i, k, l;
    int ibit, jbit, xbit[31];
    int ix, iy;

    define_constants();

    if (level == 1)
        pr = 202;
    else if (level == 2)
        pr = 397;
    else
        local_error(1);

    i = seed + MPI_PID;

    for (k = 0; k < 31; k++)
    {
        xbit[k] = i % 2;
        i /= 2;
    }

    if ((seed <= 0) || (i != 0))
        local_error(2);

    ibit = 0;
    jbit = 18;

    for (i = 0; i < 4; i++)
    {
        for (k = 0; k < 24; k++)
        {
            ix = 0;

            for (l = 0; l < 24; l++)
            {
                iy = xbit[ibit];
                ix = 2 * ix + iy;

                xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2;
                ibit = (ibit + 1) % 31;
                jbit = (jbit + 1) % 31;
            }

            if ((k % 4) != i)
                ix = 16777215 - ix;

            x.num[4 * k + i] = (float)(ldexp((double)(ix), -24));
        }
    }

    carry.c1 = 0.0f;
    carry.c2 = 0.0f;
    carry.c3 = 0.0f;
    carry.c4 = 0.0f;

    ir = 0;
    jr = 7;
    is = 91;
    is_old = 0;
    prm = pr % 12;
    init = 1;
}

void ranlxd(double r[], int n)
{
    int k;

    if (init == 0)
        local_error(6);

    for (k = 0; k < n; k++)
    {
        is = next[is];
        if (is == is_old)
            update();
        r[k] = (double)(x.num[is + 4]) + (double)(one_bit.c1 * x.num[is]);
    }
}

int rlxd_size(void)
{
    return (105);
}

void rlxd_get(int state[])
{
    int k;
    float base;

    if (init == 0)
        local_error(3);

    base = (float)(ldexp(1.0, 24));
    state[0] = rlxd_size();

    for (k = 0; k < 96; k++)
        state[k + 1] = (int)(base * x.num[k]);

    state[97] = (int)(base * carry.c1);
    state[98] = (int)(base * carry.c2);
    state[99] = (int)(base * carry.c3);
    state[100] = (int)(base * carry.c4);

    state[101] = pr;
    state[102] = ir;
    state[103] = jr;
    state[104] = is;
}

void rlxd_reset(int state[])
{
    int k;

    define_constants();

    if (state[0] != rlxd_size())
        local_error(5);

    for (k = 0; k < 96; k++)
    {
        if ((state[k + 1] < 0) || (state[k + 1] >= 167777216))
            local_error(5);

        x.num[k] = (float)(ldexp((double)(state[k + 1]), -24));
    }

    if (((state[97] != 0) && (state[97] != 1)) ||
        ((state[98] != 0) && (state[98] != 1)) ||
        ((state[99] != 0) && (state[99] != 1)) ||
        ((state[100] != 0) && (state[100] != 1)))
        local_error(5);

    carry.c1 = (float)(ldexp((double)(state[97]), -24));
    carry.c2 = (float)(ldexp((double)(state[98]), -24));
    carry.c3 = (float)(ldexp((double)(state[99]), -24));
    carry.c4 = (float)(ldexp((double)(state[100]), -24));

    pr = state[101];
    ir = state[102];
    jr = state[103];
    is = state[104];
    is_old = 8 * ir;
    prm = pr % 12;
    init = 1;

    if (((pr != 202) && (pr != 397)) ||
        (ir < 0) || (ir > 11) || (jr < 0) || (jr > 11) || (jr != ((ir + 7) % 12)) ||
        (is < 0) || (is > 91))
        local_error(5);
}

#else
#ifdef _OPENMP
#include "hr_omp.h"

#define BASE 0x1000000
#define MASK 0xffffff

typedef struct
{
    int c1, c2, c3, c4;
} vec_t;

typedef struct
{
    vec_t c1, c2;
} dble_vec_t;

typedef union
{
    dble_vec_t vec[12];
    int num[96];
} x_struct;

static int *init = NULL, *pr, *prm, *ir, *jr, *is, *is_old, **next;
static double *one_bit;
static vec_t *carry;
static x_struct *x;

#define STEP(pi, pj, ltid)                          \
    d = (*pj).c1.c1 - (*pi).c1.c1 - carry[ltid].c1; \
    (*pi).c2.c1 += (d < 0);                         \
    d += BASE;                                      \
    (*pi).c1.c1 = d & MASK;                         \
    d = (*pj).c1.c2 - (*pi).c1.c2 - carry[ltid].c2; \
    (*pi).c2.c2 += (d < 0);                         \
    d += BASE;                                      \
    (*pi).c1.c2 = d & MASK;                         \
    d = (*pj).c1.c3 - (*pi).c1.c3 - carry[ltid].c3; \
    (*pi).c2.c3 += (d < 0);                         \
    d += BASE;                                      \
    (*pi).c1.c3 = d & MASK;                         \
    d = (*pj).c1.c4 - (*pi).c1.c4 - carry[ltid].c4; \
    (*pi).c2.c4 += (d < 0);                         \
    d += BASE;                                      \
    (*pi).c1.c4 = d & MASK;                         \
    d = (*pj).c2.c1 - (*pi).c2.c1;                  \
    carry[ltid].c1 = (d < 0);                       \
    d += BASE;                                      \
    (*pi).c2.c1 = d & MASK;                         \
    d = (*pj).c2.c2 - (*pi).c2.c2;                  \
    carry[ltid].c2 = (d < 0);                       \
    d += BASE;                                      \
    (*pi).c2.c2 = d & MASK;                         \
    d = (*pj).c2.c3 - (*pi).c2.c3;                  \
    carry[ltid].c3 = (d < 0);                       \
    d += BASE;                                      \
    (*pi).c2.c3 = d & MASK;                         \
    d = (*pj).c2.c4 - (*pi).c2.c4;                  \
    carry[ltid].c4 = (d < 0);                       \
    d += BASE;                                      \
    (*pi).c2.c4 = d & MASK

static void update(void)
{
    int k, kmax, d;
    dble_vec_t *pmin, *pmax, *pi, *pj;

    int tid = omp_get_thread_num();

    kmax = pr[tid];
    pmin = &(x[tid].vec[0]);
    pmax = pmin + 12;
    pi = &(x[tid].vec[ir[tid]]);
    pj = &(x[tid].vec[jr[tid]]);

    for (k = 0; k < kmax; k++)
    {
        STEP(pi, pj, tid);

        pi += 1;
        pj += 1;
        if (pi == pmax)
            pi = pmin;
        if (pj == pmax)
            pj = pmin;
    }

    ir[tid] += prm[tid];
    jr[tid] += prm[tid];
    if (ir[tid] >= 12)
        ir[tid] -= 12;
    if (jr[tid] >= 12)
        jr[tid] -= 12;
    is[tid] = 8 * ir[tid];
    is_old[tid] = is[tid];
}

static void define_constants(void)
{
    int k;
    int tid = omp_get_thread_num();

    one_bit[tid] = ldexp(1.0, -24);

    for (k = 0; k < 96; k++)
    {
        next[tid][k] = (k + 1) % 96;
        if ((k % 4) == 3)
            next[tid][k] = (k + 5) % 96;
    }
}

void rlxd_init(int level, int seed)
{
    _OMP_PRAGMA(_omp_parallel)
    {
        int tid = omp_get_thread_num();
        int nt = omp_get_num_threads();
        _OMP_PRAGMA(single)
        {

            init = malloc(nt * sizeof(int));
            pr = malloc(nt * sizeof(int));
            prm = malloc(nt * sizeof(int));
            ir = malloc(nt * sizeof(int));
            jr = malloc(nt * sizeof(int));
            is = malloc(nt * sizeof(int));
            is_old = malloc(nt * sizeof(int));
            next = malloc(nt * sizeof(int *));

            one_bit = malloc(nt * sizeof(double));
            carry = malloc(nt * sizeof(vec_t));
            x = malloc(nt * sizeof(x_struct));
        }
        next[tid] = malloc(96 * sizeof(int));
        int i, k, l;
        int ibit, jbit, xbit[31];
        int ix, iy;

        if ((INT_MAX < 2147483647) || (FLT_RADIX != 2) || (FLT_MANT_DIG < 24) ||
            (DBL_MANT_DIG < 48))
            local_error(0);

        define_constants();

        if (level == 1)
            pr[tid] = 202;
        else if (level == 2)
            pr[tid] = 397;
        else
            local_error(1);

        i = seed + tid + nt * MPI_PID;

        for (k = 0; k < 31; k++)
        {
            xbit[k] = i % 2;
            i /= 2;
        }

        if ((seed <= 0) || (i != 0))
            local_error(2);

        ibit = 0;
        jbit = 18;

        for (i = 0; i < 4; i++)
        {
            for (k = 0; k < 24; k++)
            {
                ix = 0;

                for (l = 0; l < 24; l++)
                {
                    iy = xbit[ibit];
                    ix = 2 * ix + iy;

                    xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2;
                    ibit = (ibit + 1) % 31;
                    jbit = (jbit + 1) % 31;
                }

                if ((k % 4) != i)
                    ix = 16777215 - ix;

                x[tid].num[4 * k + i] = ix;
            }
        }

        carry[tid].c1 = 0;
        carry[tid].c2 = 0;
        carry[tid].c3 = 0;
        carry[tid].c4 = 0;

        ir[tid] = 0;
        jr[tid] = 7;
        is[tid] = 91;
        is_old[tid] = 0;
        prm[tid] = pr[tid] % 12;
        init[tid] = 1;
    }
}

void ranlxd(double r[], int n)
{
    int k;
    int tid = omp_get_thread_num();

    if (init == NULL)
        local_error(6);

    for (k = 0; k < n; k++)
    {

        is[tid] = next[tid][is[tid]];

        if (is[tid] == is_old[tid])
            update();
        r[k] = one_bit[tid] * ((double)(x[tid].num[is[tid] + 4]) + one_bit[tid] * (double)(x[tid].num[is[tid]]));
    }
}

int rlxd_size(void)
{
    return (105);
}

void rlxd_get(int state[])
{
    int k;
    int tid = omp_get_thread_num();

    if (init[tid] == 0)
        local_error(3);

    state[0] = rlxd_size();

    for (k = 0; k < 96; k++)
        state[k + 1] = x[tid].num[k];

    state[97] = carry[tid].c1;
    state[98] = carry[tid].c2;
    state[99] = carry[tid].c3;
    state[100] = carry[tid].c4;

    state[101] = pr[tid];
    state[102] = ir[tid];
    state[103] = jr[tid];
    state[104] = is[tid];
}

void rlxd_reset(int state[])
{
    int k;
    int tid = omp_get_thread_num();
    if ((INT_MAX < 2147483647) || (FLT_RADIX != 2) || (FLT_MANT_DIG < 24) ||
        (DBL_MANT_DIG < 48))
        local_error(4);

    define_constants();

    if (state[0] != rlxd_size())
        local_error(5);

    for (k = 0; k < 96; k++)
    {
        if ((state[k + 1] < 0) || (state[k + 1] >= 167777216))
            local_error(5);

        x[tid].num[k] = state[k + 1];
    }

    if (((state[97] != 0) && (state[97] != 1)) ||
        ((state[98] != 0) && (state[98] != 1)) ||
        ((state[99] != 0) && (state[99] != 1)) ||
        ((state[100] != 0) && (state[100] != 1)))
        local_error(5);

    carry[tid].c1 = state[97];
    carry[tid].c2 = state[98];
    carry[tid].c3 = state[99];
    carry[tid].c4 = state[100];

    pr[tid] = state[101];
    ir[tid] = state[102];
    jr[tid] = state[103];
    is[tid] = state[104];
    is_old[tid] = 8 * ir[tid];
    prm[tid] = pr[tid] % 12;
    init[tid] = 1;

    if (((pr[tid] != 202) && (pr[tid] != 397)) ||
        (ir[tid] < 0) || (ir[tid] > 11) || (jr[tid] < 0) || (jr[tid] > 11) || (jr[tid] != ((ir[tid] + 7) % 12)) ||
        (is[tid] < 0) || (is[tid] > 91))
        local_error(5);
}
#else

#define BASE 0x1000000
#define MASK 0xffffff

typedef struct
{
    int c1, c2, c3, c4;
} vec_t;

typedef struct
{
    vec_t c1, c2;
} dble_vec_t;

static int init = 0, pr, prm, ir, jr, is, is_old, next[96];
static double one_bit;
static vec_t carry;

static union
{
    dble_vec_t vec[12];
    int num[96];
} x;

#define STEP(pi, pj)                          \
    d = (*pj).c1.c1 - (*pi).c1.c1 - carry.c1; \
    (*pi).c2.c1 += (d < 0);                   \
    d += BASE;                                \
    (*pi).c1.c1 = d & MASK;                   \
    d = (*pj).c1.c2 - (*pi).c1.c2 - carry.c2; \
    (*pi).c2.c2 += (d < 0);                   \
    d += BASE;                                \
    (*pi).c1.c2 = d & MASK;                   \
    d = (*pj).c1.c3 - (*pi).c1.c3 - carry.c3; \
    (*pi).c2.c3 += (d < 0);                   \
    d += BASE;                                \
    (*pi).c1.c3 = d & MASK;                   \
    d = (*pj).c1.c4 - (*pi).c1.c4 - carry.c4; \
    (*pi).c2.c4 += (d < 0);                   \
    d += BASE;                                \
    (*pi).c1.c4 = d & MASK;                   \
    d = (*pj).c2.c1 - (*pi).c2.c1;            \
    carry.c1 = (d < 0);                       \
    d += BASE;                                \
    (*pi).c2.c1 = d & MASK;                   \
    d = (*pj).c2.c2 - (*pi).c2.c2;            \
    carry.c2 = (d < 0);                       \
    d += BASE;                                \
    (*pi).c2.c2 = d & MASK;                   \
    d = (*pj).c2.c3 - (*pi).c2.c3;            \
    carry.c3 = (d < 0);                       \
    d += BASE;                                \
    (*pi).c2.c3 = d & MASK;                   \
    d = (*pj).c2.c4 - (*pi).c2.c4;            \
    carry.c4 = (d < 0);                       \
    d += BASE;                                \
    (*pi).c2.c4 = d & MASK

static void update(void)
{
    int k, kmax, d;
    dble_vec_t *pmin, *pmax, *pi, *pj;

    kmax = pr;
    pmin = &x.vec[0];
    pmax = pmin + 12;
    pi = &x.vec[ir];
    pj = &x.vec[jr];

    for (k = 0; k < kmax; k++)
    {
        STEP(pi, pj);
        pi += 1;
        pj += 1;
        if (pi == pmax)
            pi = pmin;
        if (pj == pmax)
            pj = pmin;
    }

    ir += prm;
    jr += prm;
    if (ir >= 12)
        ir -= 12;
    if (jr >= 12)
        jr -= 12;
    is = 8 * ir;
    is_old = is;
}

static void define_constants(void)
{
    int k;

    one_bit = ldexp(1.0, -24);

    for (k = 0; k < 96; k++)
    {
        next[k] = (k + 1) % 96;
        if ((k % 4) == 3)
            next[k] = (k + 5) % 96;
    }
}

void rlxd_init(int level, int seed)
{
    int i, k, l;
    int ibit, jbit, xbit[31];
    int ix, iy;

    if ((INT_MAX < 2147483647) || (FLT_RADIX != 2) || (FLT_MANT_DIG < 24) ||
        (DBL_MANT_DIG < 48))
        local_error(0);

    define_constants();

    if (level == 1)
        pr = 202;
    else if (level == 2)
        pr = 397;
    else
        local_error(1);

    i = seed + MPI_PID;

    for (k = 0; k < 31; k++)
    {
        xbit[k] = i % 2;
        i /= 2;
    }

    if ((seed <= 0) || (i != 0))
        local_error(2);

    ibit = 0;
    jbit = 18;

    for (i = 0; i < 4; i++)
    {
        for (k = 0; k < 24; k++)
        {
            ix = 0;

            for (l = 0; l < 24; l++)
            {
                iy = xbit[ibit];
                ix = 2 * ix + iy;

                xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2;
                ibit = (ibit + 1) % 31;
                jbit = (jbit + 1) % 31;
            }

            if ((k % 4) != i)
                ix = 16777215 - ix;

            x.num[4 * k + i] = ix;
        }
    }

    carry.c1 = 0;
    carry.c2 = 0;
    carry.c3 = 0;
    carry.c4 = 0;

    ir = 0;
    jr = 7;
    is = 91;
    is_old = 0;
    prm = pr % 12;
    init = 1;
}

void ranlxd(double r[], int n)
{
    int k;

    if (init == 0)
        local_error(6);

    for (k = 0; k < n; k++)
    {
        is = next[is];

        if (is == is_old)
            update();

        r[k] = one_bit * ((double)(x.num[is + 4]) + one_bit * (double)(x.num[is]));
    }
}

int rlxd_size(void)
{
    return (105);
}

void rlxd_get(int state[])
{
    int k;

    if (init == 0)
        local_error(3);

    state[0] = rlxd_size();

    for (k = 0; k < 96; k++)
        state[k + 1] = x.num[k];

    state[97] = carry.c1;
    state[98] = carry.c2;
    state[99] = carry.c3;
    state[100] = carry.c4;

    state[101] = pr;
    state[102] = ir;
    state[103] = jr;
    state[104] = is;
}

void rlxd_reset(int state[])
{
    int k;

    if ((INT_MAX < 2147483647) || (FLT_RADIX != 2) || (FLT_MANT_DIG < 24) ||
        (DBL_MANT_DIG < 48))
        local_error(4);

    define_constants();

    if (state[0] != rlxd_size())
        local_error(5);

    for (k = 0; k < 96; k++)
    {
        if ((state[k + 1] < 0) || (state[k + 1] >= 167777216))
            local_error(5);

        x.num[k] = state[k + 1];
    }

    if (((state[97] != 0) && (state[97] != 1)) ||
        ((state[98] != 0) && (state[98] != 1)) ||
        ((state[99] != 0) && (state[99] != 1)) ||
        ((state[100] != 0) && (state[100] != 1)))
        local_error(5);

    carry.c1 = state[97];
    carry.c2 = state[98];
    carry.c3 = state[99];
    carry.c4 = state[100];

    pr = state[101];
    ir = state[102];
    jr = state[103];
    is = state[104];
    is_old = 8 * ir;
    prm = pr % 12;
    init = 1;

    if (((pr != 202) && (pr != 397)) ||
        (ir < 0) || (ir > 11) || (jr < 0) || (jr > 11) || (jr != ((ir + 7) % 12)) ||
        (is < 0) || (is > 91))
        local_error(5);
}
#endif

#endif
