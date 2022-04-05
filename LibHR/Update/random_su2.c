/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
 *
 * File random_su2.c
 *
 * Random SU(2) matrices functions
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random.h"
#include "update.h"
#include "hr_omp.h"

#define NVEC (32)
#define NRAN (2 * NVEC)
#define PI_HALF 1.5707964f
#define PI 3.1415927f
#define TWO_PI 6.2831854f

static int *i_vec = NULL, *i_y, *i_v;
static double **vec1, **vec2, **vec3;
static double **r, **u, **v, **y;

static void init_su2_rand()
{
   int ntd = 1;
#ifdef _OPENMP
   ntd = omp_get_num_threads();
   _OMP_PRAGMA(master)
   {
      lprintf("init_su2_rand", 0, "Init rand for threads, performing %d init", ntd);
   }
#endif
   _OMP_PRAGMA(single)
   {

      i_vec = malloc(sizeof(int) * ntd);
      i_y = malloc(sizeof(int) * ntd);
      i_v = malloc(sizeof(int) * ntd);

      vec1 = malloc(sizeof(double *) * ntd);
      vec2 = malloc(sizeof(double *) * ntd);
      vec3 = malloc(sizeof(double *) * ntd);
      r = malloc(sizeof(double *) * ntd);
      u = malloc(sizeof(double *) * ntd);
      v = malloc(sizeof(double *) * ntd);
      y = malloc(sizeof(double *) * ntd);

      for (int i = 0; i < ntd; i++)
      {
         i_vec[i] = NVEC;
         i_y[i] = NRAN;
         i_v[i] = NRAN;

         vec1[i] = malloc(sizeof(double) * NVEC);
         vec2[i] = malloc(sizeof(double) * NVEC);
         vec3[i] = malloc(sizeof(double) * NVEC);
         r[i] = malloc(sizeof(double) * NRAN);
         u[i] = malloc(sizeof(double) * NRAN);
         v[i] = malloc(sizeof(double) * NRAN);
         y[i] = malloc(sizeof(double) * NRAN);

      }
   }
   _OMP_BARRIER
}
static void update_vec(int tid)
{
   int i;

   double r1, r2, rsq;

   ranlxd(r[tid], NRAN);

   for (i = 0; i < NVEC; i++)
   {
      r1 = 2.0 * r[tid][i] - 1.0;
      r2 = TWO_PI * r[tid][NVEC + i] - PI;
      rsq = sqrt(1.0 - r1 * r1);

      vec1[tid][i] = r1;
      vec2[tid][i] = rsq * sin(r2);
      vec3[tid][i] = rsq * cos(r2);
   }

   i_vec[tid] = 0;
}

static void update_y(int tid)
{
   int i;
   double r1, r2, r3, r4, s, c;

   ranlxd(y[tid], NRAN);
   ranlxd(u[tid], NRAN);
   ranlxd(r[tid], NRAN);

   for (i = 0; i < NVEC; i++)
   {
      r1 = -log(1.0 - y[tid][i]);
      r2 = PI_HALF * y[tid][NVEC + i];
      r3 = log(1.0 - u[tid][i]);
      r4 = log(1.0 - u[tid][NVEC + i]);

      s = sin(r2);
      s *= s;
      c = 1.0 - s;

      y[tid][i] = r1 * s - r3;
      y[tid][NVEC + i] = r1 * c - r4;

      r1 = r[tid][i] * r[tid][i];
      r2 = r[tid][NVEC + i] * r[tid][NVEC + i];
      u[tid][i] = r1 + r1;
      u[tid][NVEC + i] = r2 + r2;
   }

   i_y[tid] = 0;
}

void random_su2(double rho, double s[])
/*
 *  Computes a random vector s[4] with probability density
 *  proportional to exp(rho*s[0])*delta(1-s^2) assuming rho>=0
 */
{
   double rhoinv, s0p1, ut, rt;
   double s0, s1, s2, s3, sq;
   int tid = 0;
#ifdef _OPENMP

   tid = omp_get_thread_num();
#endif
   if (i_vec == NULL)
   {
      _OMP_BARRIER
      init_su2_rand();
   }
   if (i_vec[tid] == NVEC)
      update_vec(tid);

   if (rho > 1.5)
   {
      rhoinv = 1.0 / rho;

      for (;;)
      {
         if (i_y[tid] == NRAN)
            update_y(tid);

         s0p1 = 2.0 - rhoinv * y[tid][i_y[tid]];
         ut = u[tid][i_y[tid]++];

         if (ut <= s0p1)
            break;
      }
   }
   else if (rho > 0.3)
   {
      rhoinv = 1.0 / rho;
      rt = exp(rho + rho) - 1.0;

      for (;;)
      {
         if (i_v[tid] == NRAN)
         {
            ranlxd(v[tid], NRAN);
            i_v[tid] = 0;
         }
         s0p1 = rhoinv * log(1.0 + rt * v[tid][i_v[tid]++]);
         ut = v[tid][i_v[tid]++];

         if ((ut * ut) <= (s0p1 * (2.0 - s0p1)))
            break;
      }
   }
   else
   {
      for (;;)
      {
         if (i_v[tid] == NRAN)
         {
            ranlxd(v[tid], NRAN);
            i_v[tid] = 0;
         }

         s0p1 = 2.0 * v[tid][i_v[tid]++];
         rt = exp(rho * (s0p1 - 2.0));
         ut = v[tid][i_v[tid]++];

         if ((ut * ut) <= (s0p1 * (2.0 - s0p1) * rt * rt))
            break;
      }
   }

   sq = sqrt(s0p1 * (2.0 - s0p1));
   s0 = s0p1 - 1.0;
   s1 = sq * vec1[tid][i_vec[tid]];
   s2 = sq * vec2[tid][i_vec[tid]];
   s3 = sq * vec3[tid][i_vec[tid]];

   sq = s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3;
   sq = 1.5 - 0.5 * sq;

   s[0] = sq * s0;
   s[1] = sq * s1;
   s[2] = sq * s2;
   s[3] = sq * s3;

   i_vec[tid] += 1;
}
