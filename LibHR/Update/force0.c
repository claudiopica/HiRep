/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"

#include <stdio.h>
#include <math.h>

#define _print_avect(a) printf("(%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f)\n", (a).c1, (a).c2, (a).c3, (a).c4, (a).c5, (a).c6, (a).c7, (a).c8)

#define _print_mat(a) printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",creal((a).c1_1),creal((a).c1_2),creal((a).c1_3),creal((a).c2_1),creal((a).c2_2),creal((a).c2_3),creal((a).c3_1),creal((a).c3_2),(a).c3_3)); \
  printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n", cimag((a).c1_1), cimag((a).c1_2), cimag((a).c1_3), cimag((a).c2_1), cimag((a).c2_2), cimag((a).c2_3), cimag((a).c3_1), cimag((a).c3_2), cimag((a).c3_3))

void force0(double dt, void *vpar)
{
#ifdef TIMING
  struct timeval start, end;
  struct timeval etime;

#ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
#endif
  gettimeofday(&start, 0);
#endif

  force_gauge_par *par = (force_gauge_par *)vpar;
  suNg_av_field *force = *par->momenta;
  double coeff = -dt * par->beta / NG;

  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge, *par->momenta);

#ifdef MEASURE_FORCE0
  double forcestat0 = 0.; /* used for computation of avr and max force */
  double forcestat1 = 0.; /* used for computation of avr and max force */

  _MASTER_FOR_SUM(&glattice, i, forcestat0, forcestat1)
  {
#else
  _MASTER_FOR(&glattice, i)
  {
#endif
    suNg s1, s2;
    suNg_algebra_vector f;
    for (int mu = 0; mu < 4; ++mu)
    {
      staples(i, mu, &s1);
      _suNg_times_suNg_dagger(s2, *_4FIELD_AT(u_gauge, i, mu), s1);

      /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
      _fund_algebra_project(f, s2);
      _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force, i, mu), coeff, f);

#ifdef MEASURE_FORCE0
      double nsq;
      _algebra_vector_sqnorm_g(nsq, f);
      forcestat0 += sqrt(nsq);
      for (int x = 0; x < sizeof(suNg_algebra_vector) / sizeof(double); ++x)
      {
        if (forcestat1 < fabs(f.c[x]))
          forcestat1 = fabs(f.c[x]);
      }
#endif
    }
  }

#ifdef MEASURE_FORCE0
  //  if(logger_getlevel("FORCE-STAT")>=10){
  global_sum(&forcestat0, 1);
  global_max(&forcestat1, 1);

  forcestat0 *= beta / ((4. * NG) * GLB_VOLUME);
  forcestat1 *= beta / ((double)(NG));
  //    lprintf("FORCE-STAT",10," force0 : dt= %1.8e avr |force|= %1.8e maxforce= %1.8e \n",dt,forcestat[0],forcestat[1]);
  lprintf("FORCE_STAT", 20, "GF: avr dt |force| = %1.8e dt maxforce = %1.8e, dt = %1.8e \n", forcestat0 * dt, forcestat1 * dt, dt);
  force_ave[0] += dt * forcestat0;
  force_max[0] += dt * forcestat1;
  //  }
#endif

  apply_BCs_on_momentum_field(force);

#ifdef TIMING
#ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
#endif
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("TIMING", 0, "force0 %.6f s\n", 1. * etime.tv_sec + 1.e-6 * etime.tv_usec);
#endif
}
