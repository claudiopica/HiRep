/*******************************************************************************
 *
 * Gauge covariance of the Dirac operator
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "setup.h"
#include "io.h"
#include "update.h"
#include "geometry.h"
#include "global.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "representation.h"
#include "utils.h"
#include "communications.h"

int main(int argc, char *argv[])
{

  double res1, res2, res3;
  spinor_field *s0, *s1, *s2;
  float elapsed;
  int flopsite, bytesite;
  int n_reps, n_reps_trial;
  int n_warmup = 1000;
  struct timeval start, end, etime;

  setup_process(&argc, &argv);

  setup_gauge_fields();

  lprintf("MAIN", 0, "n_warmup = %d\n", n_warmup);

  /* allocate memory */
  lprintf("MAIN", 0, "Allocating spinor field\n");
  s0 = alloc_spinor_field_f(1, &glattice);
  s1 = alloc_spinor_field_f(2, &glattice);
  s2 = s1 + 1;

  lprintf("MAIN", 0, "Randomizing spinor field...\n");
  gaussian_spinor_field(s0);
  gaussian_spinor_field(s1);

  //
  //#pragma omp parallel num_threads(1) default(shared)
  _OMP_PRAGMA(_omp_parallel num_threads(1))
  {
    res1 = spinor_field_sqnorm_f(s0);

    _OMP_BARRIER

    res2 = spinor_field_sqnorm_f(s1);
  }

  lprintf("LA_TEST", 0, "Square norm of the random spinor fields %lf and %lf\nThey must be different from zero\n", res1, res2);

  lprintf("MAIN", 0, "Generating a random gauge field... ");

  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();

  lprintf("MAIN", 0, "done.\n");

  // Check speed diracoperator

#if defined(REPR_ADJOINT)
  flopsite = 8 * NF * (7 + 8 * NF);
#else
  flopsite = 8 * NF * (7 + 16 * NF);
#endif
  bytesite = 36 * sizeof(suNf_vector) + 8 * sizeof(suNf); // 8 spinors read + 1 spinor write + 8 gauge matrices read
  bytesite = 40 * sizeof(suNf_vector) + 8 * sizeof(suNf); // 8 spinors read + 1 spinor read+write + 8 gauge matrices read

  lprintf("LA TEST", 0, "Flop per site = %d\n", flopsite);
  lprintf("LA TEST", 0, "Byte per site = %d\n", bytesite);
  lprintf("LA TEST", 0, "Dirac data movement = %d bytes\n", bytesite * VOLUME);

  // speed test Dirac operator
  lprintf("LA TEST", 0, "Warmup application of the fused Diracoperator %d times.\n", n_warmup);
  elapsed = 0;
  gettimeofday(&start, 0);
  _OMP_PRAGMA(_omp_parallel)
  {
    for (int i = 0; i < n_warmup; ++i)
    {
      Dphi_fused_(s1, s0);
    }
  }
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  elapsed = etime.tv_sec * 1000. + etime.tv_usec * 0.001;
  n_reps = (int)((double)(n_warmup * 200) / elapsed);

  lprintf("LA TEST", 0, "\nEvaluating the massless fused Diracoperator.\n");

  elapsed = 0;
  n_reps_trial = n_reps;
  bcast_int(&n_reps_trial, 1);

  do
  {
    gettimeofday(&start, 0);

    _OMP_PRAGMA(_omp_parallel)
    {
      for (int i = 0; i < n_reps_trial; ++i)
      {
        Dphi_fused_(s1, s0);
      }
    }
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    elapsed = etime.tv_sec * 1000. + etime.tv_usec * 0.001;
    n_reps = n_reps_trial;
    n_reps_trial = (int)((double)(n_reps_trial * 2200) / elapsed);
    bcast_int(&n_reps_trial, 1);
  } while (elapsed < 2000);

  lprintf("LA TEST", 0, "Massless fused Diracoperator reps: %d , data: %lf kb, time: %lf msec, GFLOPS: %1.6g , BAND: %1.6g GB/s\n",
          n_reps,
          ((double)(bytesite * VOLUME)) / 1024,
          elapsed,
          ((double)n_reps * VOLUME * flopsite) / elapsed / 1.e6,
          ((double)n_reps * VOLUME * bytesite / elapsed / 1.e6));

  lprintf("LA TEST", 0, "\nEvaluating the massless Diracoperator.\n");
  elapsed = 0.;
  n_reps /= 10;
  n_reps_trial = n_reps;
  do
  {

    gettimeofday(&start, 0);

    for (int i = 0; i < n_reps_trial; ++i)
    {
      Dphi_(s1, s0);
    }

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    elapsed = etime.tv_sec * 1000. + etime.tv_usec * 0.001;
    n_reps = n_reps_trial;
    n_reps_trial = (int)((double)(n_reps_trial * 2200) / elapsed);
    bcast_int(&n_reps_trial, 1);
  } while (elapsed < 2000);

  lprintf("LA TEST", 0, "Massless Diracoperator reps: %d , data: %lf kb, time: %lf msec, GFLOPS: %1.6g , BAND: %1.6g GB/s\n",
          n_reps,
          ((double)(bytesite * VOLUME)) / 1024,
          elapsed,
          ((double)n_reps * VOLUME * flopsite) / elapsed / 1.e6,
          ((double)n_reps * VOLUME * bytesite) / elapsed / 1.e6);

  _OMP_PRAGMA(_omp_parallel)
  {
    Dphi_fused_(s1, s0);
  }

  Dphi_(s2, s0);

  spinor_field_sub_assign_f(s2, s1);

  res1 = spinor_field_sqnorm_f(s2);

  res2 = spinor_field_sqnorm_f(s1);

  res3 = spinor_field_sqnorm_f(s0);

  lprintf("LA_TEST", 0, "\nCheck consistency of the two Dirac formulations: %.12e (should be zero or so)\nNorm Dphi input %.12e and output %.12e \n", res1, res3, res2);

  free_spinor_field_f(s0);
  free_spinor_field_f(s1);

  finalize_process();
  return 0;
}
