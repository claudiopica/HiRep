/*******************************************************************************
 *
 * NOCOMPILE= WITH_GPU
 *
 * Speed test of dirac operator for single precision
 *
 *******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[])
{
  spinor_field_flt *s0, *s1;
  float elapsed;
  int flopsite, bytesite;
  int n_reps, n_reps_trial;
  int n_warmup = 1000;
  struct timeval start, end, etime;

  setup_process(&argc, &argv);

  setup_gauge_fields();

  lprintf("MAIN", 0, "n_warmup = %d\n", n_warmup);

  u_gauge_f_flt = (suNf_field_flt *)u_gauge_flt;

  /* allocate memory */
  lprintf("MAIN", 0, "Allocating spinor field\n");
  s0 = alloc_spinor_field_f_flt(2, &glattice);
  s1 = s0 + 1;

  lprintf("MAIN", 0, "Randomizing spinor field...\n");
  gaussian_spinor_field_flt(s0);
  gaussian_spinor_field_flt(s1);

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  fflush(stdout);
  random_u(u_gauge);
  // assign_ud2u();
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
  assign_ud2u_f();

  lprintf("MAIN", 0, "done.\n");

  // Check speed diracoperator

#if defined(REPR_ADJOINT)
  flopsite = 8 * NF * (7 + 8 * NF);
#else
  flopsite = 8 * NF * (7 + 16 * NF);
#endif
  bytesite = 36 * sizeof(suNf_vector_flt) + 16 * sizeof(suNf_flt); // add integers for geometry indexes?
  bytesite = 40 * sizeof(suNf_vector) + 8 * sizeof(suNf);          // 8 spinors read + 1 spinor read+write + 8 gauge matrices read

  lprintf("LA TEST", 0, "Flop per site = %d\n", flopsite);
  lprintf("LA TEST", 0, "Byte per site = %d\n", bytesite);
  lprintf("LA TEST", 0, "Dirac data movement = %d bytes\n", bytesite * VOLUME);

  /// speed test Dirac operator
  lprintf("LA TEST", 0, "Warmup application of the Diracoperator %d times.\n", n_warmup);
  elapsed = 0;
  gettimeofday(&start, 0);
  _OMP_PRAGMA(_omp_parallel)
  {
    for (int i = 0; i < n_warmup; ++i)
    {
      Dphi_flt_(s1, s0);
    }
  }
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  elapsed = etime.tv_sec * 1000. + etime.tv_usec * 0.001;
  n_reps = (int)((double)(n_warmup * 200) / elapsed);

  lprintf("LA TEST", 0, "\nEvaluating the massless Diracoperator.\n");
  elapsed = 0.;
  n_reps /= 10;
  n_reps_trial = n_reps;
  do
  {

    gettimeofday(&start, 0);

    for (int i = 0; i < n_reps_trial; ++i)
    {
      Dphi_flt_(s1, s0);
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

  free_spinor_field_f_flt(s0);

  finalize_process();
  exit(0);
}
