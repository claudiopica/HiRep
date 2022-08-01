/******************************************************************************
*
* Test of hermiticity
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "setup.h"

#include "communications.h"

static double hmass = 0.1;

void MM_cpu(spinor_field *out, spinor_field *in)
{
#ifdef UPDATE_EO
  g5Dphi_eopre_sq_cpu(-hmass, out, in);
#else
  g5Dphi_sq_cpu(-hmass, out, in);
#endif
}


void DD_cpu(spinor_field *out, spinor_field *in)
{
  Dphi_cpu(hmass, out, in);
}


void HH_cpu(spinor_field *out, spinor_field *in)
{
  g5Dphi_cpu(-hmass, out, in);
}


void MM_gpu(spinor_field *out, spinor_field *in)
{
#ifdef UPDATE_EO
  lprintf("RESULT", 0, "In MM_gpu eopre\n");
  g5Dphi_eopre_sq(-hmass, out, in);
#else
  lprintf("RESULT", 0, "In MM_gpu\n");
  g5Dphi_sq(-hmass, out, in);
#endif
}


void DD_gpu(spinor_field *out, spinor_field *in)
{
  Dphi(hmass, out, in);
}


void HH_gpu(spinor_field *out, spinor_field *in)
{
  g5Dphi(-hmass, out, in);
}


int test_herm_cpu(spinor_operator S, const char *name, spinor_field *s1, spinor_field *s2, const int size)
{
  int return_val;
  double tau;
  spinor_field *s3, *s4;

#ifdef UPDATE_EO
  s3 = alloc_spinor_field_f(size, &glat_even);
  s4 = alloc_spinor_field_f(size, &glat_even);
#else
  s3 = alloc_spinor_field_f(size, &glattice);
  s4 = alloc_spinor_field_f(size, &glattice);
#endif
  for (int i=0; i < size; i++){
    S(&s3[i], &s1[i]);
    S(&s4[i], &s2[i]);

    lprintf("RESULT", 0, "s1 NORM %f on CPU\n", sqrt(spinor_field_sqnorm_f_cpu(&s1[i])));
    lprintf("RESULT", 0, "s2 NORM %f on CPU\n", sqrt(spinor_field_sqnorm_f_cpu(&s2[i])));
    lprintf("RESULT", 0, "s3 NORM %f on CPU\n", sqrt(spinor_field_sqnorm_f_cpu(&s3[i])));
    lprintf("RESULT", 0, "s4 NORM %f on CPU\n", sqrt(spinor_field_sqnorm_f_cpu(&s4[i])));
  }

  return_val = 0;
  return return_val;
}


int test_herm_gpu(spinor_operator S, const char *name, spinor_field *s1, spinor_field *s2, const int size)
{
  int return_val;
  double tau;
  spinor_field *s3, *s4;

#ifdef UPDATE_EO
  s3 = alloc_spinor_field_f(size, &glat_even);
  s4 = alloc_spinor_field_f(size, &glat_even);
#else
  s3 = alloc_spinor_field_f(size, &glattice);
  s4 = alloc_spinor_field_f(size, &glattice);
#endif
  for (int i=0; i < size; i++){
    spinor_field_copy_to_gpu_f(&s1[i]);
    spinor_field_copy_to_gpu_f(&s2[i]);
    spinor_field_copy_to_gpu_f(&s3[i]);
    spinor_field_copy_to_gpu_f(&s4[i]);

    S(&s3[i], &s1[i]);
    S(&s4[i], &s2[i]);

    lprintf("RESULT", 0, "s1 NORM %f on GPU\n", sqrt(spinor_field_sqnorm_f(&s1[i])));
    lprintf("RESULT", 0, "s2 NORM %f on GPU\n", sqrt(spinor_field_sqnorm_f(&s2[i])));
    lprintf("RESULT", 0, "s3 NORM %f on GPU\n", sqrt(spinor_field_sqnorm_f(&s3[i])));
    lprintf("RESULT", 0, "s4 NORM %f on GPU\n", sqrt(spinor_field_sqnorm_f(&s4[i])));
  }

  return_val = 0;
  return return_val;
}


int main(int argc, char *argv[])
{
  int return_value_cpu, return_value_gpu;
  const int sfsize = 1;
  spinor_field *sf1, *sf2, *sf3, *sf4;
  /* setup process id and communications */
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  random_u(u_gauge);
  lprintf("MAIN", 0, "done.\n");

  start_gf_sendrecv(u_gauge);

  represent_gauge_field();

#ifdef UPDATE_EO
  sf1 = alloc_spinor_field_f(sfsize, &glat_even);
  sf2 = alloc_spinor_field_f(sfsize, &glat_even);
#else
  sf1 = alloc_spinor_field_f(sfsize, &glattice);
  sf2 = alloc_spinor_field_f(sfsize, &glattice);
#endif
  gaussian_spinor_field(sf1);
  gaussian_spinor_field(sf2);
  return_value_cpu = test_herm_cpu(&MM_cpu, "M", sf1, sf2, sfsize);
  return_value_gpu = test_herm_gpu(&MM_gpu, "M", sf1, sf2, sfsize);
/*
#ifdef UPDATE_EO
  sf1 = alloc_spinor_field_f(sfsize, &glat_even);
  sf2 = alloc_spinor_field_f(sfsize, &glat_even);
#else
  sf1 = alloc_spinor_field_f(sfsize, &glattice);
  sf2 = alloc_spinor_field_f(sfsize, &glattice);
#endif
  gaussian_spinor_field(sf1);
  gaussian_spinor_field(sf2);
  return_value_cpu = test_herm_cpu(&DD_cpu, "D", sf1, sf2, sfsize);
  return_value_gpu = test_herm_gpu(&DD_gpu, "D", sf1, sf2, sfsize);

#ifdef UPDATE_EO
  sf1 = alloc_spinor_field_f(sfsize, &glat_even);
  sf2 = alloc_spinor_field_f(sfsize, &glat_even);
#else
  sf1 = alloc_spinor_field_f(sfsize, &glattice);
  sf2 = alloc_spinor_field_f(sfsize, &glattice);
#endif
  gaussian_spinor_field(sf1);
  gaussian_spinor_field(sf2);
  return_value_cpu = test_herm_cpu(&HH_cpu, "H", sf1, sf2, sfsize);
  return_value_gpu = test_herm_gpu(&HH_gpu, "H", sf1, sf2, sfsize);
*/
  return return_value_cpu + return_value_gpu;
}
