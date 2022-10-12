/******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Test of hermiticity on GPU
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

void MM_gpu(spinor_field *out, spinor_field *in)
{
  #ifdef UPDATE_EO
    g5Dphi_eopre_sq(-hmass, out, in);
  #else
    g5Dphi_sq(-hmass, out, in);
  #endif
}

void II_cpu(spinor_field *out, spinor_field *in) 
{
  spinor_field_mul_f_cpu(out, 1, in);
}

void II_gpu(spinor_field *out, spinor_field *in) 
{
  spinor_field_mul_f(out, 1, in);
}

int test_herm_cpu(spinor_operator S, char *name)
{
  lprintf("RESULT", 0, "Test if %s is hermitian on CPU: \n", name);

  spinor_field *s1, *s2, *s3, *s4;
  double tau;
  int return_val = 0;

  // Prepare initial spinor fields
  #ifdef UPDATE_EO
    s1 = alloc_spinor_field_f(4, &glat_even);
  #else
    s1 = alloc_spinor_field_f(4, &glattice);
  #endif
  s2 = s1 + 1;
  s3 = s2 + 1;
  s4 = s3 + 1;
  gaussian_spinor_field(s1);
  gaussian_spinor_field(s2);

  // Apply operator
  S(s3, s1);
  S(s4, s2);

  // Spinor field sanity checks
  lprintf("RESULT", 0, "s1 NORM %f on CPU\n", sqrt(spinor_field_sqnorm_f_cpu(s1)));
  lprintf("RESULT", 0, "s2 NORM %f on CPU\n", sqrt(spinor_field_sqnorm_f_cpu(s2)));
  lprintf("RESULT", 0, "s3 NORM %f on CPU\n", sqrt(spinor_field_sqnorm_f_cpu(s3)));
  lprintf("RESULT", 0, "s4 NORM %f on CPU\n", sqrt(spinor_field_sqnorm_f_cpu(s4)));

  // Difference tau is 0 for a hermitian operator
  tau = spinor_field_prod_re_f_cpu(s2, s3);
  tau -= spinor_field_prod_re_f_cpu(s4, s1);
  tau += spinor_field_prod_im_f_cpu(s2, s3);
  tau -= spinor_field_prod_im_f_cpu(s4, s1);
  tau /= sqrt(spinor_field_sqnorm_f_cpu(s1));
  tau /= sqrt(spinor_field_sqnorm_f_cpu(s2));

  // Print test result info
  if (fabs(tau) > 1.e-14)
  {
    lprintf("RESULT", 0, "FAILED \n");
    return_val = 1;
  }
  else
    lprintf("RESULT", 0, "OK \n");
  lprintf("RESULT", 0, "[norm = %e]\n", tau);

  // Free and return
  free_spinor_field_f(s1);
  return return_val;
}

int test_herm_gpu(spinor_operator S, char *name)
{
  lprintf("RESULT", 0, "Test if %s is hermitian on GPU: \n", name);

  spinor_field *s1, *s2, *s3, *s4;
  double tau;
  int return_val = 0;

  // Prepare inital spinor fields
  #ifdef UPDATE_EO
    s1 = alloc_spinor_field_f(4, &glat_even);
  #else
    s1 = alloc_spinor_field_f(4, &glattice);
  #endif
  s2 = s1 + 1;
  s3 = s2 + 1;
  s4 = s3 + 1;
  gaussian_spinor_field(s1);
  gaussian_spinor_field(s2);
  spinor_field_copy_to_gpu_f(s1);
  spinor_field_copy_to_gpu_f(s2);

  // Apply operator
  S(s3, s1);
  S(s4, s2);

  // Spinor field sanity checks
  lprintf("RESULT", 0, "s1 NORM %f on GPU\n", sqrt(spinor_field_sqnorm_f(s1)));
  lprintf("RESULT", 0, "s2 NORM %f on GPU\n", sqrt(spinor_field_sqnorm_f(s2)));
  lprintf("RESULT", 0, "s3 NORM %f on GPU\n", sqrt(spinor_field_sqnorm_f(s3)));
  lprintf("RESULT", 0, "s4 NORM %f on GPU\n", sqrt(spinor_field_sqnorm_f(s4)));

  // Difference tau is 0 for a hermitian operator
  tau = spinor_field_prod_re_f(s2, s3);
  tau -= spinor_field_prod_re_f(s4, s1);
  tau += spinor_field_prod_im_f(s2, s3);
  tau -= spinor_field_prod_im_f(s4, s1);
  tau /= sqrt(spinor_field_sqnorm_f(s1));
  tau /= sqrt(spinor_field_sqnorm_f(s2));

  // Print test result info
  if (fabs(tau) > 1.e-14)
  {
    lprintf("RESULT", 0, "FAILED \n");
    return_val = 1;
  }
  else
    lprintf("RESULT", 0, "OK \n");
  lprintf("RESULT", 0, "[norm = %e]\n", tau);

  // Free and return
  free_spinor_field_f(s1);
  return return_val;
}

int main(int argc, char *argv[])
{
  int return_value_cpu, return_value_gpu;
  int return_value_cpu_unit, return_value_gpu_unit;

  // setup process id and communications
  logger_map("DEBUG", "debug");
  setup_process(&argc, &argv);

  // Setup gauge field
  setup_gauge_fields();
  lprintf("MAIN", 0, "Generating a random gauge field... ");
  random_u(u_gauge);
  lprintf("MAIN", 0, "done.\n");
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();

  // Test block
  
    // Q^2
    return_value_cpu=test_herm_cpu(&MM_cpu, "M");
    return_value_gpu=test_herm_gpu(&MM_gpu, "M");

    // Unit operator
    return_value_cpu_unit=test_herm_cpu(&II_cpu, "I");
    return_value_gpu_unit=test_herm_gpu(&II_gpu, "I");

  // Finalize
  finalize_process();
  return return_value_cpu + return_value_gpu;
}
