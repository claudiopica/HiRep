/*******************************************************************************
*
* coherence of the dirac float with the dirac op
*
*******************************************************************************/

#include "libhr.h"

static double hmass = 0.1;

static void loc_D(spinor_field *out, spinor_field *in)
{
  Dphi(hmass, out, in);
}

static void loc_D_flt(spinor_field_flt *out, spinor_field_flt *in)
{
  Dphi_flt(hmass, out, in);
}

int main(int argc, char *argv[])
{
  double sig, tau;
  spinor_field *s0, *s1;
  spinor_field_flt *f0, *f1;
  int return_value = 0;

  /* setup process id and communications */
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();
  u_gauge_f_flt = alloc_gfield_f_flt(&glattice);

  s0 = alloc_spinor_field_f(2, &glattice);
  s1 = s0 + 1;
  f0 = alloc_spinor_field_f_flt(2, &glattice);
  f1 = f0 + 1;

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  fflush(stdout);

  random_u(u_gauge);

  //start_gf_sendrecv(u_gauge); // Is this necessary? (SAM)

  represent_gauge_field();
  lprintf("MAIN", 0, "done.\n");

  lprintf("MAIN", 0, "Generating a random spinor field... ");
  fflush(stdout);

  spinor_field_zero_f(s0);
  gaussian_spinor_field(s0);
  lprintf("MAIN", 0, "done.\n");

  tau = 1. / sqrt(spinor_field_sqnorm_f(s0));
  spinor_field_mul_f(s0, tau, s0);
  assign_sd2s(f0, s0);

  assign_ud2u_f();

  #ifdef WITH_GPU
    copy_to_gpu_spinor_field_f(s0);
    copy_to_gpu_spinor_field_f(s1);
    copy_to_gpu_spinor_field_f_flt(f0);
    copy_to_gpu_spinor_field_f_flt(f1);
    copy_to_gpu_gfield_f(u_gauge_f);
  #endif

  #ifdef WITH_MPI
    start_sendrecv_gpu_gfield_f(u_gauge_f);
    complete_sendrecv_gpu_gfield_f(u_gauge_f);

    start_sendrecv_gpu_gfield_f_flt(u_gauge_f_flt);
    complete_sendrecv_gpu_gfield_f_flt(u_gauge_f_flt);
  #endif

  loc_D(s1, s0);
  loc_D_flt(f1, f0);

  #ifdef WITH_GPU
    copy_from_gpu_spinor_field_f(s0);
    copy_from_gpu_spinor_field_f(s1);
    copy_from_gpu_spinor_field_f_flt(f0);
    copy_from_gpu_spinor_field_f_flt(f1);
  #endif

  lprintf("INFO", 0, "Spinor sqnorm result double precision: %0.2e\n", spinor_field_sqnorm_f_cpu(s1));
  lprintf("INFO", 0, "Spinor sqnorm result single precision: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(f1));

  assign_sd2s(f0, s1);

  spinor_field_mul_add_assign_f_flt(f0, -1.0, f1);
  sig = spinor_field_sqnorm_f_flt(f0);

  lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", sqrt(sig));
  lprintf("MAIN", 0, "(should be around 1*10^(-8) or so)\n\n");

  if (sqrt(sig) > 10.e-7)
    return_value = 1;

  free_spinor_field_f(s0);
  free_spinor_field_f_flt(f0);

  finalize_process();
  return return_value;
}
