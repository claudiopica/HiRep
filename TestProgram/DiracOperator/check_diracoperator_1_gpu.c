/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Gauge covariance of the Dirac operator on GPU
*
*******************************************************************************/

// TODO: This test fails for unknown reasons. (SAM)

#include "libhr.h"

static double hmass = 0.1;
static suNg_field *g;

static void loc_D(spinor_field *out, spinor_field *in)
{
  Dphi(hmass, out, in);
}

static void random_g(void)
{
  _MASTER_FOR(&glattice, ix)
  {
    random_suNg(_FIELD_AT(g, ix));
  }
}

static void transform_u(void)
{
  _MASTER_FOR(&glattice, ix)
  {
    suNg v;
    for (int mu = 0; mu < 4; mu++)
    {
      int iy = iup(ix, mu);
      suNg *u = pu_gauge(ix, mu);
      _suNg_times_suNg_dagger(v, *u, *_FIELD_AT(g, iy));
      _suNg_times_suNg(*u, *_FIELD_AT(g, ix), v);
    }
  }

  start_sendrecv_gfield(u_gauge);
  complete_sendrecv_gfield(u_gauge);
  represent_gauge_field();
  smear_gauge_field();

  copy_to_gpu_gfield_f(u_gauge_f);
  start_sendrecv_gpu_gfield_f(u_gauge_f);
  complete_sendrecv_gpu_gfield_f(u_gauge_f);
}

static void transform_s(spinor_field *out, spinor_field *in)
{
  _MASTER_FOR(&glattice, ix)
  {
    suNf_spinor *s = _FIELD_AT(in, ix);
    suNf_spinor *r = _FIELD_AT(out, ix);
    suNf gfx;

    _group_represent2(&gfx, _FIELD_AT(g, ix));

    _suNf_multiply(r->c[0], gfx, s->c[0]);
    _suNf_multiply(r->c[1], gfx, s->c[1]);
    _suNf_multiply(r->c[2], gfx, s->c[2]);
    _suNf_multiply(r->c[3], gfx, s->c[3]);
  }
}

int main(int argc, char *argv[])
{

  int return_value = 0;
  double sig, tau;
  spinor_field *s0, *s1, *s2, *s3;

  // Setup process
  logger_map("DEBUG", "debug");
  setup_process(&argc, &argv);

  // Setup gauge fields
  setup_gauge_fields();
  g = alloc_gtransf(&glattice); /* allocate additional memory */
  random_u(u_gauge);
  represent_gauge_field();
  copy_to_gpu_gfield_f(u_gauge_f);
  start_sendrecv_gpu_gfield_f(u_gauge_f);
  complete_sendrecv_gpu_gfield_f(u_gauge_f);

  // Generate random gauge transformation to apply
  random_g();
  copy_to_gpu_gtransf(g);
  start_sendrecv_gpu_gtransf(g);
  complete_sendrecv_gpu_gtransf(g);

  // Setup initial gauge fields
  s0 = alloc_spinor_field_f(1, &glattice);
  s1 = alloc_spinor_field_f(1, &glattice);
  s2 = alloc_spinor_field_f(1, &glattice);
  s3 = alloc_spinor_field_f(1, &glattice);
  gaussian_spinor_field(s0);
  
  // Normalize s0 + Sanity check
  tau = 1. / sqrt(spinor_field_sqnorm_f_cpu(s0));
  spinor_field_mul_f_cpu(s0, tau, s0);
  sig = spinor_field_sqnorm_f_cpu(s0);
  lprintf("MAIN", 0, "Normalized norm = %.2e\n", sqrt(sig));

  // Apply Gauge TF
  lprintf("MAIN", 0, "Gauge covariance of the Dirac operator:\n");

  // Apply dirac operator on GPU + copy back and forth
  copy_to_gpu_spinor_field_f(s0);
  loc_D(s1, s0);
  copy_from_gpu_spinor_field_f(s1);
  copy_from_gpu_spinor_field_f(s0);
  spinor_field_zero_f(s2);
  spinor_field_zero_f(s3);

  // Gauge transformation on CPU
  transform_s(s2, s1);
  transform_s(s3, s0);
  transform_u();

  // Copy results to GPU, apply Dirac operator again
  spinor_field_zero_f(s1);
  copy_to_gpu_spinor_field_f(s2);
  copy_to_gpu_spinor_field_f(s3);
  loc_D(s1, s3);
  spinor_field_mul_add_assign_f(s1, -1.0, s2);
  sig = spinor_field_sqnorm_f(s1);

  // Print test results
  lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", sqrt(sig));
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");

  return_value += check_finiteness(sqrt(sig));
  return_value += check_diff_norm(sqrt(sig), 1e-14);

  // Free and return
  free_spinor_field_f(s0);
  free_spinor_field_f(s1);
  free_spinor_field_f(s2);
  free_spinor_field_f(s3);
  free_gtransf(g);
  finalize_process();
  return return_value;
}
