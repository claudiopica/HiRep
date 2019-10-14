/******************************************************************************
*
* Test of Schwarz Alternating Procedure
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
#include "communications.h"
#include "setup.h"

int nhb, nor, nit, nth, nms, level, seed;
double beta;
spinor_field *tmp;
static double hmass = 0.1;

void M(spinor_field *out, spinor_field *in)
{
  {
    empty_buffers(tmp);
    tmp->type = in->type;
    g5Dphi(-hmass, tmp, in);
    g5Dphi(-hmass, out, tmp);
  }
}
int main(int argc, char *argv[])
{

  setup_process(&argc, &argv);

  setup_gauge_fields();
 
#ifdef WITH_UNTESTED

  int i;
  double tau;
  spinor_field *s1, *s2;
  spinor_field *res;

  mshift_par par;

  int cgiters;
 lprintf("MAIN", 0, "Generating a random gauge field... ");
  fflush(stdout);
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);
  lprintf("MAIN", 0, "done.\n");
  represent_gauge_field();

  lprintf("REV TEST", 0, "Initial plaquette: %1.8e\n", avr_plaquette());

  par.n = 1;
  par.shift = (double *)malloc(sizeof(double) * (par.n));
  par.err2 = 1.e-10;
  par.max_iter = 0;
  res = alloc_spinor_field_f(par.n + 3, &glattice);
  s1 = res + par.n;
  s2 = s1 + 1;
  tmp = s2 + 1;

  par.shift[0] = 0.0;
  /*
   lprintf("res->type",0,"\ninner_master_pieces: %d\n",res->type->inner_master_pieces);
   lprintf("res->type",0,"local_master_pieces: %d\n",res->type->local_master_pieces);
   lprintf("res->type",0,"total_spinor_master_pieces: %d\n",res->type->total_spinor_master_pieces);
   lprintf("res->type",0,"total_gauge_master_pieces: %d\n",res->type->total_gauge_master_pieces);
   lprintf("res->type",0,"master_start[0]: %d\n",res->type->master_start[0]);
   lprintf("res->type",0,"master_start[1]: %d\n\n",res->type->master_start[1]);
   */

  gaussian_spinor_field(s1);

  g5Dphi(0.1, s2, s1);

  /* TEST CG_M */

  /*s1->type = &glat_red;
  s2->type = &glat_red;*/

  /*spinor_field_zero_f(s1);*/

  lprintf("CGTEST", 0, "spinor_field_sqnorm_f(s1)=%e\n", spinor_field_sqnorm_f(s1)); // ULRIK

  g5Dphi(0.1, s2, s1);

  lprintf("CGTEST", 0, "spinor_field_sqnorm_f(s1)=%e\n", spinor_field_sqnorm_f(s1)); // ULRIK
  lprintf("CGTEST", 0, "spinor_field_sqnorm_f(s2)=%e\n", spinor_field_sqnorm_f(s2)); // ULRIK

  s1->type = &glattice;
  s2->type = &glattice;

  lprintf("SAP TEST", 0, "Testing SAP\n");
  lprintf("SAP TEST", 0, "---------------------\n");
  spinor_field_zero_f(res);
  cgiters = cg_mshift(&par, &M, s1, res);
  lprintf("SAP TEST", 0, "CG_mshift converged in = %d steps\n", cgiters);
  cgiters = 0;
  spinor_field_zero_f(res);
  SAP_prec(5, &cg_mshift, &par, &M, s1, res);

  lprintf("SAP TEST", 0, "Converged in %d iterations\n", cgiters);
  for (i = 0; i < par.n; ++i)
  {
    M(s2, &res[i]);
    spinor_field_mul_add_assign_f(s2, -par.shift[i], &res[i]);
    spinor_field_sub_assign_f(s2, s1);
    tau = spinor_field_sqnorm_f(s2) / spinor_field_sqnorm_f(s1);
    lprintf("SAP TEST", 0, "test SAP = %e (req. %e)\n", tau, par.err2);
  }

  free_spinor_field_f(res);
  free(par.shift);

#else
 lprintf("MAIN",0,"Warning no test performed, the routine SAP_prec is in UNTESTED status");
#endif

  finalize_process();
  return 0;
}
