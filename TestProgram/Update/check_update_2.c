/*******************************************************************************
* NOCOMPILE= !WITH_EXPCLOVER
* NOCOMPILE= !NG=2 !NG=3
* NOCOMPILE= !REPR_FUNDAMENTAL
* Check the horner scheme for the implementation of the exp csw.
* This test is performed only for fundamental SU(2) and SU(3) fermions
* By: Fernando Romero-López 
* 
*******************************************************************************/

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
#include "../HMC/hmc_utils.h"
#include "setup.h"

int nhb, nor, nit, nth, nms, level, seed;
double beta;

//extern suNg_av_field *momenta;

hmc_flow flow = init_hmc_flow(flow);

/*
### this is essentially a copy of the function update_ghmc which takes in input the flipped momenta.
### this is to allow for a reversibility test.
### this might have to be changed  if update_ghmc is modified.
*/

void clover_field_copy(suNf_field *g1, suNf_field *g2)
{
#ifdef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(g1, g2);
#endif
  memcpy(g1->ptr, g2->ptr, 6 * g1->type->gsize_gauge * sizeof(*(g1->ptr)));
}

int main(int argc, char *argv[])
{

  int return_value = 0;
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);
  setup_gauge_fields();

  /* Init Monte Carlo */
  init_mc_ghmc(&flow, get_input_filename());

#ifdef WITH_EXPCLOVER
  // This is the check of the exponential clover term

  spinor_field *Xl = alloc_spinor_field_f(1, &glattice);
  spinor_field *Yl = alloc_spinor_field_f(1, &glattice);

  create_z2_volume_source(Xl);
  create_z2_volume_source(Yl);

  fermion_force_begin();
  force_clover_fermion_taylor(Xl, Yl, 1.0);

  suNf_field *cl_force2;
  cl_force2 = alloc_clover_force(&glattice);
  clover_field_copy(cl_force2, cl_force);

  fermion_force_begin();
  force_clover_fermion(Xl, Yl, 1.0);

  complex double aux;
  double norm = 0.;

  _MASTER_FOR(&glattice, ix)
  {
    for (int mu = 0; mu < 6; mu++)
    {
      for (int i = 0; i < NF * NF; i++)
      {

        aux = (_6FIELD_AT(cl_force2, ix, mu))->c[i] - (_6FIELD_AT(cl_force, ix, mu))->c[i];

        norm += (creal(aux) * creal(aux) + cimag(aux) * cimag(aux));
      }
    }
  }

  lprintf("TEST CLOVER FORCE", 0, "Difference = %2.20e \n", sqrt(norm / GLB_T / GLB_X / GLB_X / GLB_X / 6 / NF / NF / 2));

  if (sqrt(norm / (VOLUME*6.0* NF * NF)) > 1e-14)
  {
    lprintf("TEST CLOVER FORCE", 0, "FAILED!!! Precision is not enough.  \n");
    return_value += 1;
  }

#endif

  finalize_process();
  return return_value;
}
