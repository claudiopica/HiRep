/*******************************************************************************
* Check that the molecular dynamics evolution is reversible
*
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

int main(int argc, char *argv[])
{

  int return_value = 0;
  double orig_plaq, new_plaq, diff;
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);
  setup_gauge_fields();

  //  return 1;

  /* Init Monte Carlo */
  init_mc_ghmc(&flow, get_input_filename());

  lprintf("MAIN", 0, "MVM during (R)HMC initialzation: %ld\n", getMVM());
  lprintf("MAIN", 0, "Initial plaquette: %1.8e\n", avr_plaquette());

  /*Vincent */
  lprintf("REVERSIBILITY TEST", 0, "Plaquette before update: %1.8e\n", avr_plaquette());
  orig_plaq = avr_plaquette();
  int rr = update_ghmc();
  /*Vincent */
  if (rr < 0)
  {
    lprintf("REVERSIBILITY TEST", 0, "Error in updating the gauge field!!\n");
    return 1;
  }

  if (rr == 1)
  {

    lprintf("REVERSIBILITY TEST", 0, "Plaquette after update: %1.8e\n", avr_plaquette());

    rr = reverse_update_ghmc();
    /*Vincent */
    if (rr < 0)
    {
      lprintf("REVERSIBILITY TEST", 0, "Error in updating the gauge field!!\n");
      return 1;
    }

    lprintf("REVERSIBILITY TEST", 0, "Plaquette after reverse update: %1.8e\n", avr_plaquette());
    new_plaq = avr_plaquette();

    diff = fabs(new_plaq - orig_plaq);
    lprintf("REVERSIBILITY TEST", 0, "diff: %1.8e\n Should be 10^-12 or so\n", diff);
    if (diff > 1e-10)
      return_value++;
  }
  else
    lprintf("REVERSIBILITY TEST", 0, "Skipped the comparison as the configuration was not accepted\n");

  /* finalize Monte Carlo */
  end_mc();

  finalize_process();
  return return_value;
}
