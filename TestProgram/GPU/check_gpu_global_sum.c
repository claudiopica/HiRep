/*******************************************************************************
*
*  NOCOMPILE= WITH_MPI
*
* Test of modules
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
#include "global.h"
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "setup.h"
#include "logger.h"

static double EPSILON=1.e-12;

int main(int argc, char *argv[])
{

      double *vector;
      int n=1024;
      int return_value = 0;
      /* setup process id and communications */
      logger_map("DEBUG", "debug");

      setup_process(&argc, &argv);

      cudaMallocManaged((void **) &vector,n*sizeof(double),cudaMemAttachGlobal);
      //lprintf("GPU TEST",0,"difference: %.2e\n",d);
      /*if (d > EPSILON){
            lprintf("GPU",0,"Test failed ?\n");
            return_value +=1;
      }*/
      
      

      finalize_process();
      return return_value;
}
