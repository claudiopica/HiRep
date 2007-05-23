/*******************************************************************************
*
* File random_fields.c
*
* Pseudorandom generation of fields
*
*******************************************************************************/

#include <stdlib.h>

#include "random.h"
#include "error.h"
#include "global.h"

void random_u(void)
{
   int ix,mu;

   error(u_gauge==NULL,1,"random_u [random_fields.c]",
         "Attempt to access unallocated memory space");   

   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
         random_suNg(pu_gauge(ix,mu));
   }
}
