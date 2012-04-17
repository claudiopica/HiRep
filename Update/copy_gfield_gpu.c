/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/
#ifdef WITH_GPU
#include "spinor_field.h"
#include <string.h>

/* g1=g2 */
void suNg_field_copy(suNg_field *g1, suNg_field *g2)
{
  _TWO_SPINORS_MATCHING(g1,g2);
  
  cudaMemcpy(g1->gpu_ptr,g2->gpu_ptr,4*g1->type->gsize*sizeof(*(g1->ptr)),cudaMemcpyDeviceToDevice);
}
#endif

