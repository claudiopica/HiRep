/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "spinor_field.h"
#include <string.h>

/* g1=g2 */
void suNg_field_copy(suNg_field *g1, suNg_field *g2)
{
#ifdef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(g1,g2);
#endif
  memcpy(g1->ptr,g2->ptr,4*g1->type->gsize*sizeof(*(g1->ptr)));
}

