/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include <string.h>
#include "spinor_field.h"
#include "update.h"

/* g1=g2 */
void suNg_field_copy(suNg_field *g1, suNg_field *g2)
{
#ifdef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(g1,g2);
#endif
  memcpy(g1->ptr,g2->ptr,4*g1->type->gsize_gauge*sizeof(*(g1->ptr)));
}

void suNf_field_copy(suNf_field *g1, suNf_field *g2)
{
#ifdef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(g1,g2);
#endif
  memcpy(g1->ptr,g2->ptr,4*g1->type->gsize_gauge*sizeof(*(g1->ptr)));
}

void suNg_scalar_field_copy(suNg_scalar_field *g1, suNg_scalar_field *g2)
{
#ifdef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(g1,g2);
#endif
  memcpy(g1->ptr,g2->ptr,g1->type->gsize_spinor*sizeof(*(g1->ptr)));
}

