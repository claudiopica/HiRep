/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "spinor_field.h"
#include "suN_repr_func.h"
#include "random.h"
#include <math.h>

void gaussian_momenta(suNg_av_field *momenta) {
  double *dptr;
  _DECLARE_INT_ITERATOR(ix);
  geometry_descriptor *gd=momenta->type;
  
  const double c3=1./sqrt(_FUND_NORM2);
  const int ngen=NG*NG-1;
  
  _PIECE_FOR(gd,ix) {
    int start=gd->master_start[_PIECE_INDEX(ix)];
    int len=ngen*4*(gd->master_end[_PIECE_INDEX(ix)]-start+1); /* lenght in doubles */
    dptr=(double*)(momenta->ptr+4*start);
    gauss(dptr,len);
    for (ix=0; ix<len; ++ix, ++dptr) {
      *(dptr)*=c3;
    }
  }
}
