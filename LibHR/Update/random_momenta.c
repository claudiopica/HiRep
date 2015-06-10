/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "spinor_field.h"
#include "suN_repr_func.h"
#include "random.h"
#include "update.h"
#include "utils.h"
#include <math.h>

void gaussian_momenta(suNg_av_field *momenta) {
  geometry_descriptor *gd=momenta->type;
  
  const double c3=1./sqrt(_FUND_NORM2);
  const int ngen=NG*NG-1;
  //const int ngen=NG*(NG-1)/2;
  
//Avoid OMP parallel region in PIECE_FOR
#undef _OMP_PRAGMA
#define _OMP_PRAGMA(s)
  
  _PIECE_FOR(gd,ixp) {
    int start=gd->master_start[ixp];
    int len=ngen*4*(gd->master_end[ixp]-start+1); /* lenght in doubles */
    double *dptr=(double*)(momenta->ptr+4*(start-gd->master_shift));
    gauss(dptr,len);
    for (int ix=0; ix<len; ++ix, ++dptr) {
      *(dptr)*=c3;
    }
  }
  
  apply_BCs_on_momentum_field(momenta);
}
