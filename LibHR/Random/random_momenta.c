/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

// #include "spinor_field.h"
// #include "suN_repr_func.h"
// #include "update.h"
// #include "utils.h"
// #include <math.h>
// #include <stdio.h>

#include "random.h"
#include "libhr_core.h"
#include "Utils/boundary_conditions.h"

//Avoid OMP parallel region in PIECE_FOR
#undef _OMP_PRAGMA
#define _OMP_PRAGMA(s)

void gaussian_momenta(suNg_av_field *momenta)
{
    geometry_descriptor *gd = momenta->type;

    const double c3 = 1. / sqrt(_FUND_NORM2);
    const int ngen = NG * NG - 1;
    //const int ngen=NG*(NG-1)/2; //TODO: check SON

    _PIECE_FOR(gd, ixp) {
        int start = gd->master_start[ixp];
        int len = ngen * 4 * (gd->master_end[ixp] - start + 1); /* length in doubles */
        double *dptr = (double *)(momenta->ptr + 4 * (start - gd->master_shift));
        gauss(dptr, len);
        for (int ix = 0; ix < len; ++ix, ++dptr) { *(dptr) *= c3; }
    }

    apply_BCs_on_momentum_field(momenta);
}

void gaussian_scalar_momenta(suNg_scalar_field *momenta)
{
    double c2 = 1. / sqrt(2.);
#ifdef WITH_FUSE_MASTER_FOR
    _FUSE_MASTER_FOR(momenta->type, ix) {
        _FUSE_IDX(momenta->type, ix);
#else
    _MASTER_FOR(momenta->type, ix) {
#endif
        suNg_vector *dptr = _FIELD_AT(momenta, ix);
        gauss((double *)dptr, sizeof(suNg_vector) / sizeof(double));
        _vector_mul_g(*dptr, c2, *dptr);
    }

    //apply_BCs_on_scalar_momentum_field(momenta);
}
