/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "random.h"
#include "memory.h"
#include "libhr_core.h"
#include "Utils/boundary_conditions.h"

//Avoid OMP parallel region in PIECE_FOR
#undef _OMP_PRAGMA
#define _OMP_PRAGMA(s)

void gaussian_momenta(suNg_av_field *momenta) {
    geometry_descriptor *gd = momenta->type;

    const double c3 = 1. / sqrt(_FUND_NORM2);
    const int ngen = NG * NG - 1;
    //const int ngen=NG*(NG-1)/2; //TODO: check SON

    _PIECE_FOR(gd, ixp) {
        int start = gd->master_start[ixp];
        int len = ngen * 4 * (gd->master_end[ixp] - start + 1); /* length in doubles */
        double *dptr = (double *)(momenta->ptr + 4 * (start - gd->master_shift));
        gauss(dptr, len);
        for (int ix = 0; ix < len; ++ix, ++dptr) {
            *(dptr) *= c3;
        }
    }

#ifdef WITH_GPU
    copy_to_gpu_suNg_av_field(momenta);
#endif

    apply_BCs_on_momentum_field(momenta);
    start_sendrecv_suNg_av_field(momenta);
    complete_sendrecv_suNg_av_field(momenta);
}

void gaussian_scalar_momenta(suNg_scalar_field *momenta) {
    double c2 = 1. / sqrt(2.);

    _MASTER_FOR(momenta->type, ix) {
        suNg_vector *dptr = _FIELD_AT(momenta, ix);
        gauss((double *)dptr, sizeof(suNg_vector) / sizeof(double));
        _vector_mul_g(*dptr, c2, *dptr);
    }

#ifdef WITH_GPU
    copy_to_gpu_suNg_scalar_field(momenta);
#endif

    //apply_BCs_on_scalar_momentum_field(momenta); // TODO: this does not seem to exist. It probably should.
    start_sendrecv_suNg_scalar_field(momenta);
    complete_sendrecv_suNg_scalar_field(momenta);
}
