/*************************************************************************** \
 * Copyright (c) 2008, Claudio Pica, Agostino Patella                        *
 * All rights reserved.                                                      *
\***************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "Utils/single_double_utils.h"
#include "Utils/boundary_conditions.h"
#include <math.h>
#include "memory.h"
#include "geometry.h"

void represent_gauge_field() {
#ifdef WITH_SMEARING
    smear_gauge_field();
#endif

#ifdef ALLOCATE_REPR_GAUGE_FIELD
#ifndef WITH_GPU
    /* loop on local lattice first */
    /* loop on the rest of master sites */
    _OMP_PRAGMA(_omp_parallel)
    for (int ip = 0; ip < glattice.local_master_pieces; ip++) {
        _OMP_PRAGMA(_omp_for)
        for (int ix = glattice.master_start[ip]; ix <= glattice.master_end[ip]; ix++) {
            for (int mu = 0; mu < 4; mu++) {
#ifdef WITH_SMEARING
                suNg *u = _4FIELD_AT(u_gauge_s, ix, mu);
#else
                suNg *u = pu_gauge(ix, mu);
#endif //WITH_SMEARING
                suNf *Ru = pu_gauge_f(ix, mu);
#ifdef UNROLL_GROUP_REPRESENT
                _group_represent(*Ru, *u);
#else
                _group_represent2(Ru, u);
#endif
            }
        }
    }

    /* wait gauge field transfer */
    complete_sendrecv_gfield(u_gauge);

    /* loop on the rest of master sites */
    _OMP_PRAGMA(_omp_parallel)
    for (int ip = glattice.local_master_pieces; ip < glattice.total_gauge_master_pieces; ip++) {
        _OMP_PRAGMA(_omp_for)
        for (int ix = glattice.master_start[ip]; ix <= glattice.master_end[ip]; ix++) {
            for (int mu = 0; mu < 4; mu++) {
#ifdef WITH_SMEARING
                suNg *u = _4FIELD_AT(u_gauge_s, ix, mu);
#else
                suNg *u = pu_gauge(ix, mu);
#endif
                suNf *Ru = pu_gauge_f(ix, mu);
#ifdef UNROLL_GROUP_REPRESENT
                _group_represent(*Ru, *u);
#else
                _group_represent2(Ru, u);
#endif
            }
        }
    }
#else
    represent_gauge_field_gpu();
#endif

    apply_BCs_on_represented_gauge_field();
#else //ALLOCATE_REPR_GAUGE_FIELD
    static int first_time = 1;
    /* wait gauge field transfer */
    complete_sendrecv_gfield(u_gauge);

    if (first_time) {
        first_time = 0;
#ifdef WITH_SMEARING
        u_gauge_f = (suNf_field *)((void *)u_gauge_s);
#else
        u_gauge_f = (suNf_field *)((void *)u_gauge);
#endif
        //    apply_BCs_on_represented_gauge_field(); //Already applied when configuration read or initialized
    }
#endif //ALLOCATE_REPR_GAUGE_FIELD
#ifdef DPHI_FLT
    assign_ud2u_f();
#endif

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    compute_clover_term();
#endif
}
