/*******************************************************************************
*
* File staples_gpu.hpp
*
* Computation of the "staples" (on GPU)
* port of staples.c
*
*******************************************************************************/

#include "update.h"
#include "libhr_core.h"

__device__ __forceinline__ void staples_dev(int ix, int mu, suNg *v, suNg *gauge, int *iup_gpu, int *idn_gpu,
                                            double *plaq_weight) {
    suNg staple, tr1, tr2;
    suNg p1, p2, p3;
    int nu, i, ixpmu, ixpnu, ixmnu, ixpmumnu;

    ixpmu = iup_gpu[4 * ix + mu];
    _suNg_zero(*v);

    for (i = 1; i < 4; i++) {
        nu = (mu + i) & 0x3;
        ixpnu = iup_gpu[4 * ix + nu];
        ixmnu = idn_gpu[4 * ix + nu];
        ixpmumnu = idn_gpu[4 * ixpmu + nu];

        //Up Staple
        read_gpu<double>(0, &p1, gauge, ix, nu, 4);
        read_gpu<double>(0, &p2, gauge, ixpnu, mu, 4);
        read_gpu<double>(0, &p3, gauge, ixpmu, nu, 4);

        _suNg_times_suNg(tr2, p1, p2);

        _suNg_dagger(tr1, p3);
        _suNg_times_suNg(staple, tr2, tr1);

#ifdef PLAQ_WEIGHTS
        if (plaq_weight != NULL) { _suNg_mul(staple, plaq_weight[ix * 16 + nu * 4 + mu], staple); }
#endif
        _suNg_add_assign(*v, staple);

        //Down Staple
        read_gpu<double>(0, &p1, gauge, ixmnu, mu, 4);
        read_gpu<double>(0, &p2, gauge, ixpmumnu, nu, 4);
        read_gpu<double>(0, &p3, gauge, ixmnu, nu, 4);

        _suNg_times_suNg(tr2, p1, p2);
        _suNg_dagger(tr1, p3);
        _suNg_times_suNg(staple, tr1, tr2);

#ifdef PLAQ_WEIGHTS
        if (plaq_weight != NULL) { _suNg_mul(staple, plaq_weight[ixmnu * 16 + mu * 4 + nu], staple); }
#endif
        _suNg_add_assign(*v, staple);
    }
}