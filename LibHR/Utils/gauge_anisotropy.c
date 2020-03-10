/***************************************************************************\
* Copyright (c) 2008, Antonio Rago                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#include "global.h"
#include "logger.h"

void init_pure_gauge_anisotropy(double *chi)
{
#ifdef PURE_GAUGE_ANISOTROPY

    error(plaq_weight == NULL, 0, "init_pure_gauge_anisotropy", "In order to use anisotropic lattice you must compile with PLAQ_WEIGHTS enabled");
    error(*chi <= 0., 0, "init_pure_gauge_anisotropy", "The anisotropy factor must be positive");

    int ix, iy, iz, it;
    int mu, nu, index;

    for (ix = 0; ix < X_EXT; ++ix)
        for (iy = 0; iy < Y_EXT; ++iy)
            for (iz = 0; iz < Z_EXT; ++iz)
                for (it = 0; it < T_EXT; ++it)
                {
                    index = ipt_ext(it, ix, iy, iz);
                    if (index != -1)
                    {
                        mu = 0;
                        for (nu = mu + 1; nu < 4; nu++)
                        {
                            plaq_weight[index * 16 + mu * 4 + nu] *= *chi;
                            plaq_weight[index * 16 + nu * 4 + mu] *= *chi;
                        }
                        for (mu = 1; mu < 3; mu++)
                            for (nu = mu + 1; nu < 4; nu++)
                            {
                                plaq_weight[index * 16 + mu * 4 + nu] /= *chi;
                                plaq_weight[index * 16 + nu * 4 + mu] /= *chi;
                            }
                    }
                }
#else
    error(0 == 0, 0, "init_pure_gauge_anisotropy", "In order to use anisotropic lattice you must compile with PURE_GAUGE_ANISOTROPY enabled");
#endif
}