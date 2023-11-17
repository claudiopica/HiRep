/***************************************************************************\
 * Copyright (c) 2017, Agostino Patella, Martin Hansen                     *
 * All rights reserved.                                                    *
\***************************************************************************/

/*

Main functions:

===  double lw_action(double beta, double c0, double c1);
  Returns the LW action, already summed over the global volume.

===  void lw_force(suNg_av_field *force, double beta, double c0, double c1);
  It sets force equal to the LW force.


Core functions:

===  void calculate_stfld(int comm);
  It calculates the staples in the local lattice and organizes them in the field stfld. If comm==true then it also
  retrieves the value of the staples in the communication buffers from neighbouring processes. Notice that not all
  communication buffers are actually filled but only the ones that are needed for the force calculation. IMPORTANT:
  This function assume that the communication buffers of the gauge field have been already filled.

===  double lw_action_density(int ix, double beta, double c0, double c1);
  It calculates the LW action density. IMPORTANT: This function assumes that the staples have been already calculated.


*/

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "io.h"
#include <string.h>
#include "utils.h"
#include "geometry.h"
#include "inverters.h"

#define COMM (1 == 1)
#define NOCOMM (1 == 0)

/*********************************************************************
mu = (nu+i+1)&0x3;
stfld[2*nu+0]->ptr[3*ix+i] = ix -> ix-nu -> ix-nu+mu -> ix+mu
stfld[2*nu+1]->ptr[3*ix+i] = ix -> ix+nu -> ix+nu+mu -> ix+mu
*********************************************************************/
static staple_field *stfld[8] = { NULL };

void calculate_stfld(int comm) {
    suNg wu1;

    if (stfld[0] == NULL) {
        for (int k = 0; k < 8; k++) {
            stfld[k] = alloc_staple_field(&glattice);
        }
    }

    for (int k = 0; k < 8; k++) {
        memset(stfld[k]->ptr, 0, glattice.gsize_gauge * sizeof(suNg) * 3);
    }

    _MASTER_FOR(&glattice, ix) {
        for (int nu = 0; nu < 4; nu++) {
            int ixpnu = iup(ix, nu);
            int ixmnu = idn(ix, nu);

            for (int i = 0; i < 3; i++) {
                int mu = (nu + i + 1) & 0x3;
                int ixpmu = iup(ix, mu);
                int ixpmumnu = idn(ixpmu, nu);

                // *---
                // |
                // *---

                _suNg_times_suNg(wu1, *pu_gauge(ixmnu, mu), *pu_gauge(ixpmumnu, nu));
                _suNg_dagger_times_suNg(*_3FIELD_AT(stfld[2 * nu + 0], ix, i), *pu_gauge(ixmnu, nu), wu1);

                // ---*
                //    |
                // ---*
                _suNg_times_suNg(wu1, *pu_gauge(ix, nu), *pu_gauge(ixpnu, mu));
                _suNg_times_suNg_dagger(*_3FIELD_AT(stfld[2 * nu + 1], ix, i), wu1, *pu_gauge(ixpmu, nu));
            }
        }
    }

    if (comm) {
        for (int k = 0; k < 8; k++) {
            start_sendrecv_staple_field(stfld[k]);
            complete_sendrecv_staple_field(stfld[k]);
        }
    }
}

double lw_action_density(int ix, double beta, double c0, double c1) {
    double plaqs = 0;
    double rects = 0;
    double p;
    suNg w1;

    for (int nu = 0; nu < 3; nu++) {
        for (int mu = nu + 1; mu < 4; mu++) {
            int i = mu - nu - 1;
            _suNg_times_suNg_dagger(w1, *_3FIELD_AT(stfld[2 * nu + 1], ix, i), *pu_gauge(ix, mu));
            _suNg_trace_re(p, w1);
#ifdef PLAQ_WEIGHTS
            p *= plaq_weight[16 * ix + 4 * mu + nu];
#endif
            plaqs -= p;
        }
    }

    for (int nu = 0; nu < 4; nu++) {
        for (int i = 0; i < 3; i++) {
            int ixpnu = ix;
            _suNg_times_suNg_dagger(w1, *_3FIELD_AT(stfld[2 * nu + 1], ixpnu, i), *_3FIELD_AT(stfld[2 * nu + 0], ixpnu, i));
            _suNg_trace_re(p, w1);
#ifdef PLAQ_WEIGHTS
            int mu = (nu + i + 1) & 0x3;
            ixpnu = idn(ix, nu);
            p *= rect_weight[16 * ixpnu + 4 * mu + nu];
#endif
            rects -= p;
        }
    }

    return (beta / NG) * (c0 * plaqs + c1 * rects);
}

double lw_action(double beta, double c0, double c1) {
    double s = 0;
    calculate_stfld(NOCOMM);

    _MASTER_FOR(&glattice, ix) {
        s += lw_action_density(ix, beta, c0, c1);
    }

    global_sum(&s, 1);
    return s;
}

void lw_local_action(scalar_field *loc_action, double beta, double c0, double c1) {
    calculate_stfld(COMM);
    int iy = ipt(2, 0, 0, 0);
    _MASTER_FOR(&glattice, ix) {
        *_FIELD_AT(loc_action, iy) += lw_action_density(ix, beta, c0, c1);
    }
}

void lw_force(double dt, void *vpar) {
    suNg ws[4], wu1, wu2;
    suNg_algebra_vector wf1;

    force_gauge_par *par = (force_gauge_par *)vpar;
    suNg_av_field *force = *par->momenta;
    double beta = par->beta;
    double c0 = par->c0;
    double c1 = par->c1;

    calculate_stfld(COMM);

    // Calculation of the force in (ix,mu).
    // In the drawings below, mu is the direction of the missing link.
    // The index nu labels the directions orthogonal to mu.
    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            _suNg_zero(ws[mu]);
        }

        for (int nu = 0; nu < 4; nu++) {
            for (int i = 0; i < 3; i++) {
                int mu = (nu + i + 1) & 0x3;

                // *---
                // |
                // *---
#ifdef PLAQ_WEIGHTS
                int ixmnu = idn(ix, nu);
                _suNg_mul(wu1, plaq_weight[16 * ixmnu + 4 * mu + nu], *_3FIELD_AT(stfld[2 * nu + 0], ix, i));
                _suNg_add_assign(ws[mu], wu1);
#else
                _suNg_add_assign(ws[mu], *_3FIELD_AT(stfld[2 * nu + 0], ix, i));
#endif

                // ---*
                //    |
                // ---*
#ifdef PLAQ_WEIGHTS
                _suNg_mul(wu1, plaq_weight[16 * ix + 4 * mu + nu], *_3FIELD_AT(stfld[2 * nu + 1], ix, i));
                _suNg_add_assign(ws[mu], wu1);
#else
                _suNg_add_assign(ws[mu], *_3FIELD_AT(stfld[2 * nu + 1], ix, i));
#endif
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            _suNg_times_suNg_dagger(wu1, *pu_gauge(ix, mu), ws[mu]);
            _fund_algebra_project(wf1, wu1);
            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force, ix, mu), dt * (-beta * c0 / NG), wf1);
            _suNg_zero(ws[mu]);
        }

        for (int nu = 0; nu < 4; nu++) {
            int ixpnu = iup(ix, nu);
            int ixmnu = idn(ix, nu);

            for (int i = 0; i < 3; i++) {
                int mu = (nu + i + 1) & 0x3;
                int ixpmu = iup(ix, mu);
                int ixpmunnu = idn(ixpmu, nu);

                // *---*---
                // |
                // *---*---
                _suNg_dagger_times_suNg(wu1, *pu_gauge(ixmnu, nu), *_3FIELD_AT(stfld[2 * nu + 0], ixmnu, i));
                _suNg_times_suNg(wu2, wu1, *pu_gauge(ixpmunnu, nu));
#ifdef PLAQ_WEIGHTS
                int ixmnumnu = idn(ixmnu, nu);
                _suNg_mul(wu2, rect_weight[16 * ixmnumnu + 4 * mu + nu], wu2);
#endif
                _suNg_add_assign(ws[mu], wu2);

                // ---*---*
                //        |
                // ---*---*
                _suNg_times_suNg(wu1, *pu_gauge(ix, nu), *_3FIELD_AT(stfld[2 * nu + 1], ixpnu, i));
                _suNg_times_suNg_dagger(wu2, wu1, *pu_gauge(ixpmu, nu));
#ifdef PLAQ_WEIGHTS
                _suNg_mul(wu2, rect_weight[16 * ix + 4 * mu + nu], wu2);
#endif
                _suNg_add_assign(ws[mu], wu2);
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            int ixpmu = iup(ix, mu);
#ifdef PLAQ_WEIGHTS
            int ixmmu = idn(ix, mu);
#endif
            for (int i = 0; i < 3; i++) {
                int nu = (mu + i + 1) & 0x3;
                int ixpnu = iup(ix, nu);
                int ixmnu = idn(ix, nu);
                int ixmnupmu = iup(ixmnu, mu);
#ifdef PLAQ_WEIGHTS
                int ixmnummu = idn(ixmnu, mu);
#endif
                // *---*
                // |   |
                // *   *
                // |
                // *---*
                _suNg_dagger_times_suNg(wu1, *pu_gauge(ixmnu, nu), *pu_gauge(ixmnu, mu));
                _suNg_times_suNg(wu2, wu1, *_3FIELD_AT(stfld[2 * mu + 1], ixmnupmu, i));
#ifdef PLAQ_WEIGHTS
                _suNg_mul(wu2, rect_weight[16 * ixmnu + 4 * nu + mu], wu2);
#endif
                _suNg_add_assign(ws[mu], wu2);

                // *---*
                // |   |
                // *   *
                //     |
                // *---*
                _suNg_times_suNg(wu1, *pu_gauge(ix, nu), *pu_gauge(ixpnu, mu));
                _suNg_times_suNg_dagger(wu2, wu1, *_3FIELD_AT(stfld[2 * mu + 1], ixpmu, i));
#ifdef PLAQ_WEIGHTS
                _suNg_mul(wu2, rect_weight[16 * ix + 4 * nu + mu], wu2);
#endif
                _suNg_add_assign(ws[mu], wu2);

                // *---*
                // |
                // *   *
                // |   |
                // *---*
                _suNg_dagger_times_suNg(wu1, *_3FIELD_AT(stfld[2 * mu + 0], ixmnu, i), *pu_gauge(ixmnu, mu));
                _suNg_times_suNg(wu2, wu1, *pu_gauge(ixmnupmu, nu));
#ifdef PLAQ_WEIGHTS
                _suNg_mul(wu2, rect_weight[16 * ixmnummu + 4 * nu + mu], wu2);
#endif
                _suNg_add_assign(ws[mu], wu2);

                // *---*
                //     |
                // *   *
                // |   |
                // *---*
                _suNg_times_suNg(wu1, *_3FIELD_AT(stfld[2 * mu + 0], ix, i), *pu_gauge(ixpnu, mu));
                _suNg_times_suNg_dagger(wu2, wu1, *pu_gauge(ixpmu, nu));
#ifdef PLAQ_WEIGHTS
                _suNg_mul(wu2, rect_weight[16 * ixmmu + 4 * nu + mu], wu2);
#endif
                _suNg_add_assign(ws[mu], wu2);
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            _suNg_times_suNg_dagger(wu1, *pu_gauge(ix, mu), ws[mu]);
            _fund_algebra_project(wf1, wu1);
            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force, ix, mu), dt * (-beta * c1 / NG), wf1);
        }
    }

    apply_BCs_on_momentum_field(force);
}
