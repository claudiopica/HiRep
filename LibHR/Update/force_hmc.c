/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Martin Hansen                           *
* All rights reserved.                                                      *
\***************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "inverters.h"
#include "io.h"

spinor_field *Xs = NULL;
spinor_field *Ys = NULL;
spinor_field *eta = NULL;

#ifdef UPDATE_EO
static spinor_field *xi = NULL;
#endif

void free_force_hmc() {
    free_spinor_field(Xs);
#ifdef UPDATE_EO
    free_spinor_field(eta);
    free_spinor_field(xi);
#endif
}

void init_force_hmc() {
    static int init = 0;
    if (init == 0) {
#ifndef UPDATE_EO
        Xs = alloc_spinor_field(3, &glattice);
        Ys = Xs + 1;
        eta = Ys + 1;
#else
        Xs = alloc_spinor_field(2, &glattice);
        Ys = Xs + 1;
        eta = alloc_spinor_field(1, &glat_even);
        xi = alloc_spinor_field(1, &glat_odd);
        Xs->type = &glat_even;
        Ys->type = &glat_even;
#endif
        init = 1;
        atexit(free_force_hmc);
    }
}

void force_hmc(double dt, void *vpar) {
    // int n_iters = 0;
    force_hmc_par *par = (force_hmc_par *)vpar;
    suNg_av_field *force = *par->momenta;
    spinor_field *pf = par->pf;

    init_force_hmc();
    set_dirac_mass(par->mass);
    set_twisted_mass(par->mu);
    fermion_force_begin();

    double tmp;
    mshift_par mpar;
    mpar.err2 = par->inv_err2;
    mpar.max_iter = 0;
    mpar.n = 1;
    mpar.shift = &tmp;
    mpar.shift[0] = 0;

#ifndef UPDATE_EO

    if (par->mu == 0 || par->hasenbusch != 0) {
        /* X = H^{-1} pf = D^{-1} g5 pf */
        zero_spinor_field(Xs);
        g5_assign_spinor_field(pf);
        // n_iters +=
        g5QMR_mshift(&mpar, &D, pf, Xs);
        g5_assign_spinor_field(pf);

        if (par->hasenbusch == 0) {
            /* Y =  D^{-1} ( g5 X ) */
            g5_spinor_field(eta, Xs);
        } else if (par->hasenbusch == 1) {
            /* Y = H^{-1} ( g5 pf[k] + b X ) = D^{-1}( pf + b g5 X ) */
            g5_spinor_field(eta, Xs);
            mul_spinor_field(eta, par->b, eta);
            add_assign_spinor_field(eta, pf);
        } else if (par->hasenbusch == 2) {
            /* Y= -i D^{-1} ( pf[k] + imu g5 X )*/
            double mu1 = par->mu;
            double mu2 = par->mu + par->b;
            double muS = mu2 * mu2 - mu1 * mu1;
            g5_spinor_field(eta, Xs);
            mul_spinor_field(Xs, muS, Xs);
        }

        zero_spinor_field(Ys);
        // n_iters +=
        g5QMR_mshift(&mpar, &D, eta, Ys);
    } else {
        // n_iters +=
        cg_mshift(&mpar, QpQm_tm_alt, pf, Xs);
        Qtm_p_alt(Ys, Xs);
    }

#else
    if (par->mu == 0) {
        /* X_e = H^{-1} pf */
        /* X_o = D_{oe} X_e = D_{oe} H^{-1} pf */
        g5_assign_spinor_field(pf);
        mre_guess(&par->mpar, 0, Xs, &D, pf);
        // n_iters +=
        g5QMR_mshift(&mpar, &D, pf, Xs);
        mre_store(&par->mpar, 0, Xs);
        g5_assign_spinor_field(pf);

        /* Y_e = H^{-1} ( g5 pf + b X_e ) */
        /* Y_o = D_oe H^{-1} ( g5 pf + b X_e ) */
        if (par->hasenbusch != 1) {
            copy_spinor_field(eta, Xs);
        } else {
            g5_spinor_field(eta, pf);
            mul_add_assign_spinor_field(eta, par->b, Xs);
        }

        g5_assign_spinor_field(eta);
        mre_guess(&par->mpar, 1, Ys, &D, eta);
        // n_iters +=
        g5QMR_mshift(&mpar, &D, eta, Ys);
        mre_store(&par->mpar, 1, Ys);
        g5_assign_spinor_field(eta);

        if (par->hasenbusch == 2) {
            double muS = par->b * par->b;
            mul_spinor_field(Xs, muS, Xs);
        }
    } else {
        /* Ye = 1/(QpQm+mu^2) \phi */
        mre_guess(&par->mpar, 0, Ys, QpQm_tm_alt, pf);
        // n_iters += 2 *
        cg_mshift(&mpar, QpQm_tm_alt, pf, Ys);
        mre_store(&par->mpar, 0, Ys);
        Qtm_m_alt(Xs, Ys);

        if (par->hasenbusch == 2) {
            double mu1 = par->mu;
            double mu2 = par->mu + par->b;
            double muS = mu2 * mu2 - mu1 * mu1;
            mul_spinor_field(Xs, muS, Xs);
        }
    }
    // to here

#endif

    if (par->hasenbusch != 1) {
        force_fermion_core(Xs, Ys, 1, dt, 1.);
    } else {
        force_fermion_core(Xs, Ys, 1, dt, par->b);
    }

#ifdef WITH_CLOVER_EO

    if (par->logdet) {
        compute_force_logdet(par->mass, 2.); // 2 = # of flavors
    }

#endif

    fermion_force_end(dt, force);
}
