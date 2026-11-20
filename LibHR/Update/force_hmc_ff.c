/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "inverters.h"
#include "io.h"
#include "geometry.h"

#define _print_avect(a)                                                                                                   \
    printf("(%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e)\n", (a).c1, (a).c2, (a).c3, (a).c4, (a).c5, (a).c6, (a).c7, \
           (a).c8)

#define _print_mat(a)                                                                                                \
    printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n", (a).c1_1.re, (a).c1_2.re, (a).c1_3.re, \
           (a).c2_1.re, (a).c2_2.re, (a).c2_3.re, (a).c3_1.re, (a).c3_2.re, (a).c3_3.re);                            \
    printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n", (a).c1_1.im, (a).c1_2.im, (a).c1_3.im, \
           (a).c2_1.im, (a).c2_2.im, (a).c2_3.im, (a).c3_1.im, (a).c3_2.im, (a).c3_3.im)

/* we need to compute  Tr  U(x,mu) g_5*(1-g_mu) chi2 # chi1^+
 * where # indicates the tensor product and Tr is the trace on Lorentz space.
 * the strategy is the following:
 * given the form of g_5(1-g_mu) one can compute only the first two lorentz
 * components of the spinor; so we first apply g_5(1-g_mu) to chi2 to find the first
 * two components; then we multiply these two vectors by U(x,mu) and
 * store the result in p.c[0], p.c[1]; when computing the trace we can factorize p.c[0] and p.c[1]
 * as they both multiply two components of chi1^+; we store these factors in p.c[2] and p.c[3].
 * the tensor product is performed by the macro 
 * _suNf_FMAT(u,p): u = p.c[0] # p.c[2]^+ + p.c[1] # p.c[3]^+
 */

/* these macros use the variables ptmp, p */

#ifdef BC_T_THETA
#define _T_theta_mulc(r)                   \
    _vector_mulc_f(ptmp, eitheta[0], (r)); \
    (r) = ptmp
#else
#define _T_theta_mulc(r)
#endif
#ifdef BC_X_THETA
#define _X_theta_mulc(r)                   \
    _vector_mulc_f(ptmp, eitheta[1], (r)); \
    (r) = ptmp
#else
#define _X_theta_mulc(r)
#endif
#ifdef BC_Y_THETA
#define _Y_theta_mulc(r)                   \
    _vector_mulc_f(ptmp, eitheta[2], (r)); \
    (r) = ptmp
#else
#define _Y_theta_mulc(r)
#endif
#ifdef BC_Z_THETA
#define _Z_theta_mulc(r)                   \
    _vector_mulc_f(ptmp, eitheta[3], (r)); \
    (r) = ptmp
#else
#define _Z_theta_mulc(r)
#endif

#define _F_DIR0(u, chi1, chi2)                         \
    _vector_add_f(ptmp, (chi2)->c[0], (chi2)->c[2]);   \
    _suNf_multiply(p.c[0], *(pu_gauge_f(x, 0)), ptmp); \
    _T_theta_mulc(p.c[0]);                             \
    _vector_add_f(ptmp, (chi2)->c[1], (chi2)->c[3]);   \
    _suNf_multiply(p.c[1], *(pu_gauge_f(x, 0)), ptmp); \
    _T_theta_mulc(p.c[1]);                             \
    _vector_sub_f(p.c[2], (chi1)->c[0], (chi1)->c[2]); \
    _vector_sub_f(p.c[3], (chi1)->c[1], (chi1)->c[3]); \
    _suNf_FMAT((u), p)

#define _F_DIR1(u, chi1, chi2)                           \
    _vector_i_add_f(ptmp, (chi2)->c[0], (chi2)->c[3]);   \
    _suNf_multiply(p.c[0], *(pu_gauge_f(x, 1)), ptmp);   \
    _X_theta_mulc(p.c[0]);                               \
    _vector_i_add_f(ptmp, (chi2)->c[1], (chi2)->c[2]);   \
    _suNf_multiply(p.c[1], *(pu_gauge_f(x, 1)), ptmp);   \
    _X_theta_mulc(p.c[1]);                               \
    _vector_i_sub_f(p.c[2], (chi1)->c[0], (chi1)->c[3]); \
    _vector_i_sub_f(p.c[3], (chi1)->c[1], (chi1)->c[2]); \
    _suNf_FMAT((u), p)

#define _F_DIR2(u, chi1, chi2)                         \
    _vector_add_f(ptmp, (chi2)->c[0], (chi2)->c[3]);   \
    _suNf_multiply(p.c[0], *(pu_gauge_f(x, 2)), ptmp); \
    _Y_theta_mulc(p.c[0]);                             \
    _vector_sub_f(ptmp, (chi2)->c[1], (chi2)->c[2]);   \
    _suNf_multiply(p.c[1], *(pu_gauge_f(x, 2)), ptmp); \
    _Y_theta_mulc(p.c[1]);                             \
    _vector_sub_f(p.c[2], (chi1)->c[0], (chi1)->c[3]); \
    _vector_add_f(p.c[3], (chi1)->c[1], (chi1)->c[2]); \
    _suNf_FMAT((u), p)

#define _F_DIR3(u, chi1, chi2)                           \
    _vector_i_add_f(ptmp, (chi2)->c[0], (chi2)->c[2]);   \
    _suNf_multiply(p.c[0], *(pu_gauge_f(x, 3)), ptmp);   \
    _Z_theta_mulc(p.c[0]);                               \
    _vector_i_sub_f(ptmp, (chi2)->c[1], (chi2)->c[3]);   \
    _suNf_multiply(p.c[1], *(pu_gauge_f(x, 3)), ptmp);   \
    _Z_theta_mulc(p.c[1]);                               \
    _vector_i_sub_f(p.c[2], (chi1)->c[0], (chi1)->c[2]); \
    _vector_i_add_f(p.c[3], (chi1)->c[1], (chi1)->c[3]); \
    _suNf_FMAT((u), p)

//Temporary storage fields from force_hmc.c
extern spinor_field *Xs, *Ys, *eta;
extern spinor_field_flt *eta_flt;

void force_hmc_ff(double dt, void *vpar) {
#ifdef UPDATE_EO
    spinor_field Xe, Xo, Ye, Yo;
#endif

    init_force_hmc();
    force_hmc_par *par = (force_hmc_par *)vpar;
    suNg_av_field *force = *par->momenta;
    int n_iters = 0;
    (void)n_iters;

    set_ff_dirac_mass(par->mass);

    /* check input types */
    _TWO_SPINORS_MATCHING(u_gauge, force);

    for (int k = 0; k < par->n_pf; ++k) {
#ifndef UPDATE_EO
        // Not implemented
#else
        double mass = par->mass;
        spinor_field *pf = par->pf;

        /*    g5QMR_fltacc_par mpar;
    mpar.err2 = par->inv_err2;
    mpar.max_iter = 0;
    mpar.err2_flt = par->inv_err2_flt;
    mpar.max_iter_flt = 0;*/
        double tmp;
        mshift_par mpar;
        mpar.err2 = par->inv_err2;
        mpar.max_iter = 0;
        mpar.n = 1;
        mpar.shift = &tmp; // mpar.n = 1
        mpar.shift[0] = 0;

        double rho = 4. + mass;
        spinor_field *tmp_spinor_field = alloc_spinor_field(1, &glattice);

        /* Y_e = (D^ D)^{-1} pf[k]  */
        Ye = *Ys;
        Ye.type = &glat_even;
        Yo = *Ys;
        Yo.type = &glat_odd;
        Yo.ptr += glat_odd.master_shift;

        /* X_e = (D^)^(-1) pf[k] = g5 D Y_e */
        Xe = *Xs;
        Xe.type = &glat_even;
        Xo = *Xs;
        Xo.type = &glat_odd;
        Xo.ptr += glat_odd.master_shift;

        set_ff_dirac_mass(par->mass);

        if (par->hasenbusch == 2) {
            /*  Y_e = (D^ D)^{-1} (D+b) pf[k] */
            set_ff_dirac_shift(par->b);
            Dff_dagger(&Xe, &pf[k]);
            set_ff_dirac_shift(0.);
        } else {
            copy_spinor_field(&Xe, &pf[k]);
        }

        lprintf("FORCE", 0, " id %d mass %e shift %e \n", par->id, par->mass, par->b);

        zero_spinor_field(&Ye);
        mre_guess(&par->mpar, 0, &Ye, &Dff_sq, &Xe);
        n_iters += cg_mshift(&mpar, &Dff_sq, &Xe, &Ye);
        mre_store(&par->mpar, 0, &Ye);

        Dff(&Xe, &Ye);

        if (par->hasenbusch == 2) {
            /*  X_e = D (D^ D)^{-1} (D+b) pf[k] - pf[k] */
            mul_add_assign_spinor_field(&Xe, -1, &pf[k]);
        }

        g5_assign_spinor_field(&Xe);

        /* Y_o = A_o^{-1} D_oe (D^ D)^{-1} pf[k]  */
        Dphi_(&Yo, &Ye);
        spinor_sigma_pi_rho_div_assign(&Yo, ff_sigma, ff_pi, rho, &Yo);

        /* X_o = (A_o)^^{-1} D_{oe} X_e = (A_o^)^{-1} D_{oe} g5 (D^)^{-1} pf[k] */
        Dphi_(&Xo, &Xe);
        spinor_sigma_pi_dagger_rho_div_assign(&Xo, ff_sigma, ff_pi, rho, &Xo);

        //Add the force of the fermion fields on the auxiliary fields,
        //from the derivative of the Dirac operator.
        force_hmc_auxfields_fermion(dt, vpar, ff_sigma_mom, ff_pi_mom, Xs, Ys, par->hasenbusch);

        free_spinor_field(tmp_spinor_field);

#endif

        if (gauge_field_active) {
            /* reset force stat counters */
            start_sendrecv_spinor_field(Xs);
            start_sendrecv_spinor_field(Ys);

#ifdef WITH_NEW_GEOMETRY
            complete_sendrecv_spinor_field(Xs);
            complete_sendrecv_spinor_field(Ys);
#endif

            _PIECE_FOR(&glattice, xp) {
                suNg_algebra_vector f;
                suNf_vector ptmp;
                suNf_spinor p;
                suNf_FMAT s1;

                if (xp == glattice.inner_master_pieces) {
                    _OMP_PRAGMA(master) {
                        complete_sendrecv_spinor_field(Xs);
                        complete_sendrecv_spinor_field(Ys);
                    }
                    _OMP_PRAGMA(barrier)
                }

                _SITE_FOR(&glattice, xp, x) {
                    for (int mu = 0; mu < 4; ++mu) {
                        int y;
                        suNf_spinor *chi1, *chi2;
                        _suNf_FMAT_zero(s1);
                        switch (mu) {
                        case 0:
                            y = iup(x, 0);
                            chi1 = _FIELD_AT(Xs, x);
                            chi2 = _FIELD_AT(Ys, y);
                            _F_DIR0(s1, chi1, chi2);
                            chi1 = _FIELD_AT(Ys, x);
                            chi2 = _FIELD_AT(Xs, y);
                            _F_DIR0(s1, chi1, chi2);
                            break;
                        case 1:
                            y = iup(x, 1);
                            chi1 = _FIELD_AT(Xs, x);
                            chi2 = _FIELD_AT(Ys, y);
                            _F_DIR1(s1, chi1, chi2);
                            chi1 = _FIELD_AT(Ys, x);
                            chi2 = _FIELD_AT(Xs, y);
                            _F_DIR1(s1, chi1, chi2);
                            break;
                        case 2:
                            y = iup(x, 2);
                            chi1 = _FIELD_AT(Xs, x);
                            chi2 = _FIELD_AT(Ys, y);
                            _F_DIR2(s1, chi1, chi2);
                            chi1 = _FIELD_AT(Ys, x);
                            chi2 = _FIELD_AT(Xs, y);
                            _F_DIR2(s1, chi1, chi2);
                            break;
                        default: /* DIR 3 */
                            y = iup(x, 3);
                            chi1 = _FIELD_AT(Xs, x);
                            chi2 = _FIELD_AT(Ys, y);
                            _F_DIR3(s1, chi1, chi2);
                            chi1 = _FIELD_AT(Ys, x);
                            chi2 = _FIELD_AT(Xs, y);
                            _F_DIR3(s1, chi1, chi2);
                        }

                        _algebra_project_FMAT(f, s1);

#ifdef UPDATE_EO
                        if (par->hasenbusch != 2) {
                            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force, x, mu), -dt * (_REPR_NORM2 / _FUND_NORM2), f);
                        } else {
                            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force, x, mu), -dt * (_REPR_NORM2 / _FUND_NORM2), f);
                        }
#else
                        if (par->hasenbusch != 2) {
                            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force, x, mu), dt * (_REPR_NORM2 / _FUND_NORM2), f);
                        } else {
                            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force, x, mu), dt * (_REPR_NORM2 / _FUND_NORM2), f);
                        }
#endif

                    } //directions for
                } //SITE_FOR
            } //PIECE FOR
        } //if
    }
}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3
