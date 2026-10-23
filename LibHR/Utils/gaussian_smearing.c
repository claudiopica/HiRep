/******************************************************************************
* File: gaussian_smearing.c
*
* Gaussian smearing for fermions
*
* Author: Pietro Butti
* Started: Nov. 2025
******************************************************************************/

#include "utils.h"
#include "libhr_core.h"
#include "memory.h"
#include "IO/logger.h"
#include "inverters.h"
#include "update.h"


/*******************************************\
    out <= F * in = (1 + aH) * in

    H(n,m)psi(m) = sum_i [U_i(n) psi(n+i)   +   U_i(n-i)^* psi(n-i)]
\*******************************************/
void gaussian_smearing(spinor_field *restrict out, spinor_field *restrict in, suNf_field *gauge_f, double alpha) {
#ifdef CHECK_SPINOR_MATCHING
    error((in == NULL) || (out == NULL), 1, "gaussian_smearing [gaussian_smearing.c]", "Attempt to access unallocated memory space");
    error(in == out, 1, "gaussian_smearing [gaussian_smearing.c]", "Input and output fields must be different");
    error(out->type == &glat_even && in->type == &glat_even, 1, "gaussian_smearing [gaussian_smearing.c]", "Spinors don't match! (1)");
    error(out->type == &glat_odd && in->type == &glat_odd, 1, "gaussian_smearing [gaussian_smearing.c]", "Spinors don't match! (2)");
#endif

    /************************ loop over all lattice sites *************************/
    /* start communication of input spinor field */
    start_sendrecv_spinor_field(in);

    _PIECE_FOR(out->type, ixp) {
#ifdef WITH_MPI
        if (ixp == out->type->inner_master_pieces) {
            /* wait for spinor to be transfered */
            complete_sendrecv_spinor_field(in);
        }
#endif

        _SITE_FOR(out->type, ixp, ix) {
#if defined(_OPENMP) && defined(WITH_PROBE_MPI)
            register int thread0 = hr_threadId();
#endif

            int i_dw, i_up;
            double factor = 1. / (1. + 6. * alpha);
            suNf *u, *u_dw;
            suNf_spinor *psi_up, *psi_dw, *r, *rold;
            suNf_vector chi0, chi1, chi2, chi3;

            r = _FIELD_AT(out, ix);
            rold = _FIELD_AT(in, ix);

            _spinor_zero_f(*r);

            for (int idir = 1; idir < 4; ++idir) {
                i_dw = idn(ix, idir);
                i_up = iup(ix, idir);

                u = (gauge_f->ptr) + coord_to_index(ix, idir);
                u_dw = (gauge_f->ptr) + coord_to_index(i_dw, idir);

                psi_up = _FIELD_AT(in, i_up);
                psi_dw = _FIELD_AT(in, i_dw);

                // up ------------------------------------------
                // chi_i <= U_ij(n) * psi_j(n+dir)
                _suNf_multiply(chi0, *(u), (*psi_up).c[0]);
                _suNf_multiply(chi1, *(u), (*psi_up).c[1]);
                _suNf_multiply(chi2, *(u), (*psi_up).c[2]);
                _suNf_multiply(chi3, *(u), (*psi_up).c[3]);

                // r <= r + chi_i
                _vector_add_assign_f((*r).c[0], chi0);
                _vector_add_assign_f((*r).c[1], chi1);
                _vector_add_assign_f((*r).c[2], chi2);
                _vector_add_assign_f((*r).c[3], chi3);

                // dw ------------------------------------------
                // chi_i <= U_ij^dag(n-dir) * psi_j(n-dir)
                _suNf_inverse_multiply(chi0, *(u_dw), (*psi_dw).c[0]);
                _suNf_inverse_multiply(chi1, *(u_dw), (*psi_dw).c[1]);
                _suNf_inverse_multiply(chi2, *(u_dw), (*psi_dw).c[2]);
                _suNf_inverse_multiply(chi3, *(u_dw), (*psi_dw).c[3]);

                // r <= r + chi_i
                _vector_add_assign_f((*r).c[0], chi0);
                _vector_add_assign_f((*r).c[1], chi1);
                _vector_add_assign_f((*r).c[2], chi2);
                _vector_add_assign_f((*r).c[3], chi3);
            }

            /******************************** end of loop *********************************/
            _spinor_lc_f(*r, factor, *rold, alpha * factor, *r);

#ifdef WITH_PROBE_MPI
            if (thread0 == 0) { probe_mpi(); }
#endif

        } /* SITE_FOR */
    } /* PIECE FOR */
}