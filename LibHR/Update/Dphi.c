/***************************************************************************\
* Copyright (c) 2008-2023, Claudio Pica, Sofie Martins                      *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
 *
 * File Dphi.c
 *
 * Action of the Wilson-Dirac operator D and hermitian g5D on a given
 * double-precision spinor field
 *
 *******************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "inverters.h"
#include "error.h"
#include "io.h"
#include "memory.h"
#include "utils.h"

#ifdef BC_T_SF_ROTATED
#include "update.h"
extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif /* BC_T_SF_ROTATED */

/*
 * Init of Dphi
 */

static int init_dirac = 1;
static int init_dirac_tm = 1;
static spinor_field *gtmp = NULL;
static spinor_field *etmp = NULL;
static spinor_field *otmp = NULL;
static spinor_field *otmp2 = NULL;
static spinor_field *etmp2 = NULL;

static void free_mem() {
    if (gtmp != NULL) {
        free_spinor_field(gtmp);
        gtmp = NULL;
    }
    if (etmp != NULL) {
        free_spinor_field(etmp);
        etmp = NULL;
    }
    if (etmp2 != NULL) {
        free_spinor_field(etmp2);
        etmp2 = NULL;
    }
    if (otmp != NULL) {
        free_spinor_field(otmp);
        otmp = NULL;
    }
    if (otmp2 != NULL) {
        free_spinor_field(otmp2);
        otmp2 = NULL;
    }
    init_dirac = 1;
}

static void init_Dirac() {
    if (init_dirac) {
        gtmp = alloc_spinor_field(1, &glattice);
        etmp = alloc_spinor_field(1, &glat_even);
        etmp2 = alloc_spinor_field(1, &glat_even);
        otmp = alloc_spinor_field(1, &glat_odd);
        atexit(&free_mem);
        init_dirac = 0;
    }
}
static void init_Dirac_tm() {
    if (init_dirac_tm) {
        otmp2 = alloc_spinor_field(1, &glat_odd);
        atexit(&free_mem);
        init_dirac_tm = 0;
    }
    init_Dirac();
}

/*
 * the following variable is used to keep trace of
 * matrix-vector multiplication.
 * we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter = 0;

unsigned long int getMVM_cpu() {
    unsigned long int res = MVMcounter >> 1; /* divide by two */
    MVMcounter = 0; /* reset counter */

    return res;
}

/* r=t*u*s */
#ifdef BC_T_THETA

#define _suNf_theta_T_multiply(r, u, s) \
    _suNf_multiply(vtmp[0], (u), (s));  \
    _vector_mulc_f((r), eitheta[0], vtmp[0])

#define _suNf_theta_T_inverse_multiply(r, u, s) \
    _suNf_inverse_multiply(vtmp[0], (u), (s));  \
    _vector_mulc_star_f((r), eitheta[0], vtmp[0])

#define _suNf_double_theta_T_multiply(r1, r2, u, s1, s2)      \
    _suNf_double_multiply(vtmp[0], vtmp[1], (u), (s1), (s2)); \
    _vector_mulc_f((r1), eitheta[0], vtmp[0]);                \
    _vector_mulc_f((r2), eitheta[0], vtmp[1])

#define _suNf_double_theta_T_inverse_multiply(r1, r2, u, s1, s2)      \
    _suNf_double_inverse_multiply(vtmp[0], vtmp[1], (u), (s1), (s2)); \
    _vector_mulc_star_f((r1), eitheta[0], vtmp[0]);                   \
    _vector_mulc_star_f((r2), eitheta[0], vtmp[1])

#else

#define _suNf_theta_T_multiply(r, u, s) _suNf_multiply((r), (u), (s))
#define _suNf_theta_T_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))
#define _suNf_double_theta_T_multiply(r1, r2, u, s1, s2) _suNf_double_multiply((r1), (r2), (u), (s1), (s2))
#define _suNf_double_theta_T_inverse_multiply(r1, r2, u, s1, s2) _suNf_double_inverse_multiply((r1), (r2), (u), (s1), (s2))

#endif

/* r=t*u*s */
#ifdef BC_X_THETA

#define _suNf_theta_X_multiply(r, u, s) \
    _suNf_multiply(vtmp[0], (u), (s));  \
    _vector_mulc_f((r), eitheta[1], vtmp[0])

#define _suNf_theta_X_inverse_multiply(r, u, s) \
    _suNf_inverse_multiply(vtmp[0], (u), (s));  \
    _vector_mulc_star_f((r), eitheta[1], vtmp[0])

#define _suNf_double_theta_X_multiply(r1, r2, u, s1, s2)      \
    _suNf_double_multiply(vtmp[0], vtmp[1], (u), (s1), (s2)); \
    _vector_mulc_f((r1), eitheta[1], vtmp[0]);                \
    _vector_mulc_f((r2), eitheta[1], vtmp[1])

#define _suNf_double_theta_X_inverse_multiply(r1, r2, u, s1, s2)      \
    _suNf_double_inverse_multiply(vtmp[0], vtmp[1], (u), (s1), (s2)); \
    _vector_mulc_star_f((r1), eitheta[1], vtmp[0]);                   \
    _vector_mulc_star_f((r2), eitheta[1], vtmp[1])

#else

#define _suNf_theta_X_multiply(r, u, s) _suNf_multiply((r), (u), (s))
#define _suNf_theta_X_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))
#define _suNf_double_theta_X_multiply(r1, r2, u, s1, s2) _suNf_double_multiply((r1), (r2), (u), (s1), (s2))
#define _suNf_double_theta_X_inverse_multiply(r1, r2, u, s1, s2) _suNf_double_inverse_multiply((r1), (r2), (u), (s1), (s2))
#endif

/* r=t*u*s */
#ifdef BC_Y_THETA

#define _suNf_theta_Y_multiply(r, u, s) \
    _suNf_multiply(vtmp[0], (u), (s));  \
    _vector_mulc_f((r), eitheta[2], vtmp[0])

#define _suNf_theta_Y_inverse_multiply(r, u, s) \
    _suNf_inverse_multiply(vtmp[0], (u), (s));  \
    _vector_mulc_star_f((r), eitheta[2], vtmp[0])

#define _suNf_double_theta_Y_multiply(r1, r2, u, s1, s2)      \
    _suNf_double_multiply(vtmp[0], vtmp[1], (u), (s1), (s2)); \
    _vector_mulc_f((r1), eitheta[2], vtmp[0]);                \
    _vector_mulc_f((r2), eitheta[2], vtmp[1])

#define _suNf_double_theta_Y_inverse_multiply(r1, r2, u, s1, s2)      \
    _suNf_double_inverse_multiply(vtmp[0], vtmp[1], (u), (s1), (s2)); \
    _vector_mulc_star_f((r1), eitheta[2], vtmp[0]);                   \
    _vector_mulc_star_f((r2), eitheta[2], vtmp[1])

#else

#define _suNf_theta_Y_multiply(r, u, s) _suNf_multiply((r), (u), (s))
#define _suNf_theta_Y_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))
#define _suNf_double_theta_Y_multiply(r1, r2, u, s1, s2) _suNf_double_multiply((r1), (r2), (u), (s1), (s2))
#define _suNf_double_theta_Y_inverse_multiply(r1, r2, u, s1, s2) _suNf_double_inverse_multiply((r1), (r2), (u), (s1), (s2))

#endif

/* r=t*u*s */
#ifdef BC_Z_THETA

#define _suNf_theta_Z_multiply(r, u, s) \
    _suNf_multiply(vtmp[0], (u), (s));  \
    _vector_mulc_f((r), eitheta[3], vtmp[0])

#define _suNf_theta_Z_inverse_multiply(r, u, s) \
    _suNf_inverse_multiply(vtmp[0], (u), (s));  \
    _vector_mulc_star_f((r), eitheta[3], vtmp[0])

#define _suNf_double_theta_T_multiply(r1, r2, u, s1, s2)      \
    _suNf_double_multiply(vtmp[0], vtmp[1], (u), (s1), (s2)); \
    _vector_mulc_f((r1), eitheta[3], vtmp[0]);                \
    _vector_mulc_f((r2), eitheta[3], vtmp[1])

#define _suNf_double_theta_T_inverse_multiply(r1, r2, u, s1, s2)      \
    _suNf_double_inverse_multiply(vtmp[0], vtmp[1], (u), (s1), (s2)); \
    _vector_mulc_star_f((r1), eitheta[3], vtmp[0]);                   \
    _vector_mulc_star_f((r2), eitheta[3], vtmp[1])

#else

#define _suNf_theta_Z_multiply(r, u, s) _suNf_multiply((r), (u), (s))
#define _suNf_theta_Z_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))
#define _suNf_double_theta_Z_multiply(r1, r2, u, s1, s2) _suNf_double_multiply((r1), (r2), (u), (s1), (s2))
#define _suNf_double_theta_Z_inverse_multiply(r1, r2, u, s1, s2) _suNf_double_inverse_multiply((r1), (r2), (u), (s1), (s2))

#endif

#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
#define _USING_THETA suNf_vector vtmp[2]
#else
#define _USING_THETA (void)(0)
#endif

// Macros for Dphi directions
#define DPHI_T_UP(ix, iy, in, r)                                    \
    do {                                                            \
        const suNf_spinor *sp = _FIELD_AT(in, iy);                  \
        const suNf *up = pu_gauge_f(ix, 0);                         \
        suNf_vector psi, chi, psi2, chi2;                           \
        _USING_THETA;                                               \
        _vector_add_f(psi, (*sp).c[0], (*sp).c[2]);                 \
        _vector_add_f(psi2, (*sp).c[1], (*sp).c[3]);                \
        _suNf_double_theta_T_multiply(chi, chi2, (*up), psi, psi2); \
        _vector_mul_add_assign_f((*r).c[0], -0.5, chi);             \
        _vector_mul_add_assign_f((*r).c[2], -0.5, chi);             \
        _vector_mul_add_assign_f((*r).c[1], -0.5, chi2);            \
        _vector_mul_add_assign_f((*r).c[3], -0.5, chi2);            \
    } while (0)

#define DPHI_T_DN(ix, iy, in, r)                                            \
    do {                                                                    \
        const suNf_spinor *sm = _FIELD_AT(in, iy);                          \
        const suNf *um = pu_gauge_f(iy, 0);                                 \
        suNf_vector psi, chi, psi2, chi2;                                   \
        _USING_THETA;                                                       \
        _vector_sub_f(psi, (*sm).c[0], (*sm).c[2]);                         \
        _vector_sub_f(psi2, (*sm).c[1], (*sm).c[3]);                        \
        _suNf_double_theta_T_inverse_multiply(chi, chi2, (*um), psi, psi2); \
        _vector_mul_add_assign_f((*r).c[0], -0.5, chi);                     \
        _vector_mul_sub_assign_f((*r).c[2], -0.5, chi);                     \
        _vector_mul_add_assign_f((*r).c[1], -0.5, chi2);                    \
        _vector_mul_sub_assign_f((*r).c[3], -0.5, chi2);                    \
    } while (0)

#define DPHI_X_UP(ix, iy, in, r)                                    \
    do {                                                            \
        const suNf_spinor *sp = _FIELD_AT(in, iy);                  \
        const suNf *up = pu_gauge_f(ix, 1);                         \
        suNf_vector psi, chi, psi2, chi2;                           \
        _USING_THETA;                                               \
        _vector_i_add_f(psi, (*sp).c[0], (*sp).c[3]);               \
        _vector_i_add_f(psi2, (*sp).c[1], (*sp).c[2]);              \
        _suNf_double_theta_X_multiply(chi, chi2, (*up), psi, psi2); \
        _vector_mul_add_assign_f((*r).c[0], -0.5, chi);             \
        _vector_i_mul_sub_assign_f((*r).c[3], -0.5, chi);           \
        _vector_mul_add_assign_f((*r).c[1], -0.5, chi2);            \
        _vector_i_mul_sub_assign_f((*r).c[2], -0.5, chi2);          \
    } while (0)

#define DPHI_X_DN(ix, iy, in, r)                                            \
    do {                                                                    \
        const suNf_spinor *sm = _FIELD_AT(in, iy);                          \
        const suNf *um = pu_gauge_f(iy, 1);                                 \
        suNf_vector psi, chi, psi2, chi2;                                   \
        _USING_THETA;                                                       \
        _vector_i_sub_f(psi, (*sm).c[0], (*sm).c[3]);                       \
        _vector_i_sub_f(psi2, (*sm).c[1], (*sm).c[2]);                      \
        _suNf_double_theta_X_inverse_multiply(chi, chi2, (*um), psi, psi2); \
        _vector_mul_add_assign_f((*r).c[0], -0.5, chi);                     \
        _vector_i_mul_add_assign_f((*r).c[3], -0.5, chi);                   \
        _vector_mul_add_assign_f((*r).c[1], -0.5, chi2);                    \
        _vector_i_mul_add_assign_f((*r).c[2], -0.5, chi2);                  \
    } while (0)

#define DPHI_Y_UP(ix, iy, in, r)                                    \
    do {                                                            \
        const suNf_spinor *sp = _FIELD_AT(in, iy);                  \
        const suNf *up = pu_gauge_f(ix, 2);                         \
        suNf_vector psi, chi, psi2, chi2;                           \
        _USING_THETA;                                               \
        _vector_add_f(psi, (*sp).c[0], (*sp).c[3]);                 \
        _vector_sub_f(psi2, (*sp).c[1], (*sp).c[2]);                \
        _suNf_double_theta_Y_multiply(chi, chi2, (*up), psi, psi2); \
        _vector_mul_add_assign_f((*r).c[0], -0.5, chi);             \
        _vector_mul_add_assign_f((*r).c[3], -0.5, chi);             \
        _vector_mul_add_assign_f((*r).c[1], -0.5, chi2);            \
        _vector_mul_sub_assign_f((*r).c[2], -0.5, chi2);            \
    } while (0)

#define DPHI_Y_DN(ix, iy, in, r)                                            \
    do {                                                                    \
        const suNf_spinor *sm = _FIELD_AT(in, iy);                          \
        const suNf *um = pu_gauge_f(iy, 2);                                 \
        suNf_vector psi, chi, psi2, chi2;                                   \
        _USING_THETA;                                                       \
        _vector_sub_f(psi, (*sm).c[0], (*sm).c[3]);                         \
        _vector_add_f(psi2, (*sm).c[1], (*sm).c[2]);                        \
        _suNf_double_theta_Y_inverse_multiply(chi, chi2, (*um), psi, psi2); \
        _vector_mul_add_assign_f((*r).c[0], -0.5, chi);                     \
        _vector_mul_sub_assign_f((*r).c[3], -0.5, chi);                     \
        _vector_mul_add_assign_f((*r).c[1], -0.5, chi2);                    \
        _vector_mul_add_assign_f((*r).c[2], -0.5, chi2);                    \
    } while (0)

#define DPHI_Z_UP(ix, iy, in, r)                                    \
    do {                                                            \
        const suNf_spinor *sp = _FIELD_AT(in, iy);                  \
        const suNf *up = pu_gauge_f(ix, 3);                         \
        suNf_vector psi, chi, psi2, chi2;                           \
        _USING_THETA;                                               \
        _vector_i_add_f(psi, (*sp).c[0], (*sp).c[2]);               \
        _vector_i_sub_f(psi2, (*sp).c[1], (*sp).c[3]);              \
        _suNf_double_theta_Z_multiply(chi, chi2, (*up), psi, psi2); \
        _vector_mul_add_assign_f((*r).c[0], -0.5, chi);             \
        _vector_i_mul_sub_assign_f((*r).c[2], -0.5, chi);           \
        _vector_mul_add_assign_f((*r).c[1], -0.5, chi2);            \
        _vector_i_mul_add_assign_f((*r).c[3], -0.5, chi2);          \
    } while (0)

#define DPHI_Z_DN(ix, iy, in, r)                                            \
    do {                                                                    \
        const suNf_spinor *sm = _FIELD_AT(in, iy);                          \
        const suNf *um = pu_gauge_f(iy, 3);                                 \
        suNf_vector psi, chi, psi2, chi2;                                   \
        _USING_THETA;                                                       \
        _vector_i_sub_f(psi, (*sm).c[0], (*sm).c[2]);                       \
        _vector_i_add_f(psi2, (*sm).c[1], (*sm).c[3]);                      \
        _suNf_double_theta_Z_inverse_multiply(chi, chi2, (*um), psi, psi2); \
        _vector_mul_add_assign_f((*r).c[0], -0.5, chi);                     \
        _vector_i_mul_add_assign_f((*r).c[2], -0.5, chi);                   \
        _vector_mul_add_assign_f((*r).c[1], -0.5, chi2);                    \
        _vector_i_mul_sub_assign_f((*r).c[3], -0.5, chi2);                  \
    } while (0)

/*
 * This function defines the massless Dirac operator
 * It can act on spinors defined on the whole lattice
 * or on spinors with definite parity
 */
#ifdef WITH_NEW_GEOMETRY
void Dphi_cpu_(spinor_field *restrict out, spinor_field *restrict in) {
#ifdef CHECK_SPINOR_MATCHING
    error((in == NULL) || (out == NULL), 1, "Dphi_cpu_ [Dphi.c]", "Attempt to access unallocated memory space");
    error(in == out, 1, "Dphi_cpu_ [Dphi.c]", "Input and output fields must be different");
    error(out->type == &glat_even && in->type == &glat_even, 1, "Dphi_cpu_ [Dphi.c]", "Spinors don't match! (1)");
    error(out->type == &glat_odd && in->type == &glat_odd, 1, "Dphi_cpu_ [Dphi.c]", "Spinors don't match! (2)");
#endif

    ++MVMcounter; /* count matrix calls */
    if (out->type == &glattice) { ++MVMcounter; }

    const int reps = (out->type->nbuffers_spinor > 0) ? 2 : 1;

    /* start communication of input spinor field */
#ifdef WITH_MPI
    start_sendrecv_spinor_field(in);
#endif

    /************************ loop over all lattice sites *************************/
    for (int repeat = 0; repeat < reps; repeat++) {
        // we repeat the loop over the master lattice twice
        // the second pass we invert the mask
        // this is achieved with comparing the condition to be different than repeat=0,1

        _MASTER_FOR(out->type, ix) {
            register int thread0 = hr_threadId();
            register suNf_spinor *r = _FIELD_AT(out, ix);
            if (repeat == 0) { _spinor_zero_f(*r); }

            /******************************* direction +0 *********************************/
            if ((!(imask[ix] & T_UP_MASK)) == repeat) {
                const int iy = iup(ix, 0);
                DPHI_T_UP(ix, iy, in, r);
            }
            /******************************* direction -0 *********************************/
            if ((!(imask[ix] & T_DN_MASK)) == repeat) {
                const int iy = idn(ix, 0);
                DPHI_T_DN(ix, iy, in, r);
            }
            /******************************* direction +1 *********************************/
            if ((!(imask[ix] & X_UP_MASK)) == repeat) {
                const int iy = iup(ix, 1);
                DPHI_X_UP(ix, iy, in, r);
            }
            /******************************* direction -1 *********************************/
            if ((!(imask[ix] & X_DN_MASK)) == repeat) {
                const int iy = idn(ix, 1);
                DPHI_X_DN(ix, iy, in, r);
            }
            /******************************* direction +2 *********************************/
            if ((!(imask[ix] & Y_UP_MASK)) == repeat) {
                const int iy = iup(ix, 2);
                DPHI_Y_UP(ix, iy, in, r);
            }
            /******************************* direction -2 *********************************/
            if ((!(imask[ix] & Y_DN_MASK)) == repeat) {
                const int iy = idn(ix, 2);
                DPHI_Y_DN(ix, iy, in, r);
            }
            /******************************* direction +3 *********************************/
            if ((!(imask[ix] & Z_UP_MASK)) == repeat) {
                const int iy = iup(ix, 3);
                DPHI_Z_UP(ix, iy, in, r);
            }
            /******************************* direction -3 *********************************/
            if ((!(imask[ix] & Z_DN_MASK)) == repeat) {
                const int iy = idn(ix, 3);
                DPHI_Z_DN(ix, iy, in, r);
            }
            /******************************** end of loop *********************************/
#ifdef WITH_PROBE_MPI
            if (thread0 == 0) { probe_mpi(); }
#endif

        } /* MASTER_FOR */

#ifdef WITH_MPI
        if (!repeat) {
            // lprintf("MAIN", 0, "Doing complete sendrecv repeat=%d\n",repeat);
            /* wait for spinor to be transfered */
            complete_sendrecv_spinor_field(in);
        }
#endif
    }
}

void Dphi_cpu_new_(spinor_field *restrict out, spinor_field *restrict in) {
#ifdef CHECK_SPINOR_MATCHING
    error((in == NULL) || (out == NULL), 1, "Dphi_cpu_ [Dphi.c]", "Attempt to access unallocated memory space");
    error(in == out, 1, "Dphi_cpu_ [Dphi.c]", "Input and output fields must be different");
    error(out->type == &glat_even && in->type == &glat_even, 1, "Dphi_cpu_ [Dphi.c]", "Spinors don't match! (1)");
    error(out->type == &glat_odd && in->type == &glat_odd, 1, "Dphi_cpu_ [Dphi.c]", "Spinors don't match! (2)");
#endif

    ++MVMcounter; /* count matrix calls */
    if (out->type == &glattice) { ++MVMcounter; }

    /************************ loop over all lattice sites *************************/
    /* start communication of input spinor field */
#ifdef WITH_MPI
    start_sendrecv_spinor_field(in);
#endif

    // we repeat the loop over the master lattice twice
    // the second pass we invert the mask
    // this is achieved with comparing the condition to be different than repeat=0,1

    _MASTER_FOR(out->type, ix) {
        register suNf_spinor *r = _FIELD_AT(out, ix);
        _spinor_zero_f(*r);

        /******************************* direction +0 *********************************/
        if (imask[ix] & T_UP_MASK) {
            const int iy = iup(ix, 0);
            DPHI_T_UP(ix, iy, in, r);
        }
        /******************************* direction -0 *********************************/
        if (imask[ix] & T_DN_MASK) {
            const int iy = idn(ix, 0);
            DPHI_T_DN(ix, iy, in, r);
        }
        /******************************* direction +1 *********************************/
        if (imask[ix] & X_UP_MASK) {
            const int iy = iup(ix, 1);
            DPHI_X_UP(ix, iy, in, r);
        }
        /******************************* direction -1 *********************************/
        if (imask[ix] & X_DN_MASK) {
            const int iy = idn(ix, 1);
            DPHI_X_DN(ix, iy, in, r);
        }
        /******************************* direction +2 *********************************/
        if (imask[ix] & Y_UP_MASK) {
            const int iy = iup(ix, 2);
            DPHI_Y_UP(ix, iy, in, r);
        }
        /******************************* direction -2 *********************************/
        if (imask[ix] & Y_DN_MASK) {
            const int iy = idn(ix, 2);
            DPHI_Y_DN(ix, iy, in, r);
        }
        /******************************* direction +3 *********************************/
        if (imask[ix] & Z_UP_MASK) {
            const int iy = iup(ix, 3);
            DPHI_Z_UP(ix, iy, in, r);
        }
        /******************************* direction -3 *********************************/
        if (imask[ix] & Z_DN_MASK) {
            const int iy = idn(ix, 3);
            DPHI_Z_DN(ix, iy, in, r);
        }
        /******************************** end of loop *********************************/

    } /* MASTER_FOR */

#ifdef WITH_MPI
    // lprintf("MAIN", 0, "Doing complete sendrecv repeat=%d\n",repeat);
    /* wait for spinor to be transfered */
    complete_sendrecv_spinor_field(in);
#endif

    // loop over receive buffers
    const geometry_descriptor *intype = in->type;
    for (int i = 0; i < (intype->nbuffers_spinor); ++i) {
        //since we don't use the geometry boxes yet
        //to figure out the direction of the buffer we look at the mask of a point
        const int r_start = intype->rbuf_start[i];
        const int r_last = r_start + intype->rbuf_len[i];
        char dir = imask[r_start];
        switch (dir) {
        case (char)T_DN_MASK: // => dir +0
            _OMP_PARALLEL_FOR
            for (int iy = r_start; iy < r_last; ++iy) {
                const int ix = idn(iy, 0);
                register suNf_spinor *r = _FIELD_AT(out, ix);
                DPHI_T_UP(ix, iy, in, r);
            }
            break;

        case (char)T_UP_MASK: // => dir -0
            _OMP_PARALLEL_FOR
            for (int iy = r_start; iy < r_last; ++iy) {
                const int ix = iup(iy, 0);
                register suNf_spinor *r = _FIELD_AT(out, ix);
                DPHI_T_DN(ix, iy, in, r);
            }
            break;

        case (char)X_DN_MASK: // => dir +1
            _OMP_PARALLEL_FOR
            for (int iy = r_start; iy < r_last; ++iy) {
                const int ix = idn(iy, 1);
                register suNf_spinor *r = _FIELD_AT(out, ix);
                DPHI_X_UP(ix, iy, in, r);
            }
            break;

        case (char)X_UP_MASK: // => dir -1
            _OMP_PARALLEL_FOR
            for (int iy = r_start; iy < r_last; ++iy) {
                const int ix = iup(iy, 1);
                register suNf_spinor *r = _FIELD_AT(out, ix);
                DPHI_X_DN(ix, iy, in, r);
            }
            break;

        case (char)Y_DN_MASK: // => dir +2
            _OMP_PARALLEL_FOR
            for (int iy = r_start; iy < r_last; ++iy) {
                const int ix = idn(iy, 2);
                register suNf_spinor *r = _FIELD_AT(out, ix);
                DPHI_Y_UP(ix, iy, in, r);
            }
            break;

        case (char)Y_UP_MASK: // => dir -2
            _OMP_PARALLEL_FOR
            for (int iy = r_start; iy < r_last; ++iy) {
                const int ix = iup(iy, 2);
                register suNf_spinor *r = _FIELD_AT(out, ix);
                DPHI_Y_DN(ix, iy, in, r);
            }
            break;

        case (char)Z_DN_MASK: // => dir +3
            _OMP_PARALLEL_FOR
            for (int iy = r_start; iy < r_last; ++iy) {
                const int ix = idn(iy, 3);
                register suNf_spinor *r = _FIELD_AT(out, ix);
                DPHI_Z_UP(ix, iy, in, r);
            }
            break;

        case (char)Z_UP_MASK: // => dir -3
            _OMP_PARALLEL_FOR
            for (int iy = r_start; iy < r_last; ++iy) {
                const int ix = iup(iy, 3);
                register suNf_spinor *r = _FIELD_AT(out, ix);
                DPHI_Z_DN(ix, iy, in, r);
            }
            break;

        default:
            error(1, 1, "Dphi_cpu_ [Dphi.c]", "Illegal direction in boundary mask!");
            break;
        }
    } //loop over recv buffers
}

#else
void Dphi_cpu_(spinor_field *restrict out, spinor_field *restrict in) {
#ifdef CHECK_SPINOR_MATCHING
    error((in == NULL) || (out == NULL), 1, "Dphi_cpu_ [Dphi.c]", "Attempt to access unallocated memory space");
    error(in == out, 1, "Dphi_cpu_ [Dphi.c]", "Input and output fields must be different");
    error(out->type == &glat_even && in->type == &glat_even, 1, "Dphi_cpu_ [Dphi.c]", "Spinors don't match! (1)");
    error(out->type == &glat_odd && in->type == &glat_odd, 1, "Dphi_cpu_ [Dphi.c]", "Spinors don't match! (2)");
#endif

    ++MVMcounter; /* count matrix calls */
    if (out->type == &glattice) { ++MVMcounter; }

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
            register int thread0 = hr_threadId();
            int iy;
            suNf *up, *um;
            suNf_vector psi, chi, psi2, chi2;
            suNf_spinor *r, *sp, *sm;
#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
            suNf_vector vtmp[2];
#endif

            r = _FIELD_AT(out, ix);

            /******************************* direction +0 *********************************/

            iy = iup(ix, 0);
            sp = _FIELD_AT(in, iy);
            up = pu_gauge_f(ix, 0);

            _vector_add_f(psi, (*sp).c[0], (*sp).c[2]);
            _vector_add_f(psi2, (*sp).c[1], (*sp).c[3]);
            _suNf_double_theta_T_multiply(chi, chi2, (*up), psi, psi2);

            (*r).c[0] = chi;
            (*r).c[2] = chi;
            (*r).c[1] = chi2;
            (*r).c[3] = chi2;

            /******************************* direction -0 *********************************/

            iy = idn(ix, 0);
            sm = _FIELD_AT(in, iy);
            um = pu_gauge_f(iy, 0);

            _vector_sub_f(psi, (*sm).c[0], (*sm).c[2]);
            _vector_sub_f(psi2, (*sm).c[1], (*sm).c[3]);
            _suNf_double_theta_T_inverse_multiply(chi, chi2, (*um), psi, psi2);

            _vector_add_assign_f((*r).c[0], chi);
            _vector_sub_assign_f((*r).c[2], chi);
            _vector_add_assign_f((*r).c[1], chi2);
            _vector_sub_assign_f((*r).c[3], chi2);

            /******************************* direction +1 *********************************/

            iy = iup(ix, 1);
            sp = _FIELD_AT(in, iy);
            up = pu_gauge_f(ix, 1);

            _vector_i_add_f(psi, (*sp).c[0], (*sp).c[3]);
            _vector_i_add_f(psi2, (*sp).c[1], (*sp).c[2]);
            _suNf_double_theta_X_multiply(chi, chi2, (*up), psi, psi2);

            _vector_add_assign_f((*r).c[0], chi);
            _vector_i_sub_assign_f((*r).c[3], chi);
            _vector_add_assign_f((*r).c[1], chi2);
            _vector_i_sub_assign_f((*r).c[2], chi2);

            /******************************* direction -1 *********************************/

            iy = idn(ix, 1);
            sm = _FIELD_AT(in, iy);
            um = pu_gauge_f(iy, 1);

            _vector_i_sub_f(psi, (*sm).c[0], (*sm).c[3]);
            _vector_i_sub_f(psi2, (*sm).c[1], (*sm).c[2]);
            _suNf_double_theta_X_inverse_multiply(chi, chi2, (*um), psi, psi2);

            _vector_add_assign_f((*r).c[0], chi);
            _vector_i_add_assign_f((*r).c[3], chi);
            _vector_add_assign_f((*r).c[1], chi2);
            _vector_i_add_assign_f((*r).c[2], chi2);

            /******************************* direction +2 *********************************/

            iy = iup(ix, 2);
            sp = _FIELD_AT(in, iy);
            up = pu_gauge_f(ix, 2);

            _vector_add_f(psi, (*sp).c[0], (*sp).c[3]);
            _vector_sub_f(psi2, (*sp).c[1], (*sp).c[2]);
            _suNf_double_theta_Y_multiply(chi, chi2, (*up), psi, psi2);

            _vector_add_assign_f((*r).c[0], chi);
            _vector_add_assign_f((*r).c[3], chi);
            _vector_add_assign_f((*r).c[1], chi2);
            _vector_sub_assign_f((*r).c[2], chi2);

            /******************************* direction -2 *********************************/

            iy = idn(ix, 2);
            sm = _FIELD_AT(in, iy);
            um = pu_gauge_f(iy, 2);

            _vector_sub_f(psi, (*sm).c[0], (*sm).c[3]);
            _vector_add_f(psi2, (*sm).c[1], (*sm).c[2]);
            _suNf_double_theta_Y_inverse_multiply(chi, chi2, (*um), psi, psi2);

            _vector_add_assign_f((*r).c[0], chi);
            _vector_sub_assign_f((*r).c[3], chi);
            _vector_add_assign_f((*r).c[1], chi2);
            _vector_add_assign_f((*r).c[2], chi2);

            /******************************* direction +3 *********************************/

            iy = iup(ix, 3);
            sp = _FIELD_AT(in, iy);
            up = pu_gauge_f(ix, 3);

            _vector_i_add_f(psi, (*sp).c[0], (*sp).c[2]);
            _vector_i_sub_f(psi2, (*sp).c[1], (*sp).c[3]);
            _suNf_double_theta_Z_multiply(chi, chi2, (*up), psi, psi2);

            _vector_add_assign_f((*r).c[0], chi);
            _vector_i_sub_assign_f((*r).c[2], chi);
            _vector_add_assign_f((*r).c[1], chi2);
            _vector_i_add_assign_f((*r).c[3], chi2);

            /******************************* direction -3 *********************************/

            iy = idn(ix, 3);
            sm = _FIELD_AT(in, iy);
            um = pu_gauge_f(iy, 3);

            _vector_i_sub_f(psi, (*sm).c[0], (*sm).c[2]);
            _vector_i_add_f(psi2, (*sm).c[1], (*sm).c[3]);
            _suNf_double_theta_Z_inverse_multiply(chi, chi2, (*um), psi, psi2);

            _vector_add_assign_f((*r).c[0], chi);
            _vector_i_add_assign_f((*r).c[2], chi);
            _vector_add_assign_f((*r).c[1], chi2);
            _vector_i_sub_assign_f((*r).c[3], chi2);

            /******************************** end of loop *********************************/
            _spinor_mul_f(*r, -0.5, *r);
#ifdef WITH_PROBE_MPI
            if (thread0 == 0) { probe_mpi(); }
#endif

        } /* SITE_FOR */
    } /* PIECE FOR */
}

#if (NG == 3) && defined(REPR_FUNDAMENTAL)
void Dphi_fused_(spinor_field *restrict out, spinor_field *restrict in) {
#ifdef CHECK_SPINOR_MATCHING
    error((in == NULL) || (out == NULL), 1, "Dphi_fused_ [Dphi.c]", "Attempt to access unallocated memory space");
    error(in == out, 1, "Dphi_fused_ [Dphi.c]", "Input and output fields must be different");
    error(out->type == &glat_even && in->type == &glat_even, 1, "Dphi_fused_ [Dphi.c]", "Spinors don't match! (1)");
    error(out->type == &glat_odd && in->type == &glat_odd, 1, "Dphi_fused_ [Dphi.c]", "Spinors don't match! (2)");
#endif

    ++MVMcounter; /* count matrix calls */
    if (out->type == &glattice) { ++MVMcounter; }
    //
    /************************ loop over all lattice sites *************************/
    /* start communication of input spinor field */
    _OMP_PRAGMA(master) {
        start_sendrecv_spinor_field(in);
    }

    int iy;
    suNf *up, *um;
    suNf_vector psi, chi, psi2, chi2;
    suNf_spinor *r, *sp, *sm;
#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
    suNf_vector vtmp[2];
#endif

    _OMP_PRAGMA(_omp_for nowait)
    for (int _fuse_master_for_ip_ix = 0; _fuse_master_for_ip_ix < (out->type)->fuse_inner_counter; _fuse_master_for_ip_ix++) {
        register int thread0 = hr_threadId();
        _FUSE_IDX(out->type, ix);

        r = _FIELD_AT(out, ix);

        /******************************* direction +0 *********************************/

        iy = iup(ix, 0);
        sp = _FIELD_AT(in, iy);
        up = pu_gauge_f(ix, 0);

        _vector_add_f(psi, (*sp).c[0], (*sp).c[2]);
        _vector_add_f(psi2, (*sp).c[1], (*sp).c[3]);
        _suNf_double_theta_T_multiply(chi, chi2, (*up), psi, psi2);

        (*r).c[0] = chi;
        (*r).c[2] = chi;
        (*r).c[1] = chi2;
        (*r).c[3] = chi2;

        /******************************* direction -0 *********************************/

        iy = idn(ix, 0);
        sm = _FIELD_AT(in, iy);
        um = pu_gauge_f(iy, 0);

        _vector_sub_f(psi, (*sm).c[0], (*sm).c[2]);
        _vector_sub_f(psi2, (*sm).c[1], (*sm).c[3]);
        _suNf_double_theta_T_inverse_multiply(chi, chi2, (*um), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_sub_assign_f((*r).c[2], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_sub_assign_f((*r).c[3], chi2);

        /******************************* direction +1 *********************************/

        iy = iup(ix, 1);
        sp = _FIELD_AT(in, iy);
        up = pu_gauge_f(ix, 1);

        _vector_i_add_f(psi, (*sp).c[0], (*sp).c[3]);
        _vector_i_add_f(psi2, (*sp).c[1], (*sp).c[2]);
        _suNf_double_theta_X_multiply(chi, chi2, (*up), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_i_sub_assign_f((*r).c[3], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_i_sub_assign_f((*r).c[2], chi2);

        /******************************* direction -1 *********************************/

        iy = idn(ix, 1);
        sm = _FIELD_AT(in, iy);
        um = pu_gauge_f(iy, 1);

        _vector_i_sub_f(psi, (*sm).c[0], (*sm).c[3]);
        _vector_i_sub_f(psi2, (*sm).c[1], (*sm).c[2]);
        _suNf_double_theta_X_inverse_multiply(chi, chi2, (*um), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_i_add_assign_f((*r).c[3], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_i_add_assign_f((*r).c[2], chi2);

        /******************************* direction +2 *********************************/

        iy = iup(ix, 2);
        sp = _FIELD_AT(in, iy);
        up = pu_gauge_f(ix, 2);

        _vector_add_f(psi, (*sp).c[0], (*sp).c[3]);
        _vector_sub_f(psi2, (*sp).c[1], (*sp).c[2]);
        _suNf_double_theta_Y_multiply(chi, chi2, (*up), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_add_assign_f((*r).c[3], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_sub_assign_f((*r).c[2], chi2);

        /******************************* direction -2 *********************************/

        iy = idn(ix, 2);
        sm = _FIELD_AT(in, iy);
        um = pu_gauge_f(iy, 2);

        _vector_sub_f(psi, (*sm).c[0], (*sm).c[3]);
        _vector_add_f(psi2, (*sm).c[1], (*sm).c[2]);
        _suNf_double_theta_Y_inverse_multiply(chi, chi2, (*um), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_sub_assign_f((*r).c[3], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_add_assign_f((*r).c[2], chi2);

        /******************************* direction +3 *********************************/

        iy = iup(ix, 3);
        sp = _FIELD_AT(in, iy);
        up = pu_gauge_f(ix, 3);

        _vector_i_add_f(psi, (*sp).c[0], (*sp).c[2]);
        _vector_i_sub_f(psi2, (*sp).c[1], (*sp).c[3]);
        _suNf_double_theta_Z_multiply(chi, chi2, (*up), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_i_sub_assign_f((*r).c[2], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_i_add_assign_f((*r).c[3], chi2);

        /******************************* direction -3 *********************************/

        iy = idn(ix, 3);
        sm = _FIELD_AT(in, iy);
        um = pu_gauge_f(iy, 3);

        _vector_i_sub_f(psi, (*sm).c[0], (*sm).c[2]);
        _vector_i_add_f(psi2, (*sm).c[1], (*sm).c[3]);
        _suNf_double_theta_Z_inverse_multiply(chi, chi2, (*um), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_i_add_assign_f((*r).c[2], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_i_sub_assign_f((*r).c[3], chi2);

        /******************************** end of loop *********************************/

        _spinor_mul_f(*r, -0.5, *r);

#ifdef WITH_PROBE_MPI
        if (thread0 == 0) { probe_mpi(); }
#endif

    } /* FUSE FOR */
#ifdef WITH_MPI

    _OMP_PRAGMA(master) {
        complete_sendrecv_spinor_field(in);
    }
    _OMP_BARRIER
#endif

    _OMP_PRAGMA(_omp_for nowait)
    for (int _fuse_master_for_ip_ix = (out->type)->fuse_inner_counter; _fuse_master_for_ip_ix < (out->type)->fuse_gauge_size;
         _fuse_master_for_ip_ix++) {
        _FUSE_IDX(out->type, ix);

        r = _FIELD_AT(out, ix);

        /******************************* direction +0 *********************************/

        iy = iup(ix, 0);
        sp = _FIELD_AT(in, iy);
        up = pu_gauge_f(ix, 0);

        _vector_add_f(psi, (*sp).c[0], (*sp).c[2]);
        _vector_add_f(psi2, (*sp).c[1], (*sp).c[3]);
        _suNf_double_theta_T_multiply(chi, chi2, (*up), psi, psi2);

        (*r).c[0] = chi;
        (*r).c[2] = chi;
        (*r).c[1] = chi2;
        (*r).c[3] = chi2;

        /******************************* direction -0 *********************************/

        iy = idn(ix, 0);
        sm = _FIELD_AT(in, iy);
        um = pu_gauge_f(iy, 0);

        _vector_sub_f(psi, (*sm).c[0], (*sm).c[2]);
        _vector_sub_f(psi2, (*sm).c[1], (*sm).c[3]);
        _suNf_double_theta_T_inverse_multiply(chi, chi2, (*um), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_sub_assign_f((*r).c[2], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_sub_assign_f((*r).c[3], chi2);

        /******************************* direction +1 *********************************/

        iy = iup(ix, 1);
        sp = _FIELD_AT(in, iy);
        up = pu_gauge_f(ix, 1);

        _vector_i_add_f(psi, (*sp).c[0], (*sp).c[3]);
        _vector_i_add_f(psi2, (*sp).c[1], (*sp).c[2]);
        _suNf_double_theta_X_multiply(chi, chi2, (*up), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_i_sub_assign_f((*r).c[3], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_i_sub_assign_f((*r).c[2], chi2);

        /******************************* direction -1 *********************************/

        iy = idn(ix, 1);
        sm = _FIELD_AT(in, iy);
        um = pu_gauge_f(iy, 1);

        _vector_i_sub_f(psi, (*sm).c[0], (*sm).c[3]);
        _vector_i_sub_f(psi2, (*sm).c[1], (*sm).c[2]);
        _suNf_double_theta_X_inverse_multiply(chi, chi2, (*um), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_i_add_assign_f((*r).c[3], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_i_add_assign_f((*r).c[2], chi2);

        /******************************* direction +2 *********************************/

        iy = iup(ix, 2);
        sp = _FIELD_AT(in, iy);
        up = pu_gauge_f(ix, 2);

        _vector_add_f(psi, (*sp).c[0], (*sp).c[3]);
        _vector_sub_f(psi2, (*sp).c[1], (*sp).c[2]);
        _suNf_double_theta_Y_multiply(chi, chi2, (*up), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_add_assign_f((*r).c[3], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_sub_assign_f((*r).c[2], chi2);

        /******************************* direction -2 *********************************/

        iy = idn(ix, 2);
        sm = _FIELD_AT(in, iy);
        um = pu_gauge_f(iy, 2);

        _vector_sub_f(psi, (*sm).c[0], (*sm).c[3]);
        _vector_add_f(psi2, (*sm).c[1], (*sm).c[2]);
        _suNf_double_theta_Y_inverse_multiply(chi, chi2, (*um), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_sub_assign_f((*r).c[3], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_add_assign_f((*r).c[2], chi2);

        /******************************* direction +3 *********************************/

        iy = iup(ix, 3);
        sp = _FIELD_AT(in, iy);
        up = pu_gauge_f(ix, 3);

        _vector_i_add_f(psi, (*sp).c[0], (*sp).c[2]);
        _vector_i_sub_f(psi2, (*sp).c[1], (*sp).c[3]);
        _suNf_double_theta_Z_multiply(chi, chi2, (*up), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_i_sub_assign_f((*r).c[2], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_i_add_assign_f((*r).c[3], chi2);

        /******************************* direction -3 *********************************/

        iy = idn(ix, 3);
        sm = _FIELD_AT(in, iy);
        um = pu_gauge_f(iy, 3);

        _vector_i_sub_f(psi, (*sm).c[0], (*sm).c[2]);
        _vector_i_add_f(psi2, (*sm).c[1], (*sm).c[3]);
        _suNf_double_theta_Z_inverse_multiply(chi, chi2, (*um), psi, psi2);

        _vector_add_assign_f((*r).c[0], chi);
        _vector_i_add_assign_f((*r).c[2], chi);
        _vector_add_assign_f((*r).c[1], chi2);
        _vector_i_sub_assign_f((*r).c[3], chi2);

        /******************************** end of loop *********************************/

        _spinor_mul_f(*r, -0.5, *r);

    } /* FUSE FOR */
}
#endif
#endif
/*
 * this function takes 2 spinors defined on the whole lattice
 */
void Dphi_cpu(double m0, spinor_field *out, spinor_field *in) {
    double rho;
#ifdef BC_T_SF_ROTATED
    int ix, iy, iz, index;
    suNf_spinor *r, *sp;
    double SFrho;
    suNf_spinor tmp;
#endif /* BC_T_SF_ROTATED */

    error((in == NULL) || (out == NULL), 1, "Dphi_cpu [Dphi.c]", "Attempt to access unallocated memory space");

    error(in == out, 1, "Dphi_cpu [Dphi.c]", "Input and output fields must be different");

    apply_BCs_on_spinor_field(in);

#ifdef CHECK_SPINOR_MATCHING
    error(out->type != &glattice || in->type != &glattice, 1, "Dphi_cpu [Dphi.c]",
          "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

    Dphi_cpu_(out, in);

    rho = 4. + m0;

    //mul_add_assign(out, rho, in);
    mul_add_assign_spinor_field_cpu(out, rho, in);

#ifdef BC_T_SF_ROTATED
    SFrho = 3. * _update_par.SF_ds + _update_par.SF_zf - 4.;

    if (COORD[0] == 0) {
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    index = ipt(1, ix, iy, iz);
                    r = _FIELD_AT(out, index);
                    sp = _FIELD_AT(in, index);
                    _spinor_mul_add_assign_f(*r, SFrho, *sp);

                    _spinor_pminus_f(tmp, *sp);
                    _spinor_g5_assign_f_cpu(tmp);
                    if (_update_par.SF_sign == 1) {
                        _spinor_i_add_assign_f(*r, tmp);
                    } else {
                        _spinor_i_sub_assign_f(*r, tmp);
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    index = ipt(T - 1, ix, iy, iz);
                    r = _FIELD_AT(out, index);
                    sp = _FIELD_AT(in, index);
                    _spinor_mul_add_assign_f(*r, SFrho, *sp);

                    _spinor_pplus_f(tmp, *sp);
                    _spinor_g5_assign_f_cpu(tmp);
                    if (_update_par.SF_sign == 1) {
                        _spinor_i_add_assign_f(*r, tmp);
                    } else {
                        _spinor_i_sub_assign_f(*r, tmp);
                    }
                }
            }
        }
    }
#endif /* BC_T_SF_ROTATED */

    apply_BCs_on_spinor_field(out);
}

void g5Dphi_cpu(double m0, spinor_field *out, spinor_field *in) {
    double rho;
#ifdef BC_T_SF_ROTATED
    int ix, iy, iz, index;
    suNf_spinor *r, *sp;
    double SFrho;
    suNf_spinor tmp;
#endif /* BC_T_SF_ROTATED */

    error((in == NULL) || (out == NULL), 1, "g5Dphi_cpu [Dphi.c]", "Attempt to access unallocated memory space");

    error(in == out, 1, "g5Dphi_cpu [Dphi.c]", "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
    error(out->type != &glattice || in->type != &glattice, 1, "g5Dphi_cpu [Dphi.c]",
          "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

    apply_BCs_on_spinor_field(in);
    Dphi_cpu_(out, in);
    rho = 4. + m0;
    mul_add_assign_spinor_field_cpu(out, rho, in);

#ifdef BC_T_SF_ROTATED
    SFrho = 3. * _update_par.SF_ds + _update_par.SF_zf - 4.;

    if (COORD[0] == 0) {
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    index = ipt(1, ix, iy, iz);
                    r = _FIELD_AT(out, index);
                    sp = _FIELD_AT(in, index);
                    _spinor_mul_add_assign_f(*r, SFrho, *sp);

                    _spinor_pminus_f(tmp, *sp);
                    _spinor_g5_assign_f_cpu(tmp);
                    if (_update_par.SF_sign == 1) {
                        _spinor_i_add_assign_f_cpu(*r, tmp);
                    } else {
                        _spinor_i_sub_assign_f(*r, tmp);
                    }
                }
            }
        }
    }
    if (COORD[0] == NP_T - 1) {
        for (ix = 0; ix < X; ++ix) {
            for (iy = 0; iy < Y; ++iy) {
                for (iz = 0; iz < Z; ++iz) {
                    index = ipt(T - 1, ix, iy, iz);
                    r = _FIELD_AT(out, index);
                    sp = _FIELD_AT(in, index);
                    _spinor_mul_add_assign_f_cpu(*r, SFrho, *sp);

                    _spinor_pplus_f(tmp, *sp);
                    _spinor_g5_assign_f_cpu(tmp);
                    if (_update_par.SF_sign == 1) {
                        _spinor_i_add_assign_f_cpu(*r, tmp);
                    } else {
                        _spinor_i_sub_assign_f(*r, tmp);
                    }
                }
            }
        }
    }
#endif /* BC_T_SF_ROTATED */

    g5_assign_spinor_field_cpu(out);

    apply_BCs_on_spinor_field(out);
}

/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the even lattice
 * Dphi in = (4+m0)^2*in - D_EO D_OE in
 *
 */
void Dphi_eopre_cpu(double m0, spinor_field *out, spinor_field *in) {
    double rho;

    error((in == NULL) || (out == NULL), 1, "Dphi_eopre_cpu [Dphi.c]", "Attempt to access unallocated memory space");

    error(in == out, 1, "Dphi_eopre_cpu [Dphi.c]", "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
    error(out->type != &glat_even || in->type != &glat_even, 1, "Dphi_eopre_cpu " __FILE__,
          "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

    apply_BCs_on_spinor_field(in);

    /* alloc memory for temporary spinor field */
    if (init_dirac) { init_Dirac(); }

    Dphi_cpu_(otmp, in);
    apply_BCs_on_spinor_field(otmp);
    Dphi_cpu_(out, otmp);

    rho = 4.0 + m0;
    rho *= -rho; /* this minus sign is taken into account below */

    mul_add_assign_spinor_field_cpu(out, rho, in);
    minus_spinor_field_cpu(out, out);
    apply_BCs_on_spinor_field(out);
}

/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the odd lattice
 * Dphi in = (4+m0)^2*in - D_OE D_EO in
 *
 */
void Dphi_oepre_cpu(double m0, spinor_field *out, spinor_field *in) {
    double rho;

    error((in == NULL) || (out == NULL), 1, "Dphi_oepre_cpu [Dphi.c]", "Attempt to access unallocated memory space");

    error(in == out, 1, "Dphi_oepre_cpu [Dphi.c]", "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
    error(out->type != &glat_odd || in->type != &glat_odd, 1, "Dphi_oepre_cpu " __FILE__,
          "Spinors are not defined on odd lattice!");
#endif /* CHECK_SPINOR_MATCHING */

    apply_BCs_on_spinor_field(in);

    /* alloc memory for temporary spinor field */
    if (init_dirac) { init_Dirac(); }

    Dphi_cpu_(etmp, in);
    apply_BCs_on_spinor_field(etmp);
    Dphi_cpu_(out, etmp);

    rho = 4.0 + m0;

    rho *= -rho; /* this minus sign is taken into account below */

    mul_add_assign_spinor_field_cpu(out, rho, in);
    minus_spinor_field_cpu(out, out);

    apply_BCs_on_spinor_field(out);
}

void g5Dphi_eopre_cpu(double m0, spinor_field *out, spinor_field *in) {
    double rho;

    error((in == NULL) || (out == NULL), 1, "g5Dphi_eopre_cpu [Dphi.c]", "Attempt to access unallocated memory space");

    error(in == out, 1, "Dphi_eopre_cpu [Dphi.c]", "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
    error(out->type != &glat_even || in->type != &glat_even, 1, "g5Dphi_eopre_cpu " __FILE__,
          "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

    apply_BCs_on_spinor_field(in);

    /* alloc memory for temporary spinor field */
    if (init_dirac) { init_Dirac(); }

    Dphi_cpu_(otmp, in);
    apply_BCs_on_spinor_field(otmp);
    Dphi_cpu_(out, otmp);

    rho = 4.0 + m0;

    rho *= -rho; /* this minus sign is taken into account below */

    mul_add_assign_spinor_field_cpu(out, rho, in);
    minus_spinor_field_cpu(out, out);
    g5_assign_spinor_field_cpu(out);

    apply_BCs_on_spinor_field(out);
}

/* g5Dphi_eopre ^2 */
void g5Dphi_eopre_sq_cpu(double m0, spinor_field *out, spinor_field *in) {
    /* alloc memory for temporary spinor field */
    if (init_dirac) { init_Dirac(); }

    g5Dphi_eopre_cpu(m0, etmp, in);
    g5Dphi_eopre_cpu(m0, out, etmp);
}

/* g5Dhi ^2 */
void g5Dphi_sq_cpu(double m0, spinor_field *out, spinor_field *in) {
    /* alloc memory for temporary spinor field */
    if (init_dirac) { init_Dirac(); }

#ifdef BC_T_SF_ROTATED
    /*the switch of the SF_sign is needed to take care of the antihermiticity of the boundary term of the dirac operator*/
    g5Dphi_cpu(m0, gtmp, in);
    _update_par.SF_sign = -_update_par.SF_sign;
    g5Dphi_cpu(m0, out, gtmp);
    _update_par.SF_sign = -_update_par.SF_sign;
#else
    g5Dphi_cpu(m0, gtmp, in);
    g5Dphi_cpu(m0, out, gtmp);
#endif
}

// Twisted mass operator for even odd preconditioned case
/* g5 (M^+-_ee-M_eo {M_oo}^{-1} M_oe*/
void Qhat_eopre(double m0, double mu, spinor_field *out, spinor_field *in) {
    double norm = (4 + m0) * (4 + m0) + mu * mu;
    double rho = (4 + m0) / norm;
    hr_complex imu;
    imu = -I * mu / norm;

    error((in == NULL) || (out == NULL), 1, "Qhat_eopre [Dphi.c]", "Attempt to access unallocated memory space");

    error(in == out, 1, "Qhat_eopre [Dphi.c]", "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
    error(out->type != &glat_even || in->type != &glat_even, 1, "Qhat_eopre " __FILE__,
          "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

    apply_BCs_on_spinor_field(in);

    /* alloc memory for temporary spinor field */
    if (init_dirac) { init_Dirac(); }
    Dphi_(otmp, in);
    apply_BCs_on_spinor_field(otmp);
    mul_spinor_field(otmp2, rho, otmp);
    g5_mulc_add_assign_spinor_field(otmp2, imu, otmp);
    Dphi_(out, otmp2);

    rho = -(4 + m0);
    mul_add_assign_spinor_field(out, rho, in);
    imu = -I * mu;
    g5_mulc_add_assign_spinor_field(out, imu, in);

    minus_spinor_field(out, out);
    g5_assign_spinor_field(out);

    apply_BCs_on_spinor_field(out);
}

void Qhat_eopre_sq(double m0, double mu, spinor_field *out, spinor_field *in) {
    /* alloc memory for temporary spinor field */
    if (init_dirac) { init_Dirac(); }

#ifdef BC_T_SF_ROTATED
    /*the switch of the SF_sign is needed to take care of the antihermiticity of the boundary term of the dirac operator*/
    error(1, "Qhat_eopre_sq", __FILE__, "Not implemented\n");
#else
    Qhat_eopre(m0, -mu, etmp, in);
    Qhat_eopre(m0, mu, out, etmp);
#endif
}

#ifdef WITH_CLOVER

/*************************************************
 * Dirac operators with clover term:             *
 * Cphi = Dphi + clover                          *
 * Cphi_eopre = D_ee - D_eo D_oo^-1 D_oe         *
 * Cphi_diag = D_oo or D_ee                      *
 * Cphi_diag_inv = D_oo^-1 or D_ee^-1            *
 *************************************************/

void Cphi_cpu_(double mass, spinor_field *dptr, spinor_field *sptr, int assign) {
    compute_clover_term_cpu();
    // Correct mass term
    mass = (4. + mass);
    // Loop over local sites
    _MASTER_FOR(dptr->type, ix) {
        suNf_vector v1, v2;
        suNf_spinor *out, *in, tmp;
        suNfc *s0, *s1, *s2, *s3;

        // Field pointers
        out = _FIELD_AT(dptr, ix);
        in = _FIELD_AT(sptr, ix);
        s0 = _4FIELD_AT(cl_term, ix, 0);
        s1 = _4FIELD_AT(cl_term, ix, 1);
        s2 = _4FIELD_AT(cl_term, ix, 2);
        s3 = _4FIELD_AT(cl_term, ix, 3);

        // Component 0
        _suNfc_multiply(v1, *s0, in->c[0]);
        _suNfc_multiply(v2, *s1, in->c[1]);
        _vector_add_f(tmp.c[0], v1, v2);

        // Component 1
        _suNfc_inverse_multiply(v1, *s1, in->c[0]);
        _suNfc_multiply(v2, *s0, in->c[1]);
        _vector_sub_f(tmp.c[1], v1, v2);

        // Component 2
        _suNfc_multiply(v1, *s2, in->c[2]);
        _suNfc_multiply(v2, *s3, in->c[3]);
        _vector_add_f(tmp.c[2], v1, v2);

        // Component 3
        _suNfc_inverse_multiply(v1, *s3, in->c[2]);
        _suNfc_multiply(v2, *s2, in->c[3]);
        _vector_sub_f(tmp.c[3], v1, v2);

        // Add mass
        _spinor_mul_add_assign_f(tmp, mass, *in);

        // Store
        if (assign) {
            _spinor_add_assign_f(*out, tmp);
        } else {
            *out = tmp;
        }
    }
}

void Cphi_inv_cpu_(double mass, spinor_field *dptr, spinor_field *sptr, int assign) {
    compute_clover_term_cpu();
    int N = 2 * NF;
    mass = (4. + mass);

    // Update LDL decomposition
    compute_ldl_decomp_cpu(mass);

    // Loop over local sites
    _MASTER_FOR(dptr->type, ix) {
        hr_complex *up, *dn, *x, c;
        suNf_spinor *out, *in, tmp;
        int n;

        // Field pointers
        up = _FIELD_AT(cl_ldl, ix)->up;
        dn = _FIELD_AT(cl_ldl, ix)->dn;
        out = _FIELD_AT(dptr, ix);
        in = _FIELD_AT(sptr, ix);

        // tmp = in
        tmp = *in;
        x = (hr_complex *)&tmp;

        // Forward substitution
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < i; k++) {
                n = i * (i + 1) / 2 + k;
                _complex_mul_sub_assign(x[i], up[n], x[k]);
                _complex_mul_sub_assign(x[i + N], dn[n], x[k + N]);
            }
        }

        // Backward substitution
        for (int i = N - 1; i >= 0; i--) {
            n = i * (i + 1) / 2 + i;
            _complex_mulr(x[i], 1. / creal(up[n]), x[i]);
            _complex_mulr(x[i + N], 1. / creal(dn[n]), x[i + N]);
            for (int k = i + 1; k < N; k++) {
                n = k * (k + 1) / 2 + i;
                c = conj(up[n]);
                _complex_mul_sub_assign(x[i], c, x[k]);
                c = conj(dn[n]);
                _complex_mul_sub_assign(x[i + N], c, x[k + N]);
            }
        }

        // Store
        if (assign) {
            _spinor_add_assign_f(*out, tmp);
        } else {
            *out = tmp;
        }
    }
}
// BC for clover term in open/schrödinger functional.
void Cphi(double mass, spinor_field *dptr, spinor_field *sptr) {
    apply_BCs_on_spinor_field(sptr);
    Dphi_(dptr, sptr);
    Cphi_(mass, dptr, sptr, 1);
    apply_BCs_on_spinor_field(dptr);
}

void Cphi_flt(double mass, spinor_field_flt *dptr, spinor_field_flt *sptr) {
    //TODO: Apply boundary conditions
    //TODO: CPU implementation

#ifndef WITH_GPU
    error(1, 1, __func__, "Single precision clover-improved dirac operator not implemented for CPU\n");
#else
    Dphi_flt_(dptr, sptr);
    Cphi_flt_(mass, dptr, sptr, 1);
#endif
}

void g5Cphi(double mass, spinor_field *dptr, spinor_field *sptr) {
    Cphi(mass, dptr, sptr);
    g5_assign_spinor_field(dptr);
}

void g5Cphi_sq(double mass, spinor_field *dptr, spinor_field *sptr) {
    if (init_dirac) { init_Dirac(); }

    g5Cphi(mass, gtmp, sptr);
    g5Cphi(mass, dptr, gtmp);
}

void Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr) {
    if (init_dirac) { init_Dirac(); }

    apply_BCs_on_spinor_field(sptr);
    Dphi_(otmp, sptr);
    Cphi_inv_(mass, otmp, otmp, 0);
    apply_BCs_on_spinor_field(otmp);
    Dphi_(dptr, otmp);
    minus_spinor_field(dptr, dptr);
    Cphi_(mass, dptr, sptr, 1);
    apply_BCs_on_spinor_field(dptr);
}

void g5Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr) {
    Cphi_eopre(mass, dptr, sptr);
    g5_assign_spinor_field(dptr);
}

void g5Cphi_eopre_sq(double mass, spinor_field *dptr, spinor_field *sptr) {
    if (init_dirac) { init_Dirac(); }

    g5Cphi_eopre(mass, etmp, sptr);
    g5Cphi_eopre(mass, dptr, etmp);
}

void Cphi_diag(double mass, spinor_field *dptr, spinor_field *sptr) {
    // Here (dptr == sptr) is allowed
    apply_BCs_on_spinor_field(sptr);
    Cphi_(mass, dptr, sptr, 0);
    apply_BCs_on_spinor_field(dptr);
}

void Cphi_diag_inv(double mass, spinor_field *dptr, spinor_field *sptr) {
    // Here (dptr == sptr) is allowed
    apply_BCs_on_spinor_field(sptr);
    Cphi_inv_(mass, dptr, sptr, 0);
    apply_BCs_on_spinor_field(dptr);
}

#ifdef WITH_GPU
void Cphi_diag_flt(double mass, spinor_field_flt *dptr, spinor_field_flt *sptr) {
    // Here (dptr == sptr) is allowed
    // TODO: boundary conditions other than periodic not supported yet
    //apply_BCs_on_spinor_field(sptr);
    Cphi_flt_(mass, dptr, sptr, 0);
    //apply_BCs_on_spinor_field(dptr);
}

void Cphi_diag_inv_flt(double mass, spinor_field_flt *dptr, spinor_field_flt *sptr) {
    // Here (dptr == sptr) is allowed
    // TODO: boundary conditions other than periodic not supported yet
    //apply_BCs_on_spinor_field(sptr);
    Cphi_inv_flt_(mass, dptr, sptr, 0);
    //apply_BCs_on_spinor_field(dptr);
}
#endif

#endif // #ifdef WITH_CLOVER

#ifdef WITH_EXPCLOVER

/*************************************************
 * Dirac operators with clover term:             *
 * Cphi = Dphi + exp clover                      *
 * Cphi_eopre = D_ee - D_eo D_oo^-1 D_oe         *
 * Cphi_diag = D_oo or D_ee                      *
 * Cphi_diag_inv = D_oo^-1 or D_ee^-1            *
 *************************************************/

// Inverse: 0, then normal operator; 1, then inverse operator;
// Assume hermiticity!!
void Cphi_cpu_(double mass, spinor_field *dptr, spinor_field *sptr, int assign, int inverse) {
    compute_clover_term_cpu();
    evaluate_sw_order(&mass);

    // Correct mass term
    mass = (4. + mass);

    double invexpmass = 1.0 / mass;
    if (inverse == 1) {
        invexpmass = -1.0 / mass;
        mass = 1 / mass;
    }

    // Loop over local sites
    _MASTER_FOR(dptr->type, ix) {
        suNfc Aplus[4];
        suNfc Aminus[4];

        suNfc expAplus[4];
        suNfc expAminus[4];
        suNf_vector v1, v2;
        suNf_spinor *out, *in, tmp;
        suNfc *s0, *s1, *s2, *s3;

        // Field pointers
        out = _FIELD_AT(dptr, ix);
        in = _FIELD_AT(sptr, ix);
        s0 = _4FIELD_AT(cl_term, ix, 0);
        s1 = _4FIELD_AT(cl_term, ix, 1);
        s2 = _4FIELD_AT(cl_term, ix, 2);
        s3 = _4FIELD_AT(cl_term, ix, 3);

        // Build matrices Aplus Aminus
        // Aplus  = (s0 s1, s1^dagger -s0)
        // Aminus = (s2 s3, s3^dagger -s2)

        _suNfc_mul(Aplus[0], invexpmass, *s0);
        _suNfc_mul(Aplus[1], invexpmass, *s1);
        _suNfc_dagger(Aplus[2], Aplus[1]);
        _suNfc_mul(Aplus[3], -invexpmass, *s0);

        _suNfc_mul(Aminus[0], invexpmass, *s2);
        _suNfc_mul(Aminus[1], invexpmass, *s3);
        _suNfc_dagger(Aminus[2], Aminus[1]);
        _suNfc_mul(Aminus[3], -invexpmass, *s2);

        // Exponentiate Aplus Aminus

        clover_exp(Aplus, expAplus, get_NN());
        clover_exp(Aminus, expAminus, get_NN());

        // Correct factor (4+m)

        _suNfc_mul_assign(expAplus[0], mass);
        _suNfc_mul_assign(expAplus[1], mass);
        _suNfc_mul_assign(expAplus[2], mass);
        _suNfc_mul_assign(expAplus[3], mass);
        _suNfc_mul_assign(expAminus[0], mass);
        _suNfc_mul_assign(expAminus[1], mass);
        _suNfc_mul_assign(expAminus[2], mass);
        _suNfc_mul_assign(expAminus[3], mass);

        // Apply Aplus Aminus to spinor

        // Comp 0
        _suNfc_multiply(v1, expAplus[0], in->c[0]);
        _suNfc_multiply(v2, expAplus[1], in->c[1]);
        _vector_add_f(tmp.c[0], v1, v2);

        // Comp 1
        _suNfc_multiply(v1, expAplus[2], in->c[0]);
        _suNfc_multiply(v2, expAplus[3], in->c[1]);
        _vector_add_f(tmp.c[1], v1, v2);

        // Comp 2
        _suNfc_multiply(v1, expAminus[0], in->c[2]);
        _suNfc_multiply(v2, expAminus[1], in->c[3]);
        _vector_add_f(tmp.c[2], v1, v2);

        // Comp 3
        _suNfc_multiply(v1, expAminus[2], in->c[2]);
        _suNfc_multiply(v2, expAminus[3], in->c[3]);
        _vector_add_f(tmp.c[3], v1, v2);

        // Store
        if (assign) {
            _spinor_add_assign_f(*out, tmp);
        } else {
            *out = tmp;
        }
    }
}

void Cphi(double mass, spinor_field *dptr, spinor_field *sptr) {
    apply_BCs_on_spinor_field(sptr);
    Dphi_(dptr, sptr);
    Cphi_(mass, dptr, sptr, 1, 0);
    apply_BCs_on_spinor_field(dptr);
}

#ifdef WITH_GPU
void Cphi_flt(double mass, spinor_field_flt *dptr, spinor_field_flt *sptr) {
    //apply_BCs_on_spinor_field(sptr);
    Dphi_flt_(dptr, sptr);
    Cphi_flt_(mass, dptr, sptr, 1, 0);
    //apply_BCs_on_spinor_field(dptr);
}
#endif

void g5Cphi(double mass, spinor_field *dptr, spinor_field *sptr) {
    Cphi(mass, dptr, sptr);
    g5_assign_spinor_field(dptr);
}

void g5Cphi_sq(double mass, spinor_field *dptr, spinor_field *sptr) {
    if (init_dirac) { init_Dirac(); }

    g5Cphi(mass, gtmp, sptr);
    g5Cphi(mass, dptr, gtmp);
}

void Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr) {
    if (init_dirac) { init_Dirac(); }

    apply_BCs_on_spinor_field(sptr);
    Dphi_(otmp, sptr);
    Cphi_(mass, otmp, otmp, 0, 1);
    apply_BCs_on_spinor_field(otmp);
    Dphi_(dptr, otmp);
    minus_spinor_field(dptr, dptr);
    Cphi_(mass, dptr, sptr, 1, 0);
    apply_BCs_on_spinor_field(dptr);
}

void g5Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr) {
    Cphi_eopre(mass, dptr, sptr);
    g5_assign_spinor_field(dptr);
}

void g5Cphi_eopre_sq(double mass, spinor_field *dptr, spinor_field *sptr) {
    if (init_dirac) { init_Dirac(); }

    g5Cphi_eopre(mass, etmp, sptr);
    g5Cphi_eopre(mass, dptr, etmp);
}

void Cphi_diag(double mass, spinor_field *dptr, spinor_field *sptr) {
    // Here (dptr == sptr) is allowed
    apply_BCs_on_spinor_field(sptr);
    Cphi_(mass, dptr, sptr, 0, 0);
    apply_BCs_on_spinor_field(dptr);
}

void Cphi_diag_inv(double mass, spinor_field *dptr, spinor_field *sptr) {
    // Here (dptr == sptr) is allowed
    apply_BCs_on_spinor_field(sptr);
    Cphi_(mass, dptr, sptr, 0, 1);
    apply_BCs_on_spinor_field(dptr);
}

#ifdef WITH_GPU
void Cphi_diag_flt(double mass, spinor_field_flt *dptr, spinor_field_flt *sptr) {
    // Here (dptr == sptr) is allowed
    //apply_BCs_on_spinor_field(sptr);
    Cphi_flt_(mass, dptr, sptr, 0, 0);
    //apply_BCs_on_spinor_field(dptr);
}

void Cphi_diag_inv_flt(double mass, spinor_field_flt *dptr, spinor_field_flt *sptr) {
    // Here (dptr == sptr) is allowed
    //apply_BCs_on_spinor_field(sptr);
    Cphi_flt_(mass, dptr, sptr, 0, 1);
    //apply_BCs_on_spinor_field(dptr);
}
#endif

#endif // With expclover

/*Dirac operator with twisted mass*/
void Dxx_tw_inv(double mass, double twmass, spinor_field *out, spinor_field *in, tw_D_type tw_type) {
    hr_complex z;
    double norm;
    spinor_field *aux = NULL;
    spinor_field *aux2 = NULL;
    suNf_spinor *gtmp_ptr;
#ifndef CHECK_SPINOR_MATCHING
    error(in->type == &glattice, 1, "Dxx_tw_inv [Dphi.c]", "input spinors type not correct");
    error(in->type != out->type, 1, "Dxx_tw_inv [Dphi.c]", "input spinors and output must be identical");
#endif
    if (init_dirac_tm) { init_Dirac_tm(); }

    gtmp_ptr = _PTR(gtmp);

    if (in->type == &glat_odd) {
        aux = gtmp;
        aux->type = &glat_odd;
        _PTR(aux) = _PTR(gtmp) + glat_odd.master_shift;
        aux2 = otmp2;
    } else {
        aux = etmp;
        aux2 = gtmp;
        aux2->type = &glat_even;
    }

    norm = 1. / ((4. + mass) * (4. + mass) + twmass * twmass);

    if (tw_type == DIRECT) {
        z = -I * twmass;
    } else {
        z = I * twmass;
    }

    // in= 1/((4+m)^2+mu^2)*((4+m) in +- i mu g5 in)

    g5_spinor_field(aux, in);
    mulc_spinor_field(aux2, z, aux);
    mul_spinor_field(out, (4. + mass), in);
    add_assign_spinor_field(out, aux2);
    mul_spinor_field(out, norm, out);

    gtmp->type = &glattice;
    _PTR(gtmp) = gtmp_ptr;
}

void g5Dphi_eopre_tw(double m0, double mu, spinor_field *out, spinor_field *in, tw_D_type tw_type) {
    double rho;

#ifdef CHECK_SPINOR_MATCHING
    error((in == NULL) || (out == NULL), 1, "g5Dphi_eopre_tw [Dphi.c]", "Attempt to access unallocated memory space");

    error(in == out, 1, "g5Dphi_eopre_tw [Dphi.c]", "Input and output fields must be different");

    error(in->type != &glat_even, 1, "g5Dphi_eopre_tw " __FILE__, "in spinor is not defined on even lattice!");
    error(out->type != &glat_even, 1, "g5Dphi_eopre_tw " __FILE__, "out spinor is not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

    apply_BCs_on_spinor_field(in);

    /* alloc memory for temporary spinor field */
    if (init_dirac_tm) { init_Dirac_tm(); }

    /*************************************************
   * note rho = (4. +m0)^2 + mu^2
   * D_tw_ee_(dag/direct) = rho D_tw_ee_(direct/dag)_inv
   * D_tw_ee D_tw_ee^dag = rho
   *************************************************/

    rho = 4.0 + m0;
    rho *= rho; /* this minus sign is taken into account below */
    rho += mu * mu;

    if (tw_type == DIRECT) {
        /*************************************************
     * Operators here implemented is:
     * D_pre = rho -  D_eo D_oo_tw_inv D_oe D_tw_ee^dag
     ************************************************/

        Dxx_tw_inv(m0, mu, etmp, in, DIRECT);
        mul_spinor_field(etmp, -rho, etmp);
        Dphi_(otmp, etmp);
        Dxx_tw_inv(m0, mu, otmp, otmp, DIRECT);
        Dphi_(out, otmp);
        mul_add_assign_spinor_field(out, rho, in);
    } else {
        /*************************************************
     * Operators here implemented is:
     * D_pre^dag = rho - D_ee_tw  D_eo D_oo_tw_inv^dag D_oe
     ************************************************/
        Dphi_(otmp, in);
        Dxx_tw_inv(m0, mu, otmp, otmp, DAGGER);
        Dphi_(out, otmp);
        Dxx_tw_inv(m0, mu, out, out, DAGGER);
        mul_spinor_field(out, -rho, out);
        mul_add_assign_spinor_field(out, rho, in);
    }
    g5_assign_spinor_field(out);
    apply_BCs_on_spinor_field(out);
}

void g5Dphi_eopre_tw_sq(double m0, double mu, spinor_field *out, spinor_field *in) {
    /* alloc memory for temporary spinor field */
    if (init_dirac_tm) { init_Dirac_tm(); }

    g5Dphi_eopre_tw(m0, mu, etmp2, in, DIRECT);
    g5Dphi_eopre_tw(m0, mu, out, etmp2, DAGGER);
}

#ifndef WITH_GPU
unsigned long int (*getMVM)() = getMVM_cpu;
void (*Dphi_)(spinor_field *restrict out, spinor_field *restrict in) = Dphi_cpu_;
void (*Dphi)(double m0, spinor_field *out, spinor_field *in) = Dphi_cpu;
void (*g5Dphi)(double m0, spinor_field *out, spinor_field *in) = g5Dphi_cpu;
void (*g5Dphi_sq)(double m0, spinor_field *out, spinor_field *in) = g5Dphi_sq_cpu;
void (*Dphi_eopre)(double m0, spinor_field *out, spinor_field *in) = Dphi_eopre_cpu;
void (*Dphi_oepre)(double m0, spinor_field *out, spinor_field *in) = Dphi_oepre_cpu;
void (*g5Dphi_eopre)(double m0, spinor_field *out, spinor_field *in) = g5Dphi_eopre_cpu;
void (*g5Dphi_eopre_sq)(double m0, spinor_field *out, spinor_field *in) = g5Dphi_eopre_sq_cpu;
#ifdef WITH_CLOVER
void (*Cphi_)(double mass, spinor_field *, spinor_field *, int) = Cphi_cpu_;
void (*Cphi_inv_)(double mass, spinor_field *, spinor_field *, int) = Cphi_inv_cpu_;
#endif
#ifdef WITH_EXPCLOVER
void (*Cphi_)(double mass, spinor_field *, spinor_field *, int, int) = Cphi_cpu_;
#endif
#endif // WITH_GPU
