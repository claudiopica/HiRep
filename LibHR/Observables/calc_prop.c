/***************************************************************************\
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen, Jarno Rantaharju            *
*                                                                           *
*                                                                           *
\***************************************************************************/

#include "observables.h"
#include "libhr_core.h"
#include "io.h"
#include "memory.h"
#include "utils.h"
#include "inverters.h"
#include "update.h"
#include "geometry.h"

// Helps QMR solver find more accurate solutions
#undef GAUSSIAN_NOISE
//#define GAUSSIAN_NOISE

static double hmass_pre;
static double tw_mass;

static void D_pre(spinor_field *out, spinor_field *in) {
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    Cphi_eopre(hmass_pre, out, in);
#else
    Dphi_eopre(hmass_pre, out, in);
#endif
}

static void H_pre(spinor_field *out, spinor_field *in) {
    g5Dphi_eopre(hmass_pre, out, in);
}

static void H2_pre(spinor_field *out, spinor_field *in) {
    g5Dphi_eopre_sq(hmass_pre, out, in);
}

// using g5 D g5, not the most efficient but not really crucial
static void Ddag_pre(spinor_field *out, spinor_field *in, spinor_field *ttmp) {
    copy_spinor_field(ttmp, in);
    g5_assign_spinor_field(ttmp);
    H_pre(out, ttmp);
}

static void Q_eopre_tw(spinor_field *out, spinor_field *in) {
    g5Dphi_eopre_tw(hmass_pre, tw_mass, out, in, DIRECT);
}

static void Q_eopre_tw_dag(spinor_field *out, spinor_field *in) {
    g5Dphi_eopre_tw(hmass_pre, tw_mass, out, in, DAGGER);
}

static void Q2_eopre_tw(spinor_field *out, spinor_field *in) {
    g5Dphi_eopre_tw_sq(hmass_pre, tw_mass, out, in);
}

static int init = 0;
static int init_eig = 0;
static int neigs = 0;

static mshift_par QMR_par;
static double *shift;
static double *mass;
#ifdef GAUSSIAN_NOISE
static spinor_field *QMR_noise;
static spinor_field *QMR_resdn;
#endif
static spinor_field *resd_even;
static spinor_field *tmp_even = NULL;
static spinor_field *tmp_odd = NULL;
static spinor_field *tmp_even2 = NULL;
// static spinor_field *tmp_norm = NULL;

// EVA parameters
double *eva_val;
static spinor_field *eva_vec;
static spinor_field *tmp_sf;

enum { _g5QMR = 0, _MINRES, _CG, _CG_4F };

/* Initialises the propagator, nm is the number of masses for multimass solver,
   m is the array of masses, and acc is the inverter accuracy
*/

static void init_eva(int nevt) {
    if (init_eig == 0) {
        eva_val = malloc(sizeof(double) * nevt);
        eva_vec = alloc_spinor_field(nevt + 1, &glat_even);
        tmp_sf = eva_vec + nevt;
        init_eig = 1;
    }
}

void init_propagator_eo(int nm, double *m, double acc) {
    int i;
#ifdef GAUSSIAN_NOISE
    int cgiter = 0;
    double norm;
#endif

    if (init == 0) {
        shift = (double *)malloc(sizeof(double) * (nm));
        mass = (double *)malloc(sizeof(double) * (nm));
        hmass_pre = m[0]; /* we can put any number here!!! */
        for (i = 0; i < nm; ++i) {
            mass[i] = m[i];
            shift[i] = (4. + hmass_pre) * (4. + hmass_pre) - (4. + m[i]) * (4. + m[i]);
        }
        QMR_par.n = nm;
        QMR_par.shift = shift;
        QMR_par.err2 = .5 * acc;
        QMR_par.max_iter = 0;

        resd_even = alloc_spinor_field(QMR_par.n, &glat_even);
        tmp_even = alloc_spinor_field(1, &glat_even);
        tmp_even2 = alloc_spinor_field(1, &glat_even);
        tmp_odd = alloc_spinor_field(1, &glat_odd);

#ifdef GAUSSIAN_NOISE
        QMR_noise = alloc_spinor_field(nm + 1, &glat_even);
        QMR_resdn = QMR_noise + 1;
        /* noisy background */
        gaussian_spinor_field(QMR_noise);
        norm = sqrt(sqnorm_spinor_field(QMR_noise));
        mul_spinor_field(QMR_noise, 1. / norm, QMR_noise);
#endif
        init = 1;
    }
#ifdef GAUSSIAN_NOISE
    /* invert noise */
    for (i = 0; i < QMR_par.n; ++i) {
        zero_spinor_field(&QMR_resdn[i]);
    }
    cgiter += g5QMR_mshift(&QMR_par, &D_pre, QMR_noise, QMR_resdn);
    lprintf("Z2SEMWALL", 10, "QMR_eo MVM = %d\n", cgiter);
#endif
}

void free_propagator_eo() {
    error(init == 0, 1, "free_propagator_eo.c", "propagator not initialized!");

    free(shift);
    free(mass);

#ifdef GAUSSIAN_NOISE
    free_spinor_field(QMR_noise);
#endif
    if (tmp_even != NULL) {
        free_spinor_field(tmp_even);
        tmp_even = NULL;
    }
    if (resd_even != NULL) {
        free_spinor_field(resd_even);
        resd_even = NULL;
    }
    if (tmp_odd != NULL) {
        free_spinor_field(tmp_odd);
        tmp_odd = NULL;
    }
    if (init_eig) {
        free(eva_val);
        free_spinor_field(eva_vec);
        init_eig = 0;
    }
    init = 0;
}

/***************************************************************************\

psi = D^{-1} eta

\***************************************************************************/
#if !defined(WITH_CLOVER) && !defined(WITH_EXPCLOVER)

static void calc_propagator_eo_core(spinor_field *psi, spinor_field *eta, int solver) {
#ifdef CHECK_SPINOR_MATCHING
    error(psi->type != &glattice, 1, "calc_propagator_eo_core [calc_prop.c]", "pdi type must be glattice!");
    error(eta->type == &glat_odd, 1, "calc_propagator_eo_core [calc_prop.c]", "eta type must not be glat_odd!");
#endif /* CHECK_SPINOR_MATCHING */

    spinor_field qprop_mask;
    int i, cgiter = 0;
    error(init == 0, 1, "calc_propagator_eo_core.c", "z2semwall method not initialized!");

    /* add source */
#ifdef GAUSSIAN_NOISE
    add_spinor_field(tmp_even, eta, QMR_noise);
#else
    spinor_field eta_mask;
    eta_mask.type = &glat_even;
#ifdef WITH_GPU
    eta_mask.gpu_ptr = eta->gpu_ptr;
#endif
    eta_mask.ptr = eta->ptr;
    copy_spinor_field(tmp_even, &eta_mask);
#endif

    // if the solution vector is empty use zero guess
    if (sqnorm_spinor_field(psi) < 1e-28) {
        for (i = 0; i < QMR_par.n; ++i) {
            zero_spinor_field(&resd_even[i]);
        }
    } else {
        for (i = 0; i < QMR_par.n; ++i) {
            psi[i].type = &glat_even;
            mul_spinor_field(&resd_even[i], 1 / (4. + mass[i]), &psi[i]);
            psi[i].type = &glattice;
        }
    }

    if (solver == _CG) {
        qprop_mask.type = &glat_even;
        copy_spinor_field(&qprop_mask, tmp_even);
        g5_assign_spinor_field(&qprop_mask);
        H_pre(tmp_even, &qprop_mask);
        cgiter += cg_mshift(&QMR_par, &H2_pre, tmp_even, resd_even);
    } else if (solver == _MINRES) {
        g5_spinor_field(tmp_even, tmp_even);
        cgiter += MINRES_mshift(&QMR_par, &H_pre, tmp_even, resd_even);
    } else {
        cgiter += g5QMR_mshift(&QMR_par, &D_pre, tmp_even, resd_even);
    }

    for (i = 0; i < QMR_par.n; ++i) {
#ifdef GAUSSIAN_NOISE
        sub_assign_spinor_field(&resd_even[i], &QMR_resdn[i]);
#endif
        /* compute solution */
        qprop_mask = psi[i];
        qprop_mask.type = &glat_even;
        mul_spinor_field(&qprop_mask, (4. + mass[i]), &resd_even[i]);
        qprop_mask.type = &glat_odd;
#ifdef WITH_GPU
        qprop_mask.gpu_ptr = psi[i].gpu_ptr + glat_odd.master_shift;
#endif
        qprop_mask.ptr = psi[i].ptr + glat_odd.master_shift;
        zero_spinor_field(&qprop_mask);
        Dphi_(&qprop_mask, &resd_even[i]);
        minus_spinor_field(&qprop_mask, &qprop_mask);
        if (i & 1) { ++cgiter; /* count only half of calls. works because the number of sources is even */ }
    }
    start_sendrecv_spinor_field(psi);
    complete_sendrecv_spinor_field(psi);
    lprintf("CALC_PROP", 10, "QMR_eo MVM = %d\n", cgiter);
}
#endif

static void calc_propagator_core(spinor_field *psi, spinor_field *eta, int solver) {
#ifdef CHECK_SPINOR_MATCHING
    error(psi->type != &glattice, 1, "calc_propagator_core [calc_prop.c]", "psi type must be glattice!");
    error(eta->type != &glattice, 1, "calc_propagator_core [calc_prop.c]", "eta type must be glattice!");
#endif /* CHECK_SPINOR_MATCHING */

    start_sendrecv_spinor_field(eta);
    complete_sendrecv_spinor_field(eta);

    spinor_field qprop_mask_eta;
    spinor_field qprop_mask_psi;
    int cgiter = 0;

    /* Construct source
     eta_even' = eta_even - D_eo D_oo^-1 eta_odd
   */
    qprop_mask_eta = *eta;
    qprop_mask_eta.type = &glat_odd;
#ifdef WITH_GPU
    qprop_mask_eta.gpu_ptr = eta->gpu_ptr + glat_odd.master_shift;
#endif
    qprop_mask_eta.ptr = eta->ptr + glat_odd.master_shift;
    mul_spinor_field(tmp_odd, (1. / (4. + mass[0])), &qprop_mask_eta);
    Dphi_(tmp_even, tmp_odd);
    qprop_mask_eta.type = &glat_even;
#ifdef WITH_GPU
    qprop_mask_eta.gpu_ptr = eta->gpu_ptr;
#endif
    qprop_mask_eta.ptr = eta->ptr;
    sub_spinor_field(tmp_even, &qprop_mask_eta, tmp_even);
#ifdef GAUSSIAN_NOISE
    add_assign_spinor_field(tmp_even, QMR_noise);
#endif
    // sub_spinor_field(resd_even,&qprop_mask_eta,tmp_even);

    // if the solution vector is empty use zero guess
    if (sqnorm_spinor_field(psi) < 1e-28) {
        zero_spinor_field(resd_even);
    } else {
        psi[0].type = &glat_even;
        mul_spinor_field(resd_even, 1 / (4. + mass[0]), psi);
        psi[0].type = &glattice;
    }

#ifdef GAUSSIAN_NOISE
    cgiter += g5QMR_mshift(&QMR_par, &D_pre, tmp_even, resd_even);
#else
    if (solver == _CG) {
        qprop_mask_eta.type = &glat_even;
        copy_spinor_field(&qprop_mask_eta, tmp_even);
        g5_assign_spinor_field(&qprop_mask_eta);
        H_pre(tmp_even, &qprop_mask_eta);
        cgiter += cg_mshift(&QMR_par, &H2_pre, tmp_even, resd_even);
    } else if (solver == _MINRES) {
        g5_spinor_field(tmp_even, tmp_even);
        cgiter += MINRES_mshift(&QMR_par, &H_pre, tmp_even, resd_even);
    } else {
        cgiter += g5QMR_mshift(&QMR_par, &D_pre, tmp_even, resd_even);
    }

    // g5_spinor_field(tmp_even,tmp_even);
    // cgiter+=MINRES_mshift(&QMR_par, &H_pre, tmp_even, resd_even);

    // cgiter+=g5QMR_mshift(&QMR_par, &D_pre, tmp_even, resd_even);

    // cg stuff
    /*qprop_mask_eta.type=&glat_even;
  copy_spinor_field( &qprop_mask_eta, tmp_even );
  g5_assign_spinor_field( &qprop_mask_eta );
  H_pre(tmp_even, &qprop_mask_eta);
  cgiter+=cg_mshift(&QMR_par, &H2_pre, tmp_even, resd_even);*/
#endif

#ifdef GAUSSIAN_NOISE
    sub_assign_spinor_field(resd_even, QMR_resdn);
#endif
    /* compute solution
     psi_even = D_ee*resd_e
     psi_odd = D_oo^-1*eta_odd-D_oe resd_e
  */

    qprop_mask_psi = *psi;
    qprop_mask_psi.type = &glat_even;
    mul_spinor_field(&qprop_mask_psi, (4. + mass[0]), resd_even);

    qprop_mask_psi.type = &glat_odd;
#ifdef WITH_GPU
    qprop_mask_psi.gpu_ptr = psi->gpu_ptr + glat_odd.master_shift;
#endif
    qprop_mask_psi.ptr = psi->ptr + glat_odd.master_shift;
    Dphi_(&qprop_mask_psi, resd_even);

    sub_spinor_field(&qprop_mask_psi, tmp_odd, &qprop_mask_psi);

    ++cgiter; /* One whole call*/
    lprintf("CALC_PROP_CORE", 10, "QMR_eo MVM = %d\n", cgiter);

    start_sendrecv_spinor_field(psi);
    complete_sendrecv_spinor_field(psi);
}
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)

static void calc_propagator_clover(spinor_field *dptr, spinor_field *sptr) {
#ifdef CHECK_SPINOR_MATCHING
    error(dptr->type != &glattice, 1, "calc_propagator_clover [calc_prop.c]", "dptr type must be glattice!");
    error(sptr->type == &glat_odd, 1, "calc_propagator_clover [calc_prop.c]", "sptr type must not be glat_odd!");
#endif /* CHECK_SPINOR_MATCHING */

    static spinor_field *otmp;
    static spinor_field *etmp, *stmp;
    static int local_init = 0;

    // Inverter
    mshift_par mpar;
    double loc_tmp;

    mpar.err2 = QMR_par.err2;
    mpar.max_iter = 0;
    mpar.n = 1;
    mpar.shift = &loc_tmp;
    mpar.shift[0] = 0;

    // Allocate temporary fields
    if (local_init == 0) {
        etmp = alloc_spinor_field(1, &glat_even);
        otmp = alloc_spinor_field(1, &glat_odd);
        stmp = alloc_spinor_field(1, &glattice);
        local_init = 1;
    }

    // Destination even/odd
    spinor_field dptr_e, dptr_o;
    dptr_e = *dptr;
    dptr_e.type = &glat_even;
    dptr_o = *dptr;
#ifdef WITH_GPU
    dptr_o.gpu_ptr += glat_odd.master_shift;
#endif
    dptr_o.ptr += glat_odd.master_shift;
    dptr_o.type = &glat_odd;

    // Source even/odd
    spinor_field sptr_e, sptr_o;
    sptr_e = *stmp;
    sptr_e.type = &glat_even;
    sptr_o = *stmp;
#ifdef WITH_GPU
    sptr_o.gpu_ptr += glat_odd.master_shift;
#endif
    sptr_o.ptr += glat_odd.master_shift;
    sptr_o.type = &glat_odd;

    // Handle source
    if (sptr->type == &glat_even) {
        zero_spinor_field(stmp);
        copy_spinor_field(&sptr_e, sptr);
    } else {
        copy_spinor_field(stmp, sptr);
    }

    // etmp = sptr_e - D_eo D_oo^-1 sptr_o
    Cphi_diag_inv(hmass_pre, otmp, &sptr_o);
    Dphi_(etmp, otmp);
    minus_spinor_field(etmp, etmp);
    add_assign_spinor_field(etmp, &sptr_e);
    // Call inverter
    g5QMR_mshift(&mpar, D_pre, etmp, &dptr_e);
    // dptr_o = D_oo^-1 ( sptr_o - D_oe dptr_e )
    Dphi_(&dptr_o, &dptr_e);
    minus_spinor_field(&dptr_o, &dptr_o);
    add_assign_spinor_field(&dptr_o, &sptr_o);
    Cphi_diag_inv(hmass_pre, &dptr_o, &dptr_o);
}
#endif

void calc_propagator(spinor_field *psi, spinor_field *eta, int ndilute) {
    error(init == 0, 1, "calc_propagator", "propagator not initialized!");

    int beta, i, n_masses;
    double *m;
    m = mass;
    n_masses = QMR_par.n;
    QMR_par.n = 1;
    for (beta = 0; beta < ndilute; ++beta) {
        for (i = 0; i < n_masses; ++i) {
            lprintf("CALC_PROPAGATOR", 10, "n masses=%d, mass = %g\n", n_masses, mass[0]);
            hmass_pre = mass[0];
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
            calc_propagator_clover(&psi[beta * n_masses + i], &eta[beta]);
#else
            calc_propagator_core(&psi[beta * n_masses + i], &eta[beta], _g5QMR);
#endif
            mass++;
        }
        mass = m;
    }
    QMR_par.n = n_masses;
    hmass_pre = mass[0];
}

void calc_propagator_eo(spinor_field *psi, spinor_field *eta, int ndilute) {
    error(init == 0, 1, "calc_propagator_eo", "propagator not initialized!");

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    calc_propagator(psi, eta, ndilute);
#else
    int beta;
    lprintf("CALC_PROPAGATOR_EO", 20, "Calculating EO propagator with ndilute: %d\n", ndilute);
    for (beta = 0; beta < ndilute; ++beta) {
        calc_propagator_eo_core(&psi[beta * QMR_par.n], &eta[beta], _g5QMR);
    }
#endif
}

/*Different source for each mass. Needed in sequential propagators
  with multiple masses */
void calc_propagator_multisource(spinor_field *psi, spinor_field *eta, int ndilute) {
    error(init == 0, 1, "calc_propagator_multisource", "propagator not initialized!");

    int beta, i, n_masses;
    double *m;
    m = mass;
    n_masses = QMR_par.n;
    QMR_par.n = 1;
    for (i = 0; i < n_masses; ++i) {
        hmass_pre = mass[0];
        for (beta = 0; beta < ndilute; ++beta) {
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
            calc_propagator_clover(&psi[beta * n_masses + i], &eta[beta * n_masses + i]);
#else
            calc_propagator_core(&psi[beta * n_masses + i], &eta[beta * n_masses + i], _g5QMR);
#endif
            mass++;
        }
    }
    QMR_par.n = n_masses;
    mass = m;
    hmass_pre = mass[0];
}

void eig_init(int nev, int nevt, int kmax, int maxiter, double lbnd, double omega1, double omega2) {
    if (!init_eig) { init_eva(nevt); }
    hmass_pre = mass[0];

    double max, mupp;
    int status, ie, n;
    mupp = fabs(hmass_pre + 4) + 4;
    mupp *= mupp;

    // Eigen Stuff
    int MVM = 0; // counter for matrix-vector multiplications

    max_eigval(&H2_pre, &glat_even, &max);
    // lprintf("MAIN",0,"MAXCHECK: cnfg=%e  uppbound=%e diff=%e %s\n",max,mupp,mupp-max,(mupp-max)<0?"[FAILED]":"[OK]");
    max = 1.1 * max;

    ie = eva_tuned(nev, nevt, 0, kmax, maxiter, lbnd, max, omega1, omega2, &H2_pre, eva_vec, eva_val, &status);
    MVM += status;
    while (ie != 0) { // if failed restart EVA
        lprintf("MAIN", 0, "Restarting EVA!\n");
        ie = eva_tuned(nev, nevt, 2, kmax, maxiter, lbnd, max, omega1, omega2, &H2_pre, eva_vec, eva_val, &status);
        MVM += status;
    }
    lprintf("MAIN", 0, "EVA MVM = %d\n", MVM);
    neigs = nev;

    for (n = 0; n < nev; ++n) {
        H2_pre(tmp_sf, &eva_vec[n]);
        lprintf("RESULT", 0, "Eig %d = %.15e %.15e\n", n, eva_val[n],
                prod_re_spinor_field(tmp_sf, &eva_vec[n]) / sqnorm_spinor_field(&eva_vec[n]));
    }
}

void copy_evec(int n, spinor_field *psi1, double *eval) {
    copy_spinor_field(psi1, &eva_vec[n]);
    *eval = eva_val[n];
}

void calc_deflated_propagator(spinor_field *psi, spinor_field *eta, int ndilute, int Nuse) {
    int beta, i, n_masses, n;
    double *m;
    m = mass;
    n_masses = QMR_par.n;
    QMR_par.n = 1;
    if (Nuse < 0 || Nuse > neigs) { Nuse = neigs; }
    lprintf("CALC_DEFLATED_PROPAGATOR", 10, "n masses=%d, mass = %g, neigs = %d\n", n_masses, mass[0], Nuse);
    for (i = 0; i < n_masses; ++i) {
        hmass_pre = mass[0];
        for (beta = 0; beta < ndilute; ++beta) {
            psi[beta * n_masses + i].type = &glat_even; // even guy
            eta[beta].type = &glat_even; // even guy
            zero_spinor_field(&psi[beta * n_masses + i]);
            Ddag_pre(tmp_even, &eta[beta], tmp_sf);
            mul_spinor_field(tmp_even, (4. + m[0]), tmp_even);
            for (n = 0; n < Nuse; ++n) {
                hr_complex p = prod_spinor_field(&eva_vec[n], tmp_even);
                _complex_mulr(p, (1. / eva_val[n]), p);
                mulc_add_assign_spinor_field(&psi[beta * n_masses + i], p, &eva_vec[n]);
            }
            calc_propagator_core(&psi[beta * n_masses + i], &eta[beta], _MINRES);
        }
        mass++;
    }
    QMR_par.n = n_masses;
    mass = m;
    hmass_pre = mass[0];
}

static void calc_propagator_eo_tw_core(spinor_field *psi, spinor_field *eta, int solver) {
#ifndef CHECK_SPINOR_MATCHING
    error(eta->type == &glattice, 1, "calc_propagator_eo_tw_core [calc_prop.c]", "incorrect type for the input (eta) spinor");
    error(psi->type == &glattice, 1, "calc_propagator_eo_tw_core [calc_prop.c]", "incorrect type for the input (psi) spinor");
#endif

    spinor_field qprop_mask_eta;
    spinor_field qprop_mask_psi;
    int cgiter = 0;

    /* Construct source
     eta_even' = eta_even - D_eo D_oo^-1 eta_odd
   */

    qprop_mask_eta = *eta;
    qprop_mask_eta.type = &glat_odd;
#ifdef WITH_GPU
    qprop_mask_eta.gpu_ptr = eta->gpu_ptr + glat_odd.master_shift;
#endif
    qprop_mask_eta.ptr = eta->ptr + glat_odd.master_shift;

    Dxx_tw_inv(mass[0], tw_mass, tmp_odd, &qprop_mask_eta, DIRECT);

    Dphi_(tmp_even, tmp_odd);
    qprop_mask_eta.type = &glat_even;
#ifdef WITH_GPU
    qprop_mask_eta.gpu_ptr = eta->gpu_ptr;
#endif
    qprop_mask_eta.ptr = eta->ptr;

    sub_assign_spinor_field(tmp_even, &qprop_mask_eta);
    // Note tmp_even = - eta_even'
    // if the solution vector is empty use zero guess

    if (sqnorm_spinor_field(psi) < 1e-28) {
        zero_spinor_field(resd_even);
    } else {
        psi[0].type = &glat_even;
        Dxx_tw_inv(mass[0], tw_mass, resd_even, psi, DAGGER);
        psi[0].type = &glattice;
    }

    if (solver == _CG) {
        qprop_mask_eta.type = &glat_even;

        g5_spinor_field(tmp_even2, tmp_even);

        Q_eopre_tw_dag(tmp_even, tmp_even2);

        cgiter += cg_mshift(&QMR_par, &Q2_eopre_tw, tmp_even, resd_even);
    } else if (solver == _MINRES) {
        g5_spinor_field(tmp_even, tmp_even);
        cgiter += MINRES_mshift(&QMR_par, &Q_eopre_tw, tmp_even, resd_even);
    } else {
        error(1 == 1, 0, "[calc_propagator_eo_tw_core]", "Invereter not implemented");
    }

    /* compute solution
     psi_even = D_ee^dag(-resd_e)
     psi_odd = D_oo^-1*eta_odd-Doo^-1 D_oe psi_even
  */

    Dxx_tw_inv(mass[0], tw_mass, resd_even, resd_even, DIRECT);

    qprop_mask_psi = *psi;
    qprop_mask_psi.type = &glat_even;
    mul_spinor_field(&qprop_mask_psi, -((4. + mass[0]) * (4. + mass[0]) + tw_mass * tw_mass), resd_even);

    qprop_mask_psi.type = &glat_odd;
#ifdef WITH_GPU
    qprop_mask_psi.gpu_ptr = psi->gpu_ptr + glat_odd.master_shift;
#endif
    qprop_mask_psi.ptr = psi->ptr + glat_odd.master_shift;
    qprop_mask_eta.type = &glat_odd;
#ifdef WITH_GPU
    qprop_mask_eta.gpu_ptr = eta->gpu_ptr + glat_odd.master_shift;
#endif
    qprop_mask_eta.ptr = eta->ptr + glat_odd.master_shift;

    Dxx_tw_inv(mass[0], tw_mass, &qprop_mask_psi, &qprop_mask_eta, DIRECT);

    qprop_mask_psi.type = &glat_even;
#ifdef WITH_GPU
    qprop_mask_psi.gpu_ptr = psi->gpu_ptr;
#endif
    qprop_mask_psi.ptr = psi->ptr;
    Dphi_(tmp_odd, &qprop_mask_psi);
    Dxx_tw_inv(mass[0], tw_mass, tmp_odd, tmp_odd, DIRECT);

    qprop_mask_psi.type = &glat_odd;
#ifdef WITH_GPU
    qprop_mask_psi.gpu_ptr = psi->gpu_ptr + glat_odd.master_shift;
#endif
    qprop_mask_psi.ptr = psi->ptr + glat_odd.master_shift;

    sub_assign_spinor_field(&qprop_mask_psi, tmp_odd);

    ++cgiter; /* One whole call*/
    lprintf("CALC_PROPAGATOR_EO_TW_CORE", 0, "QMR_eo MVM = %d\n", cgiter);

    start_sendrecv_spinor_field(psi);
    complete_sendrecv_spinor_field(psi);

#ifndef CHECK_SPINOR_MATCHING
    error(psi->type == &glattice, 1, "calc_propagator_eo_tw_core [calc_prop.c]", "incorrect type for the input (psi) spinor");
#endif
}
void calc_propagator_tw(double *lmass, double mu, spinor_field *psi, spinor_field *eta, int ndilute) {
    error(init == 0, 1, "calc_propagator_tw", "propagator not initialized!");

    int beta, i, n_masses;
    double *m;
    m = lmass;
    n_masses = QMR_par.n;
    QMR_par.n = 1;
    tw_mass = mu;

#ifndef CHECK_SPINOR_MATCHING
    error(eta->type == &glattice, 1, "calc_propagator_tw [calc_prop.c]", "incorrect type for the input (eta) spinor");
    error(psi->type == &glattice, 1, "calc_propagator_tw [calc_prop.c]", "incorrect type for the input (psi) spinor");
#endif

#ifdef WITH_CLOVER
    lprintf("CALC_PROPAGATOR_TW", 0,
            "Calc propagator has no csw improved implementation. The TM Dirac Operator is automatically O(a^2) improved ");
#endif
#ifdef GAUSSIAN_NOISE
    error(1 == 1, 1, "calc_propagator_tw [calc_prop.c]", "Gaussian noise has not been implemented");
#endif

    for (beta = 0; beta < ndilute; ++beta) {
        start_sendrecv_spinor_field(eta + beta);
        complete_sendrecv_spinor_field(eta + beta);
        for (i = 0; i < n_masses; ++i) {
            lprintf("CALC_PROPAGATOR_TW", 10, "n masses=%d, mass = %g\n", n_masses, m[i]);
            hmass_pre = m[i];
            calc_propagator_eo_tw_core(&psi[beta * n_masses + i], &eta[beta], _CG);

            m++;
        }
        m = lmass;
    }
    QMR_par.n = n_masses;
    hmass_pre = m[0];
}
#undef GAUSSIAN_NOISE
