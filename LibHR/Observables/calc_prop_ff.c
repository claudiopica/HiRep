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

//Helps QMR solver find more accurate solutions
#undef GAUSSIAN_NOISE
//#define GAUSSIAN_NOISE

static int init = 0;
static int init_odd = 0;

static mshift_par QMR_par;
static double *shift;
static double *mass;
static spinor_field *resd;
static spinor_field *tmp;
static spinor_field *tmp_odd;

//hopping = order of hopping parameter expansion
static int n_hopping = 0;

enum { _g5QMR = 0, _MINRES, _CG, _CG_4F };

void init_propagator_ff_eo(int nm, double *m, double acc) {
    int i;

    if (init == 0) {
        shift = (double *)malloc(sizeof(double) * (nm));
        mass = (double *)malloc(sizeof(double) * (nm));
        for (i = 0; i < nm; ++i) {
            mass[i] = m[i];
            shift[i] = 0;
        }
        QMR_par.n = nm;
        QMR_par.shift = shift;
        QMR_par.err2 = .5 * acc;
        QMR_par.max_iter = 0;

        set_ff_dirac_mass(mass[0]);

        resd = alloc_spinor_field(QMR_par.n, &glat_even);
        tmp = alloc_spinor_field(1, &glat_even);

        init = 1;
    }
}

void free_propagator_ff_eo() {
    error(init == 0, 1, "calc_prop_ff.c", "propagator not initialized!");

    free_spinor_field(tmp);
    free_spinor_field(resd);

    free(shift);
    free(mass);

    if (init_odd) {
        free_spinor_field(tmp_odd);
        init_odd = 0;
    }
    init = 0;
}

/***************************************************************************\

psi = D^{-1} eta

\***************************************************************************/

static void calc_propagator_ff_eo_core(spinor_field *psi, spinor_field *eta, int solver) {
    spinor_field qprop_mask;
    int i, cgiter = 0;
    error(init == 0, 1, "calc_prop_ff.c", "z2semwall method not initialized!");

    if (init_odd == 0) {
        tmp_odd = alloc_spinor_field(1, &glat_odd);
        init_odd = 1;
    }

    /* add source */
    copy_spinor_field(tmp, eta);

    for (i = 0; i < QMR_par.n; ++i) {
        zero_spinor_field(&resd[i]);
    }

    if (solver == _CG_4F) {
        spinor_field *etmp = alloc_spinor_field(1, &glat_even);
        Dff_dagger(etmp, tmp);
        cgiter += cg_mshift(&QMR_par, &Dff_sq, etmp, resd);
        free_spinor_field(etmp);
    } else {
        //The fermion matrix is not g5-hermitian, use congrad
        error(1, 1, "calc_prop.c", "Solver undefined in calc_propagator_eo_core");
    }

    for (i = 0; i < QMR_par.n; ++i) {
        /* compute solution */
        qprop_mask = psi[i];
        qprop_mask.type = &glat_even;
        copy_spinor_field(&qprop_mask, &resd[i]);
        qprop_mask.type = &glat_odd;
        qprop_mask.ptr = psi[i].ptr + glat_odd.master_shift;
        //zero_spinor_field(&qprop_mask);
        Dphi_(&qprop_mask, &resd[i]);

        spinor_sigma_pi_rho_div_assign(&qprop_mask, ff_sigma, ff_pi, (4. + mass[i]), &qprop_mask);
        minus_spinor_field(&qprop_mask, &qprop_mask);

        spinor_field *gtmp = alloc_spinor_field(1, &glattice);
        for (int k = 0; k < n_hopping; k++) {
            Dphi_(gtmp, &psi[i]);
            spinor_sigma_pi_rho_div_assign(&psi[i], ff_sigma, ff_pi, (4. + mass[i]), gtmp);
        }
        free_spinor_field(gtmp);

        if (i & 1) { ++cgiter; /* count only half of calls. Works because the number of sources is even */ }
    }
    start_sendrecv_spinor_field(psi);
    complete_sendrecv_spinor_field(psi);
    lprintf("CALC_PROP", 10, "QMR_eo MVM = %d\n", cgiter);
}

static void calc_propagator_ff_oe_core(spinor_field *psi, spinor_field *eta, int solver) {
    spinor_field qprop_mask;
    spinor_field qprop_mask_eta;
    int i, cgiter = 0;
    error(init == 0, 1, "calc_prop_ff.c", "z2semwall method not initialized!");

    if (init_odd == 0) {
        tmp_odd = alloc_spinor_field(1, &glat_odd);
        init_odd = 1;
    }

    /* add source */
    qprop_mask_eta = eta[0];
    qprop_mask_eta.type = &glat_odd;
    qprop_mask_eta.ptr = eta[0].ptr + glat_odd.master_shift;

    spinor_field *etmp = alloc_spinor_field(1, &glat_even);
    spinor_sigma_pi_rho_div_assign(tmp_odd, ff_sigma, ff_pi, (4. + mass[0]), &qprop_mask_eta);
    Dphi_(etmp, tmp_odd);

    qprop_mask_eta = eta[0];
    qprop_mask_eta.type = &glat_even;
    sub_spinor_field(tmp, &qprop_mask_eta, etmp);

    for (i = 0; i < QMR_par.n; ++i) {
        zero_spinor_field(&resd[i]);
    }

    if (solver == _CG_4F) {
        Dff_dagger(etmp, tmp);
        cgiter += cg_mshift(&QMR_par, &Dff_sq, etmp, resd);
    } else {
        //The fermion matrix is not g5-hermitian, use congrad
        error(1, 1, "calc_prop_ff.c", "Solver undefined in calc_propagator_ff_eo_core");
    }
    free_spinor_field(etmp);

    for (i = 0; i < QMR_par.n; ++i) {
        /* compute solution */
        zero_spinor_field(&psi[i]);
        qprop_mask = psi[i];
        qprop_mask.type = &glat_even;
        copy_spinor_field(&qprop_mask, &resd[i]);

        qprop_mask.type = &glat_odd;
        qprop_mask.ptr = psi[i].ptr + glat_odd.master_shift;
        Dphi_(&qprop_mask, &resd[i]);
        spinor_sigma_pi_rho_div_assign(&qprop_mask, ff_sigma, ff_pi, (4. + mass[i]), &qprop_mask);
        sub_spinor_field(&qprop_mask, tmp_odd, &qprop_mask);

        spinor_field *gtmp = alloc_spinor_field(1, &glattice);
        for (int k = 0; k < n_hopping; k++) {
            Dphi_(gtmp, &psi[i]);
            spinor_sigma_pi_rho_div_assign(&psi[i], ff_sigma, ff_pi, (4. + mass[i]), gtmp);
            minus_spinor_field(&psi[i], &psi[i]);
        }
        free_spinor_field(gtmp);
    }
    cgiter += QMR_par.n * (1 + n_hopping) / 2;

    start_sendrecv_spinor_field(psi);
    complete_sendrecv_spinor_field(psi);
    lprintf("CALC_PROP", 10, "QMR_eo MVM = %d\n", cgiter);
}

static void calc_propagator_ff_hopping_series_core(spinor_field *psi, spinor_field *eta) {
    int i, cgiter = 0;
    error(init == 0, 1, "calc_prop_ff.c", "z2semwall method not initialized!");

    for (i = 0; i < QMR_par.n; ++i) {
        zero_spinor_field(&psi[i]);
    }

    if (n_hopping > 0) {
        spinor_field *gtmp = alloc_spinor_field(1, &glattice);
        spinor_field *gtmp2 = alloc_spinor_field(1, &glattice);

        for (i = 0; i < QMR_par.n; ++i) {
            zero_spinor_field(gtmp);
            zero_spinor_field(gtmp2);
            double rho = 4. + mass[i];
            if (n_hopping > 0) {
                spinor_sigma_pi_rho_div_assign(gtmp, ff_sigma, ff_pi, rho, eta);
                copy_spinor_field(&psi[i], gtmp);
                for (int k = 0; k < n_hopping - 1; k++) {
                    Dphi_(gtmp2, gtmp);
                    spinor_sigma_pi_rho_div_assign(gtmp, ff_sigma, ff_pi, rho, gtmp2);
                    minus_spinor_field(gtmp, gtmp);
                    add_assign_spinor_field(&psi[i], gtmp);
                }
            }
        }
        cgiter += QMR_par.n * (n_hopping - 1) / 2;

        free_spinor_field(gtmp);
        free_spinor_field(gtmp2);
    }

    lprintf("CALC_PROP", 10, " MVM = %d\n", cgiter);

    start_sendrecv_spinor_field(psi);
    complete_sendrecv_spinor_field(psi);
}

static void calc_propagator_ff_core(spinor_field *psi, spinor_field *eta, int solver) {
    start_sendrecv_spinor_field(eta);
    complete_sendrecv_spinor_field(eta);

    spinor_field qprop_mask_eta;
    spinor_field qprop_mask_psi;
    int cgiter = 0;
    if (init_odd == 0) {
        tmp_odd = alloc_spinor_field(1, &glat_odd);
        init_odd = 1;
    }
    error(init == 0, 1, "calc_prop.c", "calc_propagator_core method not initialized!");

    /* Construct source
     eta_even' = eta_even - D_eo eta_odd
     eta_odd' = A_o eta_odd
   */

    qprop_mask_eta = *eta;
    qprop_mask_eta.type = &glat_odd;
    qprop_mask_eta.ptr = eta->ptr + glat_odd.master_shift;
    minus_spinor_field(&qprop_mask_eta, &qprop_mask_eta);
    Dphi_(tmp, &qprop_mask_eta);
    minus_spinor_field(&qprop_mask_eta, &qprop_mask_eta);

    spinor_sigma_pi_dagger_rho_div_assign(tmp, ff_sigma, ff_pi, (4. + mass[0]), &qprop_mask_eta);

    //if the solution vector is empty use zero guess
    if (sqnorm_spinor_field(psi) < 1e-28) {
        zero_spinor_field(resd);
    } else {
        copy_spinor_field(resd, psi);
    }

    if (solver == _CG_4F) {
        spinor_field *etmp = alloc_spinor_field(1, &glat_even);
        Dff_dagger(etmp, tmp);
        cgiter += cg_mshift(&QMR_par, &Dff_sq, etmp, resd);
        free_spinor_field(etmp);
    } else {
        //The fermion matrix is not g5-hermitian, use congrad
        error(1, 1, "calc_prop.c", "Solver undefined in calc_propagator_eo_core (4f)");
    }

    /* compute solution 
     psi_even = D_ee*resd_e
     psi_odd = D_oo^-1*eta_odd-D_oe resd_e
  */

    qprop_mask_psi = *psi;
    qprop_mask_psi.type = &glat_even;
    copy_spinor_field(&qprop_mask_psi, resd);

    qprop_mask_psi.type = &glat_odd;
    qprop_mask_psi.ptr = psi->ptr + glat_odd.master_shift;
    Dphi_(&qprop_mask_psi, resd);

    spinor_sigma_pi_rho_div_assign(&qprop_mask_psi, ff_sigma, ff_pi, (4. + mass[0]), &qprop_mask_psi);
    minus_spinor_field(&qprop_mask_psi, &qprop_mask_psi);

    ++cgiter; /* One whole call*/
    lprintf("CALC_PROP_CORE", 10, "QMR_eo MVM = %d\n", cgiter);

    start_sendrecv_spinor_field(psi);
    complete_sendrecv_spinor_field(psi);
}

void calc_propagator_ff(spinor_field *psi, spinor_field *eta, int ndilute) {
    int beta, i, n_masses;
    double *m;
    m = mass;
    n_masses = QMR_par.n;
    QMR_par.n = 1;
    for (beta = 0; beta < ndilute; ++beta) {
        for (i = 0; i < n_masses; ++i) {
            lprintf("CALC_PROPAGATOR", 10, "n masses=%d, mass = %g\n", n_masses, mass[0]);
            set_ff_dirac_mass(mass[0]);
            calc_propagator_ff_core(&psi[beta * n_masses + i], &eta[beta], _CG_4F);
            mass++;
        }
        mass = m;
    }
    QMR_par.n = n_masses;
}

void calc_propagator_ff_eo(spinor_field *psi, spinor_field *eta, int ndilute) {
    int beta;
    n_hopping = 0;
    lprintf("CALC_PROPAGATOR_EO", 20, "Calculating EO propagator with ndilute: %d\n", ndilute);
    for (beta = 0; beta < ndilute; ++beta) {
        calc_propagator_ff_eo_core(&psi[beta * QMR_par.n], &eta[beta], _CG_4F);
    }
}

void calc_propagator_ff_oe(spinor_field *psi, spinor_field *eta, int ndilute) {
    int beta;
    n_hopping = 0;
    lprintf("CALC_PROPAGATOR_EO", 20, "Calculating EO propagator with ndilute: %d\n", ndilute);
    for (beta = 0; beta < ndilute; ++beta) {
        calc_propagator_ff_oe_core(&psi[beta * QMR_par.n], &eta[beta], _CG_4F);
    }
}

void calc_propagator_ff_hopping_eo(spinor_field *psi, spinor_field *eta, int hopping, int ndilute) {
    int beta;
    n_hopping = hopping;
    lprintf("CALC_PROPAGATOR_EO", 20, "Calculating EO propagator with ndilute: %d,hopping %d\n", ndilute, hopping);
    for (beta = 0; beta < ndilute; ++beta) {
        calc_propagator_ff_eo_core(&psi[beta * QMR_par.n], &eta[beta], _CG_4F);
    }
}

void calc_propagator_ff_hopping_oe(spinor_field *psi, spinor_field *eta, int hopping, int ndilute) {
    int beta;
    n_hopping = hopping;
    lprintf("CALC_PROPAGATOR_EO", 20, "Calculating EO propagator with ndilute: %d,hopping %d\n", ndilute, hopping);
    for (beta = 0; beta < ndilute; ++beta) {
        calc_propagator_ff_oe_core(&psi[beta * QMR_par.n], &eta[beta], _CG_4F);
    }
}

void calc_propagator_ff_hopping_series(spinor_field *psi, spinor_field *eta, int hopping, int ndilute) {
    int beta;
    n_hopping = hopping;
    lprintf("CALC_PROPAGATOR", 20, "Calculating propagator with ndilute: %d, hopping %d\n", ndilute, hopping);
    for (beta = 0; beta < ndilute; ++beta) {
        calc_propagator_ff_hopping_series_core(&psi[beta * QMR_par.n], &eta[beta]);
    }
}
