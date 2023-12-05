/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica, Jarno Rantaharju      *
 * All rights reserved.                                                   *
 \***************************************************************************/

//Monomial defining a fermion field that interacts with the gauge, if it exists,
//and trough a four fermions interaction. The monomial four_fermion needs to be
//included to define the auxiliary fields.

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "random.h"
#include "Inverters/linear_algebra.h"

static spinor_field *tmp_pf = NULL;
static int mon_init = 1;

static void hmc_ff_gaussian_pf(monomial const *m) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);
    gaussian_spinor_field(par->pf);
}

static void hmc_ff_correct_pf(monomial const *m) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);

    /* compute H2^{1/2}*pf = H*pf */
    g5_spinor_field(tmp_pf, par->pf);
    set_ff_dirac_mass(par->mass);
    Dff_dagger(par->pf, tmp_pf);
}

static void hmc_ff_correct_la_pf(monomial const *m) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);
    double shift;

    mshift_par cgpar;
    shift = 0.;
    cgpar.n = 1;
    cgpar.shift = &shift;
    cgpar.max_iter = 0;
    cgpar.err2 = m->data.MT_prec;

    set_ff_dirac_mass(par->mass);
    zero_spinor_field(tmp_pf);
    cg_mshift(&cgpar, &Dff_sq, par->pf, tmp_pf);
    Dff(par->pf, tmp_pf);
    g5_assign_spinor_field(par->pf);
}

static const spinor_field *hmc_ff_pseudofermion(monomial const *m) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);
    return par->pf;
}

static void hmc_ff_add_local_action(monomial const *m, scalar_field *loc_action) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);
    pf_local_action(loc_action, par->pf);

#ifdef UPDATE_EO
    /* If EO preconditioning is used, there the odd diagonal part of the
	 * fermion determinant is not included in the pseudo-fermion action.
	 * Det(A_o) = -exp( Trlog(A_o^ A_o) ) */
    _MASTER_FOR(&glat_odd, ix) {
        double a = 0.;
        double ts = *_FIELD_AT(ff_sigma, ix);
        double tp = *_FIELD_AT(ff_pi, ix);
        double rho = 4. + par->mass + ts;
        int Nd = 4; //Number of dirac d.o.f.

        //Trace over site -> N_d*N_c. N_f == 2.
        a = -Nd * NF * log(rho * rho + tp * tp);
        *_FIELD_AT(loc_action, ix) += a;
    }
#endif
}

static void hmc_ff_free(monomial *m) {
    mon_hmc_par *par = (mon_hmc_par *)m->data.par;

    if (par->pf != NULL) { free_spinor_field(par->pf); }

    free(par);
    free(m);
}

monomial *hmc_ff_create(monomial_data const *data) {
    monomial *m = malloc(sizeof(*m));
    mon_hmc_par *par = (mon_hmc_par *)data->par;

    // Copy data structure
    m->data = *data;

    // Allocate memory for spinor field
    if (mon_init) {
        tmp_pf = alloc_spinor_field(1, &glat_default);
        mon_init = 0;
    }
    par->pf = alloc_spinor_field(1, &glat_default);

    // Setup force parameters
    par->fpar.id = data->id;
    par->fpar.n_pf = 1;
    par->fpar.pf = par->pf;
    par->fpar.inv_err2 = data->force_prec;
    par->fpar.inv_err2_flt = 1e-6;
    par->fpar.mass = par->mass;
    par->fpar.b = 0;
    par->fpar.hasenbusch = 0;
    par->fpar.mu = 0;
    par->fpar.momenta = &suN_momenta;

    // Setup chronological inverter
    mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

    // Setup pointers to update functions
    m->free = &hmc_ff_free;
    m->update_force = &force_hmc_ff;
    m->force_par = &par->fpar;
    m->update_field = 0;
    m->field_par = 0;

    m->pseudofermion = &hmc_ff_pseudofermion;
    m->gaussian_pf = &hmc_ff_gaussian_pf;
    m->correct_pf = &hmc_ff_correct_pf;
    m->correct_la_pf = &hmc_ff_correct_la_pf;
    m->add_local_action = &hmc_ff_add_local_action;

    return m;
}
