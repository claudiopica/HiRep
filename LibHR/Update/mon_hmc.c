/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica                        *
 * All rights reserved.                                                   *
 \***************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "random.h"
#include "Inverters/linear_algebra.h"

static spinor_field *tmp_pf = NULL;
static int mon_init = 1;

static void hmc_gaussian_pf(monomial const *m) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);
    gaussian_spinor_field(par->pf);
}

static void hmc_correct_pf(monomial const *m) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);

    /* compute H2^{1/2}*pf = H*pf */
    copy_spinor_field(tmp_pf, par->pf);
    set_dirac_mass(par->mass);
    H(par->pf, tmp_pf);
}

static void hmc_correct_la_pf(monomial const *m) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);
    double shift;

    mshift_par mpar;
    mpar.err2 = m->data.MT_prec;
    mpar.max_iter = 0;
    mpar.n = 1;
    mpar.shift = &shift;
    mpar.shift[0] = 0;

    /* compute H2^{-1/2}*pf = H^{-1}*pf */
    g5_spinor_field(tmp_pf, par->pf);
    set_dirac_mass(par->mass);
    zero_spinor_field(par->pf); /* mshift inverter uses this as initial guess for 1 shift */
    g5QMR_mshift(&mpar, &D, tmp_pf, par->pf);
}

static const spinor_field *hmc_pseudofermion(monomial const *m) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);
    return par->pf;
}

static void hmc_add_local_action(monomial const *m, scalar_field *loc_action) {
    mon_hmc_par *par = (mon_hmc_par *)(m->data.par);
    pf_local_action(loc_action, par->pf);
#ifdef WITH_CLOVER_EO
    clover_la_logdet(2., par->mass, loc_action);
#endif
}

static void hmc_free(monomial *m) {
    mon_hmc_par *par = (mon_hmc_par *)m->data.par;

    if (par->pf != NULL) { free_spinor_field(par->pf); }

    free(par);
    free(m);
}

monomial *hmc_create(monomial_data const *data) {
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
    par->fpar.logdet = 1;
    par->fpar.momenta = &suN_momenta;

    // Setup chronological inverter
    mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

    // Setup pointers to update functions
    m->free = &hmc_free;
    m->update_force = &force_hmc;
    m->force_par = &par->fpar;
    m->update_field = 0;
    m->field_par = 0;

    m->pseudofermion = &hmc_pseudofermion;
    m->gaussian_pf = &hmc_gaussian_pf;
    m->correct_pf = &hmc_correct_pf;
    m->correct_la_pf = &hmc_correct_la_pf;
    m->add_local_action = &hmc_add_local_action;

    return m;
}
