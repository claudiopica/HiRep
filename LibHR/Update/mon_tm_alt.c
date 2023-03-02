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

void tm_alt_gaussian_pf(monomial const *m) {
    mon_tm_par *par = (mon_tm_par *)(m->data.par);
    gaussian_spinor_field(par->pf);
}

void tm_alt_correct_pf(monomial const *m) {
    mon_tm_par *par = (mon_tm_par *)(m->data.par);

    /* psi =  (g5 D+imu) pf */
    spinor_field_copy_f(tmp_pf, par->pf);
    set_dirac_mass(par->mass);
    set_twisted_mass(par->mu);
    Qtm_p_alt(par->pf, tmp_pf);
}

void tm_alt_correct_la_pf(monomial const *m) {
    mon_tm_par *par = (mon_tm_par *)(m->data.par);
    double shift;

    mshift_par mpar;
    mpar.err2 = m->data.MT_prec;
    mpar.max_iter = 0;
    mpar.n = 1;
    mpar.shift = &shift;
    mpar.shift[0] = 0;

    /* compute H2^{-1/2}*pf = H^{-1}*pf */
    set_dirac_mass(par->mass);
    set_twisted_mass(par->mu);
    spinor_field_copy_f(tmp_pf, par->pf);
    spinor_field_zero_f(par->pf);
    tm_invert_alt(par->pf, tmp_pf, &mpar);
}

const spinor_field *tm_alt_pseudofermion(monomial const *m) {
    mon_tm_par *par = (mon_tm_par *)(m->data.par);
    return par->pf;
}

void tm_alt_add_local_action(monomial const *m, scalar_field *loc_action) {
    mon_tm_par *par = (mon_tm_par *)(m->data.par);
    pf_local_action(loc_action, par->pf);

#ifdef WITH_CLOVER_EO
    clover_la_logdet(2., par->mass, loc_action);
#endif
}

void tm_alt_free(monomial *m) {
    mon_tm_par *par = (mon_tm_par *)m->data.par;

    if (par->pf != NULL) { free_spinor_field_f(par->pf); }

    free(par);
    free(m);
}

monomial *tm_alt_create(monomial_data const *data) {
    monomial *m = malloc(sizeof(*m));
    mon_tm_par *par = (mon_tm_par *)data->par;

    // Copy data structure
    m->data = *data;

    // Allocate memory for spinor field
    if (mon_init) {
        tmp_pf = alloc_spinor_field_f(1, &glat_default);
        mon_init = 0;
    }
    par->pf = alloc_spinor_field_f(1, &glat_default);

    // Setup force parameters
    par->fpar.id = data->id;
    par->fpar.n_pf = 1;
    par->fpar.pf = par->pf;
    par->fpar.inv_err2 = data->force_prec;
    par->fpar.inv_err2_flt = 1e-6;
    par->fpar.mass = par->mass;
    par->fpar.mu = par->mu;
    par->fpar.b = 0;
    par->fpar.hasenbusch = 0;
    par->fpar.logdet = 1;
    par->fpar.momenta = &suN_momenta;

    // Setup chronological inverter
    mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

    // Setup pointers to update functions
    m->free = &tm_alt_free;
    m->update_force = &force_hmc;
    m->force_par = &par->fpar;
    m->update_field = 0;
    m->field_par = 0;

    m->pseudofermion = &tm_alt_pseudofermion;
    m->gaussian_pf = &tm_alt_gaussian_pf;
    m->correct_pf = &tm_alt_correct_pf;
    m->correct_la_pf = &tm_alt_correct_la_pf;
    m->add_local_action = &tm_alt_add_local_action;

    return m;
}
