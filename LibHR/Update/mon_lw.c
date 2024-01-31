/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica                        *
 * All rights reserved.                                                   *
 \***************************************************************************/

#include "update.h"
#include "libhr_core.h"

static void lw_gaussian_pf(monomial const *m) {
    /* empty */
}

static void lw_correct_pf(monomial const *m) {
    /* empty */
}

static void lw_correct_la_pf(monomial const *m) {
    /* empty */
}

static const spinor_field *lw_pseudofermion(monomial const *m) {
    return NULL;
}

static void lw_add_local_action(monomial const *m, scalar_field *loc_action) {
    mon_lw_par *par = (mon_lw_par *)(m->data.par);
    lw_local_action(loc_action, par->beta, par->c0, par->c1);
}

static void lw_free(monomial *m) {
    mon_lw_par *par = (mon_lw_par *)m->data.par;
    free(par);
    free(m);
}

monomial *lw_create(monomial_data const *data) {
    monomial *m = malloc(sizeof(*m));
    mon_lw_par *par = (mon_lw_par *)(data->par);

    // Copy data structure
    m->data = *data;
    par->c1 = (1. - par->c0) / 8.;
    gauge_field_active = 1;

    // Setup force parameters
    par->force_par.beta = par->beta;
    par->force_par.c0 = par->c0;
    par->force_par.c1 = par->c1;
    par->force_par.momenta = &suN_momenta;
    par->field_par.field = &u_gauge;
    par->field_par.momenta = &suN_momenta;

    // Setup pointers to update functions
    m->free = &lw_free;
    m->update_force = lw_force;
    m->force_par = &par->force_par;
    m->update_field = &update_gauge_field;
    m->field_par = &par->field_par;

    m->pseudofermion = &lw_pseudofermion;
    m->gaussian_pf = &lw_gaussian_pf;
    m->correct_pf = &lw_correct_pf;
    m->correct_la_pf = &lw_correct_la_pf;
    m->add_local_action = &lw_add_local_action;

    return m;
}
