/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica                        *
 * All rights reserved.                                                   *
 \***************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "Update/avr_plaquette.h"

static void pg_gaussian_pf(monomial const *m) {
    /* empty */
}

static void pg_correct_pf(monomial const *m) {
    /* empty */
}

static void pg_correct_la_pf(monomial const *m) {
    /* empty */
}

static const spinor_field *pg_pseudofermion(monomial const *m) {
    return NULL;
}

static void pg_add_local_action(monomial const *m, scalar_field *loc_action) {
    mon_pg_par *par = (mon_pg_par *)(m->data.par);
    _MASTER_FOR(&glattice, ix) {
        *_FIELD_AT(loc_action, ix) += -(par->beta / ((double)NG)) * local_plaq(ix);
    }
}

static void pg_free(monomial *m) {
    mon_pg_par *par = (mon_pg_par *)m->data.par;
    free(par);
    free(m);
}

monomial *pg_create(monomial_data const *data) {
    monomial *m = malloc(sizeof(*m));
    mon_pg_par *par = (mon_pg_par *)(data->par);

    // Copy data structure
    m->data = *data;
    gauge_field_active = 1;

    // Setup force parameters
    par->force_par.beta = par->beta;
    par->force_par.momenta = &suN_momenta;
    par->field_par.field = &u_gauge;
    par->field_par.momenta = &suN_momenta;

    // Setup pointers to update functions
    m->free = &pg_free;
    m->update_force = &force0;
    m->force_par = &par->force_par;
    m->update_field = &update_gauge_field;
    m->field_par = &par->field_par;

    m->pseudofermion = &pg_pseudofermion;
    m->gaussian_pf = &pg_gaussian_pf;
    m->correct_pf = &pg_correct_pf;
    m->correct_la_pf = &pg_correct_la_pf;
    m->add_local_action = &pg_add_local_action;

    return m;
}
