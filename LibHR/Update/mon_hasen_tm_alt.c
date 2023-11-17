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

static void hasen_tm_alt_gaussian_pf(monomial const *m) {
    mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par *)(m->data.par);
    gaussian_spinor_field(par->pf);
}

static void hasen_tm_alt_correct_pf(monomial const *m) {
    mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par *)(m->data.par);
    double shift;

    mshift_par mpar;
    mpar.err2 = m->data.MT_prec;
    mpar.max_iter = 0;
    mpar.n = 1;
    mpar.shift = &shift;
    mpar.shift[0] = 0;

    /* S = | (g5 D + i mu ) D^{-1} g5 psi |^2 */
    /* ( g5 D + i mu ) D^{-1} g5 psi = A */
    /* psi = g5 D ( g5 D + i mu )^{-1}  A */
    set_dirac_mass(par->mass);
    set_twisted_mass(par->mu + par->dmu);
    spinor_field_zero_f(tmp_pf);
    tm_invert_alt(tmp_pf, par->pf, &mpar);
    /*set_twisted_mass(par->mu);*/
    set_twisted_mass(par->mu);
    Qtm_p_alt(par->pf, tmp_pf);
}

static void hasen_tm_alt_correct_la_pf(monomial const *m) {
    mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par *)(m->data.par);
    double shift;

    mshift_par mpar;
    mpar.err2 = m->data.MT_prec;
    mpar.max_iter = 0;
    mpar.n = 1;
    mpar.shift = &shift;
    mpar.shift[0] = 0;

    /* S = | (g5 D + i mu ) D^{-1} g5 psi |^2 */
    spinor_field_copy_f(tmp_pf, par->pf);
    spinor_field_zero_f(par->pf);
    set_dirac_mass(par->mass);
    set_twisted_mass(par->mu);
    tm_invert_alt(par->pf, tmp_pf, &mpar);
    set_twisted_mass(par->mu + par->dmu);
    Qtm_p_alt(tmp_pf, par->pf);
    spinor_field_copy_f(par->pf, tmp_pf);
}

static const spinor_field *hasen_tm_alt_pseudofermion(monomial const *m) {
    mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par *)(m->data.par);
    return par->pf;
}

static void hasen_tm_alt_add_local_action(monomial const *m, scalar_field *loc_action) {
    mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par *)(m->data.par);
    pf_local_action(loc_action, par->pf);
}

static void hasen_tm_alt_free(monomial *m) {
    mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par *)m->data.par;

    if (par->pf != NULL) { free_spinor_field(par->pf); }

    free(par);
    free(m);
}

monomial *hasen_tm_alt_create(monomial_data const *data) {
    monomial *m = malloc(sizeof(*m));
    mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par *)data->par;

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
    par->fpar.mu = par->mu;
    par->fpar.b = par->dmu;
    par->fpar.hasenbusch = 2;
    par->fpar.logdet = 0;
    par->fpar.momenta = &suN_momenta;

    // Setup chronological inverter
    mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

    // Setup pointers to update functions
    m->free = &hasen_tm_alt_free;
    m->update_force = &force_hmc;
    m->force_par = &par->fpar;
    m->update_field = 0;
    m->field_par = 0;

    m->pseudofermion = &hasen_tm_alt_pseudofermion;
    m->gaussian_pf = &hasen_tm_alt_gaussian_pf;
    m->correct_pf = &hasen_tm_alt_correct_pf;
    m->correct_la_pf = &hasen_tm_alt_correct_la_pf;
    m->add_local_action = &hasen_tm_alt_add_local_action;

    return m;
}
