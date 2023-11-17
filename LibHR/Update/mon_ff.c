/***************************************************************************\
 * Copyright (c) 2016, Jarno Rantaharju                                    *
 * All rights reserved.                                                    *
 \***************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "Inverters/scalarfield_operations.h"
#include <string.h>

// void write_ff_scalar_fields(char filename[]);
// void read_ff_scalar_fields(char filename[]);

static void ff_gaussian_pf(monomial const *m) {
    gaussian_scalar_field(ff_sigma_mom);
    gaussian_scalar_field(ff_pi_mom);
}

static void ff_correct_pf(monomial const *m) {
    /* empty */
}

static void ff_correct_la_pf(monomial const *m) {
    /* empty */
}

static const spinor_field *ff_pseudofermion(monomial const *m) {
    return NULL;
}

static void ff_add_local_action(monomial const *m, scalar_field *loc_action) {
    mon_ff_par *par = (mon_ff_par *)(m->data.par);
    _MASTER_FOR(&glattice, ix) {
        double a = 0.;
        double ts = *_FIELD_AT(ff_sigma, ix);
        double tp = *_FIELD_AT(ff_pi, ix);
        double g = par->gamma * par->gamma * 4.0;
        a = tp * tp / g;
        a += ts * ts / g;

        ts = *_FIELD_AT(ff_sigma_mom, ix);
        tp = *_FIELD_AT(ff_pi_mom, ix);
        a += 0.5 * (tp * tp + ts * ts);
        *_FIELD_AT(loc_action, ix) += a;
    }
}

static void ff_free(monomial *m) {
    mon_ff_par *par = (mon_ff_par *)m->data.par;
    free(ff_sigma);
    free(ff_pi);
    free(ff_sigma_mom);
    free(ff_pi_mom);
    free(par);
    free(m);
}

monomial *ff_create(monomial_data const *data) {
    monomial *m = malloc(sizeof(*m));
    mon_ff_par *par = (mon_ff_par *)(data->par);

    // Copy data structure
    m->data = *data;
    four_fermion_active = 1;

    // Allocate auxiliary field
    ff_sigma = alloc_scalar_field(1, &glattice);
    ff_pi = alloc_scalar_field(1, &glattice);
    ff_sigma_mom = alloc_scalar_field(1, &glattice);
    ff_pi_mom = alloc_scalar_field(1, &glattice);

    //Read or set auxiliary field configuration or set cold value
    if (strcmp(par->start_config, "cold") == 0) {
        set_scalar_field(ff_sigma, par->start_value);
        set_scalar_field(ff_pi, 0);
    }

    // Setup force parameters
    par->fpar.gamma = par->gamma;

    // Setup pointers to update functions
    m->free = &ff_free;
    m->update_force = &force_hmc_auxfields;
    m->force_par = &par->fpar;
    m->update_field = &update_auxfields;
    m->field_par = 0;

    m->pseudofermion = &ff_pseudofermion;
    m->gaussian_pf = &ff_gaussian_pf;
    m->correct_pf = &ff_correct_pf;
    m->correct_la_pf = &ff_correct_la_pf;
    m->add_local_action = &ff_add_local_action;

    return m;
}
