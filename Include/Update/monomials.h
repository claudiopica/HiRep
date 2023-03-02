/***************************************************************************\
* Copyright (c) 2017, Martin Hansen, Claudio Pica                           *
* All rights reserved.                                                      * 
\***************************************************************************/

/// Header file for:
/// - monomials_core.c
/// - mon_pg.c
/// - mon_lw.c
/// - mon_hmc.c
/// - mon_rhmc.c
/// - mon_tm.c
/// - mon_tm_alt.c
/// - mon_hasen.c
/// - mon_hasen_tm.c
/// - mon_hasen_tm_alt.c
/// - mon_ff.c
/// - mon_scalar.c

#ifndef MONOMIALS_H
#define MONOMIALS_H

#include "forces.h"
#include "rational_functions.h"

#ifdef __cplusplus
extern "C" {
#endif

//monomials_core.c

// List of monomials
typedef enum {
    PureGauge,
    LuscherWeisz,
    FourFermion,
    HMC,
    RHMC,
    TM,
    TM_alt,
    Hasenbusch,
    Hasenbusch_tm,
    Hasenbusch_tm_alt,
    HMC_ff,
    Hasenbusch_ff,
    Scalar
} mon_type;

typedef struct monomial_data {
    int id;
    mon_type type;
    void *par; //generic pointer to a monomial parameter structure of given type
    double MT_prec;
    double MD_prec;
    double force_prec;
} monomial_data;

typedef struct monomial {
    monomial_data data;
    void *force_par;
    void *field_par;
    void (*update_force)(double, void *);
    void (*update_field)(double, void *);
    void (*free)(struct monomial *);
    void (*gaussian_pf)(struct monomial const *);
    void (*correct_pf)(struct monomial const *);
    void (*correct_la_pf)(struct monomial const *);
    const spinor_field *(*pseudofermion)(struct monomial const *);
    void (*add_local_action)(struct monomial const *, scalar_field *);
} monomial;

monomial const *add_mon(monomial_data *mon_dat);
monomial const *mon_n(int i);
int num_mon(void);

//mon_pg.c
typedef struct mon_pg_par {
    double beta;
    force_gauge_par force_par;
    field_gauge_par field_par;
} mon_pg_par;

monomial *pg_create(monomial_data const *data);

//mon_lw.c
typedef struct mon_lw_par {
    double beta;
    double c0;
    double c1;
    force_gauge_par force_par;
    field_gauge_par field_par;
} mon_lw_par;

monomial *lw_create(monomial_data const *data);

//mon_hmc.c
typedef struct mon_hmc_par {
    double mass;
    int mre_past;
    force_hmc_par fpar;
    spinor_field *pf;
} mon_hmc_par;

monomial *hmc_create(monomial_data const *data);

//mon_rhmc.c
typedef struct mon_rhmc_par {
    double mass;
    rational_app ratio;
    force_rhmc_par fpar;
    spinor_field *pf;
} mon_rhmc_par;

monomial *rhmc_create(monomial_data const *data);

//mon_tm.c
typedef struct mon_tm_par {
    double mass;
    double mu;
    int mre_past;
    force_hmc_par fpar;
    spinor_field *pf;
} mon_tm_par;

monomial *tm_create(monomial_data const *data);

//mon_tm_alt.c
monomial *tm_alt_create(monomial_data const *data);

//mon_hasen.c
typedef struct mon_hasenbusch_par {
    double mass;
    double dm;
    int mre_past;
    force_hmc_par fpar;
    spinor_field *pf;
} mon_hasenbusch_par;

monomial *hasen_create(monomial_data const *data);

//mon_hasen_tm.c
typedef struct mon_hasenbusch_tm_par {
    double mass;
    double mu;
    double dmu;
    int mre_past;
    force_hmc_par fpar;
    spinor_field *pf;
} mon_hasenbusch_tm_par;

monomial *hasen_tm_create(monomial_data const *data);

//mon_hasen_tm_alt.c
monomial *hasen_tm_alt_create(monomial_data const *data);

//mon_ff.c
typedef struct mon_ff_par {
    double gamma;
    char *start_config;
    double start_value;
    force_auxfield_par fpar;
} mon_ff_par;

monomial *ff_create(monomial_data const *data);

//mon_hmc_ff.c
monomial *hmc_ff_create(monomial_data const *data);

//mon_hasen_ff.c
monomial *hasen_ff_create(monomial_data const *data);

//mon_scalar.c
typedef struct mon_scalar_par {
    double mass;
    double lambda;
    force_scalar_par force_par;
    field_scalar_par field_par;
} mon_scalar_par;

monomial *scalar_create(monomial_data const *data);

#ifdef __cplusplus
}
#endif
#endif
