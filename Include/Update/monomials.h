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

typedef struct {
	int id;
	mon_type type;
	void *par; //generic pointer to a monomial parameter structure of given type
	double MT_prec;
	double MD_prec;
	double force_prec;
} monomial_data;

typedef struct _monomial {
	monomial_data data;
	void *force_par;
	void *field_par;
	void (*update_force)(double, void*);
	void (*update_field)(double, void*);
	void (*free)(struct _monomial*);
	void (*gaussian_pf)(const struct _monomial*);
	void (*correct_pf)(const struct _monomial*);
	void (*correct_la_pf)(const struct _monomial*);
	const spinor_field *(*pseudofermion)(const struct _monomial*);
	void (*add_local_action)(const struct _monomial*, scalar_field*);
} monomial;

const monomial *add_mon(monomial_data *mon_dat);
const monomial *mon_n(int i);
int num_mon(void);

//mon_pg.c
typedef struct {
	double beta;
	force_gauge_par force_par;
	field_gauge_par field_par;
} mon_pg_par;

monomial *pg_create(const monomial_data *data);

//mon_lw.c
typedef struct {
	double beta;
	double c0;
	double c1;
	force_gauge_par force_par;
	field_gauge_par field_par;
} mon_lw_par;

monomial *lw_create(const monomial_data *data);

//mon_hmc.c
typedef struct {
	double mass;
	int mre_past;
	force_hmc_par fpar;
	spinor_field *pf;
} mon_hmc_par;

monomial *hmc_create(const monomial_data *data);

//mon_rhmc.c
typedef struct {
	double mass;
	rational_app ratio;
	force_rhmc_par fpar;
	spinor_field *pf;
} mon_rhmc_par;

monomial *rhmc_create(const monomial_data *data);

//mon_tm.c
typedef struct {
	double mass;
	double mu;
	int mre_past;
	force_hmc_par fpar;
	spinor_field *pf;
} mon_tm_par;

monomial *tm_create(const monomial_data *data);

//mon_tm_alt.c
monomial *tm_alt_create(const monomial_data *data);

//mon_hasen.c
typedef struct {
	double mass;
	double dm;
	int mre_past;
	force_hmc_par fpar;
	spinor_field *pf;
} mon_hasenbusch_par;

monomial *hasen_create(const monomial_data *data);

//mon_hasen_tm.c 
typedef struct {
	double mass;
	double mu;
	double dmu;
	int mre_past;
	force_hmc_par fpar;
	spinor_field *pf;
} mon_hasenbusch_tm_par;

monomial *hasen_tm_create(const monomial_data *data);

//mon_hasen_tm_alt.c
monomial *hasen_tm_alt_create(const monomial_data *data);

//mon_ff.c
typedef struct {
	double gamma;
	char *start_config;
	double start_value;
	force_auxfield_par fpar;
} mon_ff_par;

monomial *ff_create(const monomial_data *data);

//mon_hmc_ff.c
monomial *hmc_ff_create(const monomial_data *data);

//mon_hasen_ff.c
monomial *hasen_ff_create(const monomial_data *data);

//mon_scalar.c
typedef struct {
	double mass;
	double lambda;
	force_scalar_par force_par;
	field_scalar_par field_par;
} mon_scalar_par;

monomial *scalar_create(const monomial_data *data);


#ifdef __cplusplus
	}
#endif
#endif
