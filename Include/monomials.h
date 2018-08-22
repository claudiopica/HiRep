/***************************************************************************\
* Copyright (c) 2017, Martin Hansen, Claudio Pica                           *
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef MONOMIALS_H
#define MONOMIALS_H

//
// This header file is included in update.h
// It should not be included directly!
//

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
	double beta;
	force_gauge_par force_par;
	field_gauge_par field_par;
} mon_pg_par;

typedef struct {
	double beta;
	double c0;
	double c1;
	force_gauge_par force_par;
	field_gauge_par field_par;
} mon_lw_par;

typedef struct {
	double gamma;
	char *start_config;
	double start_value;
	force_auxfield_par fpar;
} mon_ff_par;

typedef struct {
	double mass;
	int mre_past;
	force_hmc_par fpar;
	spinor_field *pf;
} mon_hmc_par;

typedef struct {
	double mass;
	rational_app ratio;
	force_rhmc_par fpar;
	spinor_field *pf;
} mon_rhmc_par;

typedef struct {
	double mass;
	double mu;
	int mre_past;
	force_hmc_par fpar;
	spinor_field *pf;
} mon_tm_par;

typedef struct {
	double mass;
	double dm;
	int mre_past;
	force_hmc_par fpar;
	spinor_field *pf;
} mon_hasenbusch_par;

typedef struct {
	double mass;
	double mu;
	double dmu;
	int mre_past;
	force_hmc_par fpar;
	spinor_field *pf;
} mon_hasenbusch_tm_par;

typedef struct {
	double mass;
	double lambda;
	force_scalar_par force_par;
	field_scalar_par field_par;
} mon_scalar_par;

typedef struct {
	int id;
	mon_type type;
	void *par;
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

struct _monomial* pg_create(const monomial_data*);
struct _monomial* lw_create(const monomial_data*);
struct _monomial* hmc_create(const monomial_data*);
struct _monomial* rhmc_create(const monomial_data*);
struct _monomial* tm_create(const monomial_data*);
struct _monomial* tm_alt_create(const monomial_data*);
struct _monomial* hasen_create(const monomial_data*);
struct _monomial* hasen_tm_create(const monomial_data*);
struct _monomial* hasen_tm_alt_create(const monomial_data*);
struct _monomial* ff_create(const monomial_data*);
struct _monomial* hmc_ff_create(const monomial_data*);
struct _monomial* hasen_ff_create(const monomial_data*);
struct _monomial* scalar_create(const monomial_data*);

const monomial *add_mon(monomial_data*);
const monomial *mon_n(int);
int num_mon();

#endif
