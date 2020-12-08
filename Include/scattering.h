
/*******************************************************************************
*
* File scattering.h
* 
* Functions for Finite Size Method 
*
*******************************************************************************/

#ifndef SCATTERING_H
#define SCATTERING_H

#include "suN.h"
#include "meson_observables.h"
#include <stdio.h>

void measure_pion_scattering(double* m, int nhits,int conf_num, double precision,int ts);
void measure_pion_scattering_I2(double* m, int numsources, double precision,char* path,char* cnfg_filename,meson_observable** mo_arr);
void measure_pion_scattering_I0(double* m, int numsources, double precision,char* path,char* cnfg_filename, int seq_prop, meson_observable** mo_arr);
void measure_pion_scattering_I0_TS(double* m, int numsources, double precision,char* path,char* cnfg_filename, int seq_prop, meson_observable** mo_arr);

/**
 * @brief Propagator sources with zero-momentum
 */
struct src_common{
	spinor_field *src_0; /**< Zero momentum, noise 1 */
	spinor_field *src_0_eta; /**< Zero momentum, noise 2 */
	spinor_field *src_0_0; /**< Sequential source from timeslice tau with zero-momentum insertion */
};

/**
 * @brief Propagator sources with momentum p
 */
struct src_p{
	int p[3]; /**< Momentum */
	spinor_field *src_p;
	spinor_field *src_mp;
	spinor_field *src_0_p;
	spinor_field *src_0_mp;
	spinor_field *src_p_0;
	spinor_field *src_mp_0;
};


/**
 * @brief Bundle of propagators with zero momentum.
 */
struct prop_common{
	spinor_field *Q_0;
	spinor_field *Q_0_eta;
	spinor_field **W_0_0;
};

/**
 * @brief Bundle of propagators with momentum p.
 */
struct prop_p{
	spinor_field *Q_p;
	spinor_field *Q_mp;
	spinor_field *W_0_p;
	spinor_field *W_0_mp;
	spinor_field *W_p_0;
	spinor_field *W_mp_0;
};


/**
 * @brief Bundle of meson_observables with momentum 0. 
 */
struct mo_0{
	meson_observable *rho[3][3];
	meson_observable *pi;
};

/**
 * @brief Bundle of meson_observables with momentum p. 
 */
struct mo_p{
	int p[3];
	meson_observable *d, *r1, *r2, *r3, *r4, *pi;
	meson_observable *t1[3], *t2[3], *rho[3][3];
};






int** getmomlist(char* momstring, int* N);
void freep(int **p, int N);
void init_mo(meson_observable* mo, char* name, int size);
void reset_mo(meson_observable* mo);
void free_mo(meson_observable* mo);
void init_src_common(struct src_common* src, int tau);
void init_src_common_point(struct src_common* src, int tau);
void init_src_p(struct src_p* srcp, struct src_common* src0, int px, int py, int pz);
void free_src_common(struct src_common* src);
void free_src_p(struct src_p* src);
void make_propagator_P(spinor_field* prop, spinor_field* src, int ndilute, int tau);
void make_propagator_PA(spinor_field* prop, spinor_field* src, int ndilute, int tau);
void make_prop_common(struct prop_common* prop, struct src_common* src0, int ndilute, int tau, char* bc);
void free_mo_0(struct mo_0* mo);
void free_mo_p(struct mo_p* mo);
void io2pt(meson_observable* mo, int pmax, int sourceno, char* path, char* name,char * cnfg_filename);
void io4pt(meson_observable* mo, int pmax, int sourceno, char* path, char* name,char * cnfg_filename);
void IOold_0(struct mo_0* molist[], int numsources, char* path, char* cnfg_filename);
void IOold_p(struct mo_p* molist[], int numsources, char* path, char* cnfg_filename	);
void IO_json_0(struct mo_0* molist[], int numsources, char* path,char * cnfg_filename);
void IO_json_p(struct mo_p* molist[], int numsources, char* path, char* cnfg_filename);
void init_mo_0(struct mo_0* mo);
void init_mo_p(struct mo_p* mo, int px, int py, int pz);
void gen_mo_0(struct mo_0* mo, struct prop_common* p0, struct src_common* s0, int tau);
void gen_mo_p(struct mo_p* mo, struct prop_common* p0, struct prop_p* pp, struct src_common* s0, int tau);
void make_prop_p(struct prop_p* prop, struct src_p* srcp, struct src_common* src0, int ndilute, int tau, char* bc);
void free_prop_common(struct prop_common* prop);
void free_prop_p(struct prop_p* prop);


#endif
