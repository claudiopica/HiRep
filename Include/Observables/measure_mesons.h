// Header file for:
// - measure_meson.c
// - measure_ff.c

#ifndef MEASURE_MESONS_H
#define MEASURE_MESONS_H

#include "spinor_field.h"
#include "meson_observables.h"

#ifdef __cplusplus
	extern "C" {
#endif

//measure_meson.c
extern meson_observable *meson_correlators;
extern meson_observable *discon_correlators;
extern meson_observable *cvc_correlators;

void init_meson_correlators(int meas_offdiag);
void init_discon_correlators(void);
void init_vcvl_correlators(void); //TODO: not defined in lib
void init_cvc_correlators(void);
void free_meson_observables(void);

void measure_mesons_core(spinor_field *psi0, spinor_field *psi1, spinor_field *eta, meson_observable *mo, int nm, int tau, int n_mom, int offset, int lt);

void measure_mesons(meson_observable *mo, spinor_field *psi0, spinor_field *eta, int nm, int tau);
void measure_diquarks(meson_observable *mo, spinor_field *psi0, spinor_field *psi1, spinor_field *eta, int nm, int tau);
void measure_conserved_currents(meson_observable *mo, spinor_field *psi0, spinor_field *eta, int nm, int tau);
void measure_mesons_ext(meson_observable *mo, spinor_field *psi0, spinor_field *eta, int nm, int tau, int begin);
void measure_point_mesons_momenta(meson_observable *mo, spinor_field *psi0, spinor_field *eta, int nm, int tau, int n_mom);
void measure_point_mesons_momenta_ext(meson_observable *mo, spinor_field *psi0, spinor_field *eta, int nm, int tau, int n_mom, int begin);
void print_mesons(meson_observable *mo, double norm, int conf, int nm, double *mass, int lt, int n_mom, char *label);

//measure_ff.c
void measure_formfactors(spinor_field *psi0, spinor_field *psi1, spinor_field *eta, int nm, int ti, int tf, int n_mom, int *pt);
void measure_formfactors_ext(spinor_field *psi0, spinor_field *psi1, spinor_field *eta, int nm, int ti, int tf, int n_mom, int begin); //TODO: not defined in lib
void print_formfactor(int conf, int nm, double *mass, int n_mom, char *label, int tf);
void print_formfactor_ext(int conf, int nm, double *mass, int n_mom, char *label, int tf);

#ifdef __cplusplus
	}
#endif
#endif //MEASURE_MESONS_H
