#ifndef CALC_PROP_H
#define CALC_PROP_H

#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

//calc_prop.h
void init_propagator_eo(int nm, double *m, double acc);
void free_propagator_eo(void);
void eig_init(int nev, int nevt, int kmax, int maxiter, double lbnd, double omega1, double omega2);
void calc_propagator(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator_eo(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator_tw(double *m, double mu, spinor_field *psi, spinor_field *eta, int ndilute);
void calc_propagator_multisource(spinor_field *psi, spinor_field *eta, int ndilute);
void calc_deflated_propagator(spinor_field *psi, spinor_field *eta, int ndilute, int Nuse);
void copy_evec(int n, spinor_field *psi1, double *eval);

#ifdef __cplusplus
	}
#endif
#endif //CALC_PROP_H
