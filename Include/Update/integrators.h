#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "monomials.h"

#ifdef __cplusplus
extern "C" {
#endif

//integrators.c
typedef struct integrator_par {
    int nsteps;
    int nmon;
    monomial const **mon_list;
    void (*integrator)(double, struct integrator_par *);
    struct integrator_par *next;
    int level;
} integrator_par;

void leapfrog_multistep(double tlen, integrator_par *int_par);
void O2MN_multistep(double tlen, integrator_par *int_par);
void O4MN_multistep(double tlen, integrator_par *int_par);

#ifdef __cplusplus
}
#endif
#endif //INTEGRATORS_H
