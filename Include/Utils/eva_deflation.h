#ifndef EVA_DEFLATION_H
#define EVA_DEFLATION_H

#include "spinor_field.h"
#include "Inverters/linear_solvers.h"
#include "Geometry/geometry_descriptor.h"

#ifdef __cplusplus
extern "C" {
#endif

//eva_deflation.c
/* EVA preconditioning */
typedef struct eva_prec {
    /* EVA parameters */
    int nevt; /* search space dimension */
    int nev; /* number of accurate eigenvalues */
    int kmax; /* max degree of polynomial */
    int maxiter; /* max number of subiterations */
    double omega1; /* absolute precision */
    double omega2; /* relative precision */
} eva_prec;
void set_def_matrix(eva_prec *e_par, spinor_operator H, geometry_descriptor *type);
void eva_def(spinor_field *out, spinor_field *in);
void eva_def_inv(spinor_field *out, spinor_field *in, double m);

#ifdef __cplusplus
}
#endif
#endif //EVA_DEFLATION_H
