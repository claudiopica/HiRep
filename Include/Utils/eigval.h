#ifndef EIGVAL_H
#define EIGVAL_H

#include "Geometry/geometry_descriptor.h"
#include "Inverters/linear_solvers.h"

#ifdef __cplusplus
extern "C" {
#endif

/* use inverse power method to find smallest eigenvalue of H */
int min_eigval(spinor_operator H, geometry_descriptor *type, double *min);

/* use power method to find largest eigenvalue of H */
int max_eigval(spinor_operator H, geometry_descriptor *type, double *max);

#ifdef __cplusplus
}
#endif
#endif
