#ifndef MAX_EIG_H
#define MAX_EIG_H

#include "Geometry/geometry_descriptor.h"
#include "Inverters/linear_solvers.h"

#ifdef __cplusplus
	extern "C" {
#endif

/* use power method to find max eigvalue of H2 */
int max_H(spinor_operator H, geometry_descriptor *type, double *max);

#ifdef __cplusplus
	}
#endif
#endif //MAX_EIG_H
