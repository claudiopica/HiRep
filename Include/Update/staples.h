#ifndef STAPLES_H
#define STAPLES_H

#include "suN_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void staples(int ix, int mu, suNg *v);
// void test_staples(); //TODO: it is commented out

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
__device__ void staples_dev(int ix, int mu, suNg *v, suNg *gauge, int *iup_gpu, int *idn_gpu, double *plaq_weight);
#endif
#endif //STAPLES_H
