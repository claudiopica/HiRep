#ifndef WILSONLOOPS_H
#define WILSONLOOPS_H

#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

//wilsonloops.c
#define _WL_3VOL_INDEX(x, y, z) ((x) + (y)*X + (z)*X * Y)

void WL_initialize();
void WL_free();
void WL_load_path(int c[3], int nsteps);
void WL_Hamiltonian_gauge(suNg_field *out, suNg_field *in);
void WL_broadcast_polyakov(suNg *poly, suNg_field *gf);
void WL_correlators(double **ret, const suNg_field *gf, const suNg *poly, const int nsteps, const int *path, const int length, const int perm[3], int sign[3]);
void WL_wilsonloops(double HYP_weight[3]);

typedef struct WL_path_t {
  int c[3];
  int* path;
  int length;
  int** perm;
  int nperms;
  int nsteps;
} WL_path_t;
//defined in wilsonloops.c also used in tests
extern WL_path_t WL_path[256];
extern int WL_npaths;
extern int WL_max_nsteps;


#ifdef __cplusplus
	}
#endif
#endif //WILSONLOOPS_H
