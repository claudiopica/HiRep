/// Headerfile for:
/// - suN_exp_group.c
/// - suN_utils.c

#ifndef SUN_MAT_UTILS_H
#define SUN_MAT_UTILS_H

#include "libhr_core.h"

#ifdef __cplusplus
extern "C" {
#endif

// suN_exp_group.c
// SUN exp matrix
visible void suNg_Exp(suNg *u, suNg *Xin); //global function pointer to the correct implementation
visible void ExpX(double dt, suNg_algebra_vector *h, suNg *u);
visible void suNg_Exp_Taylor(suNg *u, suNg *Xin);

//suN_utils.c
visible void vector_star(suNg_vector *, suNg_vector *);
visible void project_to_suNg(suNg *u);
visible void project_to_suNg_flt(suNg_flt *u);
#ifndef GAUGE_SON
visible void project_cooling_to_suNg(suNg *g_out, suNg *g_in, int cooling);
#endif
visible void covariant_project_to_suNg(suNg *u);
#ifdef GAUGE_SON
int project_to_suNg_real(suNg *out, suNg *in);
#endif

#ifdef __cplusplus
}
#endif

#endif //SUN_MAT_UTILS_H
