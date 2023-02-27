#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "suN_types.h"
#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

// boundary_conditions.c
typedef struct
{
  double gauge_boundary_improvement_cs;
  double gauge_boundary_improvement_ct;
  double chiSF_boundary_improvement_ds;
  double fermion_twisting_theta[4];
  int SF_BCs;
  suNg gauge_boundary_up;
  suNg gauge_boundary_dn;
} BCs_pars_t;

void init_BCs(BCs_pars_t *pars);
void init_plaq_open_BCs(double *plaq_weight, double *rect_weight, double ct, double cs);
void free_BCs();
void apply_BCs_on_represented_gauge_field();
void apply_BCs_on_fundamental_gauge_field();
void apply_BCs_on_momentum_field(suNg_av_field *force);
void apply_BCs_on_spinor_field(spinor_field *sp);
void apply_BCs_on_spinor_field_flt(spinor_field_flt *sp);
void apply_background_field_zdir(suNg_field *V, double Q, int n);
void apply_BCs_on_clover_term(suNfc_field *);
#if defined(BC_T_SF) || defined(BC_T_SF_ROTATED)
void SF_classical_solution();
#endif


#ifdef __cplusplus
	}
#endif
#endif //BOUNDARY_CONDITIONS_H
