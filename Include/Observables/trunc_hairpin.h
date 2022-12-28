#ifndef TRUNC_HAIRPIN_H
#define TRUNC_HAIRPIN_H

#include "hr_complex.h"

#ifdef __cplusplus
	extern "C" {
#endif

typedef enum
{
  NO_DILUTION,
  TIME_DILUTION,
  TIME_SPIN_DILUTION,
  EXACT
} dilution_mode;

typedef struct _ata_qprop_pars
{
  int n_masses;
  double mass[256];
  int n_eigenvalues;
  int eva_nevt;
  double eva_omega1;
  double eva_omega2;
  int eva_imax;
  int eva_kmax;
  int hopping_order;
  int n_truncation_steps;
  int n_sources_truncation;
  int n_sources_correction;
  dilution_mode dilution;
  double inverter_precision;
} ata_qprop_pars;

void traced_ata_qprop(hr_complex ***prop, int n_points);
void ata_qprop_init(ata_qprop_pars *p);
void ata_qprop_free(void);

#ifdef __cplusplus
	}
#endif
#endif //TRUNC_HAIRPIN_H
