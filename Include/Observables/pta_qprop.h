#ifndef PTA_QPROP_H
#define PTA_QPROP_H

#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

void pta_qprop_QMR_eo(int g0[4], spinor_field **pta_qprop, int nm, double *m, double acc);
void pta_qprop_QMR(int g0[4], spinor_field **pta_qprop, int nm, double *m, double acc);
void pta_qprop_MINRES(int g0[4], spinor_field **pta_qprop, int nm, double *m, double acc);


#ifdef __cplusplus
	}
#endif
#endif //PTA_QPROP_H
