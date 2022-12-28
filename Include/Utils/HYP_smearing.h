#ifndef HYP_SMEARING_H
#define HYP_SMEARING_H

#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

/* HYP smearing */
//void spatialHYP_smearing(suNg_field *out, suNg_field *in, double weight[3]); //TODO: why are this is commented out
void HYP_smearing(suNg_field *out, suNg_field *in, double weight[3]);
double min_tplaq(suNg_field *g);
void HYP_span_parameters(double mtp[6859]);
int HYP_best_parameters(double mtp[6859], double w[3]);

#ifdef __cplusplus
	}
#endif
#endif //HYP_SMEARING_H
