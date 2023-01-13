#ifndef MEASURE_BARYONS_H
#define MEASURE_BARYONS_H
#if NG==3

#include "spinor_field.h"
#include "Utils/data_storage.h"

#ifdef __cplusplus
	extern "C" {
#endif

void contract_baryons(spinor_field *psi0, int tau, storage_switch swc, data_storage_array **ret);

#endif

#ifdef __cplusplus
	}
#endif
#endif //MEASURE_BARYONS_H
