
/*******************************************************************************
*
* File observables.h
*
* Functions for measuring observables
*
*******************************************************************************/

#ifndef LOOP_TOOLS_H
#define LOOP_TOOLS_H

#if !defined(ROTATED_SF) && !defined(BASIC_SF) && !defined(FERMION_THETA)

#include <stdio.h>
#include "spinor_field.h"
#include "Utils/data_storage.h"

#ifdef __cplusplus
	extern "C" {
#endif

void measure_bilinear_loops_spinorfield(spinor_field* prop,spinor_field* source,int src_id, int n_mom, storage_switch swc, data_storage_array **ret);
void measure_bilinear_loops_4spinorfield(spinor_field* prop,spinor_field* source,int src_id,int tau,int col,int eo,storage_switch swc, data_storage_array **ret);
void measure_loops(double* m, int nhits,int conf_num, double precision,int source_type,int n_mom,storage_switch swc, data_storage_array **ret);

#ifdef __cplusplus
    }
#endif

#endif
#endif
