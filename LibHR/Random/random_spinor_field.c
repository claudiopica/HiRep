/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "random.h"
#include "libhr_core.h"
#include "Inverters/linear_algebra.h"
#include "Utils/boundary_conditions.h"

void gaussian_spinor_field(spinor_field *s)
{
	const double c1 = 1. / sqrt(2.);
	geometry_descriptor *type = s->type;
	for (int i = 0; i < type->local_master_pieces; i++)
		gauss((double *)(s->ptr + (type->master_start[i] - type->master_shift)), (type->master_end[i] - type->master_start[i] + 1) * sizeof(suNf_spinor) / sizeof(double));

	spinor_field_mul_f_cpu(s, c1, s);
	apply_BCs_on_spinor_field(s);
}

void gaussian_spinor_field_flt(spinor_field_flt *s)
{
	const float c1 = (float)(1. / sqrt(2.));
	geometry_descriptor *type = s->type;
	for (int i = 0; i < type->local_master_pieces; i++)
		gauss_flt((float *)(s->ptr + (type->master_start[i] - type->master_shift)), (type->master_end[i] - type->master_start[i] + 1) * sizeof(suNf_spinor_flt) / sizeof(float));

	spinor_field_mul_f_flt_cpu(s, c1, s);
	apply_BCs_on_spinor_field_flt(s);
}

void z2_spinor_field(spinor_field *s)
{
	geometry_descriptor *type = s->type;
	for (int i = 0; i < type->local_master_pieces; i++)
		ranz2((double *)(s->ptr + (type->master_start[i] - type->master_shift)), (type->master_end[i] - type->master_start[i] + 1) * sizeof(suNf_spinor) / sizeof(double));
	apply_BCs_on_spinor_field(s);
}
