/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "dirac.h"
#include "suN.h"
#include "random.h"
#include "linear_algebra.h"
#include "update.h"
#include <math.h>

void gaussian_spinor_field(spinor_field *s) {
	const double c1=1./sqrt(2.);
	unsigned int i;
	geometry_descriptor *type = s->type;
	for(i=0;i<type->local_master_pieces;i++)
 	  gauss((double*)(s->ptr+type->master_start[i]),(type->master_end[i]-type->master_start[i]+1)*sizeof(suNf_spinor)/sizeof(double));
	spinor_field_mul_f(s,c1,s);
#if defined(BASIC_SF) || defined(ROTATED_SF)
        SF_spinor_bcs(s);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */
}

void gaussian_spinor_field_flt(spinor_field_flt *s) {
	const float c1=1./sqrt(2.);
	unsigned int i;
	geometry_descriptor *type = s->type;
	for(i=0;i<type->local_master_pieces;i++)
 	  gauss_flt((float*)(s->ptr+type->master_start[i]),(type->master_end[i]-type->master_start[i]+1)*sizeof(suNf_spinor_flt)/sizeof(float));
	spinor_field_mul_f_flt(s,c1,s);
#if defined(BASIC_SF) || defined(ROTATED_SF)
        SF_spinor_bcs_flt(s);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */
}


void z2_spinor_field(spinor_field *s) {
	unsigned int i;
	geometry_descriptor *type = s->type;
	for(i=0;i<type->local_master_pieces;i++)
 	  ranz2((double*)(s->ptr+type->master_start[i]),(type->master_end[i]-type->master_start[i]+1)*sizeof(suNf_spinor)/sizeof(double));
#if defined(BASIC_SF) || defined(ROTATED_SF)
        SF_spinor_bcs(s);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */
}

