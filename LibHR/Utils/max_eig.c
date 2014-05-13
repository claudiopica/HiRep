/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "utils.h"
#include "linear_algebra.h"
#include "update.h"
#include "random.h"
#include "memory.h"
#include "geometry.h"
#include "global.h"
#include <math.h>
#include <stdlib.h>
#include "logger.h"


/* use power method to find max eigvalue */
int max_H(spinor_operator H, geometry_descriptor *type, double *max) {
  double norm, oldmax, dt;
  spinor_field *s1,*s2,*s3;
  int count;

  s1=alloc_spinor_field_f(3,type);
/* #ifdef UPDATE_EO */
/*  s1=alloc_spinor_field_f(3,&glat_even); */
/* #else */
/*   s1=alloc_spinor_field_f(3,&glattice); */
/* #endif */
  s2=s1+1;
  s3=s2+1;

  /* random spinor */
  gaussian_spinor_field(s1);
  norm=sqrt(spinor_field_sqnorm_f(s1));
  spinor_field_mul_f(s1,1./norm,s1);
  norm=1.;

  dt=1.;

  H(s3,s1); count=1;

  do {

    spinor_field_mul_f(s1,dt,s3);

    norm=sqrt(spinor_field_sqnorm_f(s1));
    spinor_field_mul_f(s1,1./norm,s1);

    oldmax=*max;
    H(s3,s1); ++count;
    *max=spinor_field_prod_re_f(s1,s3);
    
    /* printf("Iter %d: %4.5e\n",count,fabs(oldnorm-norm)); */
  } while (fabs((*max-oldmax)/(*max))>1.e-3);

  *max*=1.1; /* do not know exact bound */

  lprintf("MaxH",10,"Max_eig = %1.8e [MVM = %d]\n",*max,count); 

  free_spinor_field_f(s1);
  
  return count;
}



