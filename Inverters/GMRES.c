/***************************************************************************\
* Copyright (c) 2013, Ari Hietanen, Ulrik Soendergaard, Claudio Pica        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "linear_algebra.h"
#include "complex.h"
#include "memory.h"
#include "update.h"
#include "logger.h"
#include "communications.h"
#include "global.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int GMRES_core(short int *valid, MINRES_par *par, int kry_dim , spinor_operator M, spinor_field *in, spinor_field *out, spinor_field *trial){

int j;
int cgiter = 0;
double beta;
spinor_field * w;
spinor_field * v;
spinor_field * p;


#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  w = alloc_spinor_field_f(kry_dim+2,in->type);
  v=w+kry_dim;
  p=w+kry_dim+1;
  alloc_mem_t=std_mem_t; /* set the allocation memory type back */



  spinor_field_copy_f(v, in);  // p2 now has the rhs-spinorfield "b"
  if(trial!=NULL) {
    M.dbl(p,trial);
    ++cgiter;
    spinor_field_sub_assign_f(v,p);
    
    if(out!=trial){
      spinor_field_copy_f(out,trial);
    }
    
  } else {
    spinor_field_zero_f(out);
  }

  beta=sqrt(spinor_field_sqnorm_f(v));
  spinor_field_mul_f(v,1./beta,v);

for (j=0;j<kry_dim+1;++j){

}





  /* free memory */
  free_spinor_field_f(w);

  /* return number of cg iter */
  return cgiter;

}