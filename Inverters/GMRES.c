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

int i,j;
int cgiter = 0;
double beta;
spinor_field * w;
spinor_field * v;
spinor_field * p;
complex * h;
complex c;
double s;
double tmp;
complex * gbar;
double * ybar;
complex c_tmp1,c_tmp2;

gbar=(complex*)malloc((kry_dim+1)*sizeof(gbar));
ybar=(complex*)malloc((kry_dim)*sizeof(ybar));

h=(complex*)malloc(kry_dim*(kry_dim+1)*sizeof(h));
//h(i,j)=h(i*kry_dim+j)


#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  v = alloc_spinor_field_f(kry_dim+2,in->type);
  w=v+kry_dim;
  p=v+kry_dim+1;
  alloc_mem_t=std_mem_t; /* set the allocation memory type back */



  spinor_field_copy_f(&v[0], in);  // p2 now has the rhs-spinorfield "b"
  if(trial!=NULL) {
    M.dbl(p,trial);
    ++cgiter;
    spinor_field_sub_assign_f(&v[0],p);
    
    if(out!=trial){
      spinor_field_copy_f(out,trial);
    }
    
  } else {
    spinor_field_zero_f(out);
  }

  beta=sqrt(spinor_field_sqnorm_f(&v[0]));
  spinor_field_mul_f(&v[0],1./beta,&v[0]);
  gbar[0].re=beta;
  gbar[0].im=0;

for (j=0;j<kry_dim;++j){
    
    ++cgiter;
    M.dbl(w,&v[j])
    
    for (i=0;i<j+1;++i){	// Gram-Schmidth orthogonalization
    	h[i*kry_dim+j]=spinor_field_prod_f(w,&v[i]);
    	spinor_field_mul_add_assign_f(w,-h[i*kry_dim+j],&v[i]);
    }
    h[(j+1)*kry_dim+j].re=sqrt(spinor_field_sqnorm_f(w));
    h[(j+1)*kry_dim+j].im=0;
    
    if (h[(j+1)*kry_dim+j].re < par->err2){
    	kry_dim=j;
    	break;
    }
    
    if (j<kry_dim-1){
    	spinor_field_mul_f(&v[j+1],1./h[(j+1)*kry_dim+j].re,w);
	}

	gbar[j+1]= 0;	
}

// Now doing the Givens rotations
for (i=0;i<kry_dim-1;++i){
	tmp=sqrt(  _complex_prod_re( h[i*kry_dim+i], h[i*kry_dim+i])    +h[i*(kry_dim+1)+i].re*h[i*(kry_dim+1)+i].re)
	s=h[i*(kry_dim+1)+i] / tmp;
	_complex_mulr(c,1./tmp,h[i*kry_dim+i]);
	
	// Rotate gbar	
	c_tmp1=gbar[i];
	_complex_mul_star(gbar[i],c_tmp1,c);

	c_tmp2=gbar[i+1];
	_complex_mulr(gbar[i+1],-s,c_tmp1);
	_complex_mul_assign(gbar[i+1],c,c_tmp2);
	
	
	// Rotate h
	for (j=i;j<kry_dim+1;++j){
		c_tmp1=h[i*kry_dim+j];
		_complex_mul_star(h[i*kry_dim+j],c_tmp1,c);
		_complex_mulr_assign(h[i*kry_dim+j],s,h[(i+1)*kry_dim+j]);	
		c_tmp2=h[(1+i)*kry_dim+j];
		_complex_mulr(h[(i+1)*kry_dim+j],-s,c_tmp1);
		_complex_mul_assign(h[(i+1)*kry_dim+j],c,c_tmp2);
	}
}


// Hopefully h is now upper triangular... now to invert it.
for (i=kry_dim-1;i>=0;--i){
  ybar[i]=gbar[i].re;
  for (j=kry_dim-1;j>i;--j){
  	ybar[i]-= _complex_mul_re(h[i*kry_dim+j],ybar[j]);
  }
  ybar /= h[i*kry_dim+i].re;
}

//
for (i=0;i<kry_dim;++i){
	spinor_field_mul_add_assign_f(out,ybar[i],&v[i]);	
}

  /* free memory */
  free_spinor_field_f(v);

  /* return number of cg iter */
  return cgiter;

}