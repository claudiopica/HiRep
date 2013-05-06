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
#include "utils.h"

#define Hind(i,j) ((i)*kry_dim - (i)*((i)+1)/2+(j))

static double pythag3(double a, double b, double l){
  double r1,r2;
  if (l>a && l>b){
    r1 = a/l;
    r2 = b/l;
    return l*sqrt(1.+r1*r1+r2*r2);
  }
  else{
    return sqrt(a*a+b*b+l*l);
  }
}

static int GMRES_core(short int *valid, inverter_par *par, spinor_operator M, spinor_field *in, spinor_field *out, spinor_field *trial){
  int kry_dim = par->kry_dim;
  int i,j,k;
  int cgiter = 0;
  double beta;
  spinor_field * w;
  spinor_field * v;
  complex  h[kry_dim*(kry_dim+1)/2];
  double ss[kry_dim];
  complex cs[kry_dim];
  complex gbar[kry_dim+1];

  double tmp;
  double norm_w;
  complex c_tmp1,c_tmp2;

  double inorm = spinor_field_sqnorm_f(in);

#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  v = alloc_spinor_field_f(kry_dim+1,in->type);
  w=v+kry_dim;
  alloc_mem_t=std_mem_t; /* set the allocation memory type back */

  spinor_field_copy_f(&v[0], in);  // p2 now has the rhs-spinorfield "b"
    
  if(trial!=NULL) {
    M.dbl(w,trial);
    ++cgiter;
    spinor_field_sub_assign_f(&v[0],w);
    
    if(out!=trial){
      spinor_field_copy_f(out,trial);
    }
    
  } else {
    spinor_field_zero_f(out);
  }

  beta=sqrt(spinor_field_sqnorm_f(&v[0]));

  lprintf("INVERTER",100,"GMRES Initial norm: %e\n",beta);

  spinor_field_mul_f(&v[0],1./beta,&v[0]);
  gbar[0].re=beta; 
  for (j=0;j<kry_dim;++j){
    ++cgiter;
    M.dbl(w,&v[j]);
    
    tmp = sqrt(spinor_field_sqnorm_f(w));
    for (i=0;i<j+1;++i){	// Gram-Schmidth orthogonalization (modified)
      h[Hind(i,j)]=spinor_field_prod_f(w,&v[i]);
      _complex_minus(c_tmp1,h[Hind(i,j)]);
      spinor_field_mulc_add_assign_f(w,c_tmp1,&v[i]);
    }
    norm_w = sqrt(spinor_field_sqnorm_f(w));

    /*    if ( norm_w/tmp < 0.001){
     lprintf("INVERTER",100,"GMRES in second GR\n");
     for (i=0;i<j+1;++i){	// Gram-Schmidth orthogonalization (modified)
       c_tmp1=spinor_field_prod_f(w,&v[i]);
       _complex_add_assign(h[Hind(i,j)],c_tmp1);
       _complex_minus(c_tmp1,c_tmp1);
       spinor_field_mulc_add_assign_f(w,c_tmp1,&v[i]);
     }
      norm_w = sqrt(spinor_field_sqnorm_f(w));      
      }*/


    /*    if (norm_w*norm_w < 1e-28){
      lprintf("INVERTER",90,"GMRES_core: Lucky termination!:\n");
      break;
      }*/

    /* Do the givens rotations */
    //First rotate new column with old rot matrixes
    for (i=0;i<j;++i){
      c_tmp1=h[Hind(i,j)];
      _complex_mul_star(h[Hind(i,j)],c_tmp1,cs[i]);
      _complex_mulr_assign(h[Hind(i,j)],ss[i],h[Hind((i+1),j)]);
      c_tmp2=h[Hind((i+1),j)];
      _complex_mulr(h[Hind((i+1),j)],-ss[i],c_tmp1);
      _complex_mul_assign(h[Hind((i+1),j)],cs[i],c_tmp2);
    }

    //Calculate new sin and cos term 
    tmp=sqrt(_complex_prod_re(h[Hind(j,j)],h[Hind(j,j)])+norm_w*norm_w);
    ss[j]= norm_w/ tmp;
    _complex_mulr(cs[j],1./tmp,h[Hind(j,j)]);

    //Rotate gbar
    tmp=gbar[j].re;
    _complex_mulr_star(gbar[j],tmp,cs[j]);
    gbar[j+1].re = -ss[j]*tmp;

    //Rotate new column with the new rot matrix
    c_tmp1=h[Hind(j,j)];
    _complex_mul_star(h[Hind(j,j)],c_tmp1,cs[j]);
    h[Hind(j,j)].re += ss[j]*norm_w;

    lprintf("INVERTER",90,"GMRES Error in iteration %d: %e\n",j,gbar[j+1].re*gbar[j+1].re/inorm);
    /*Check Convergence*/
    *valid = (fabs(gbar[j+1].re*gbar[j+1].re/inorm) < par->err2 );
    if (*valid){
      j++;
      break;
    }

    /* If not converged pick new vector */
    if (j<kry_dim-1){
      spinor_field_mul_f(&v[j+1],1./norm_w,w);
    }
  }

  // h is now upper triangular... 
  // now the vector that minimizes |\beta e_1 - h y| is found by y = h^-1 g
  // reusing the variable cs (not needed anymore)
  for (i=j-1;i>=0;--i){
    cs[i]=gbar[i];
    for (k=j-1;k>i;--k){
      _complex_mul_sub_assign(cs[i],h[Hind(i,k)],cs[k]);
    }
    _complex_mulr(cs[i],1./h[Hind(i,i)].re,cs[i]);
  }

  //The solutions is (out_m = out_0 + V y)
  for (i=0;i<j;++i){
    spinor_field_mulc_add_assign_f(out,cs[i],&v[i]);
  }
  
  lprintf("INVERTER",50,"GMRES Error after GMRES_core: %e\n",gbar[j].re*gbar[j].re/inorm);
  
  //Debug the error
  /*  M.dbl(w,out);
  spinor_field_sub_assign_f(w,in);
  double tau=spinor_field_sqnorm_f(w)/inorm;
  lprintf("INVERTER",100,"test  = %e (req. %e)\n",tau,par->err2);*/

  /* free memory */
  free_spinor_field_f(v);
  /* return number of cg iter */
  return cgiter;
}

static int GMRES_core_flt(short int *valid, inverter_par *par , spinor_operator M, spinor_field_flt *in, spinor_field_flt *out, spinor_field_flt *trial){
  int kry_dim = par->kry_dim;
  int i,j,k;
  int cgiter = 0;
  double beta;
  spinor_field_flt * w;
  spinor_field_flt * v;
  complex  h[kry_dim*(kry_dim+1)/2];
  double ss[kry_dim];
  complex cs[kry_dim];
  complex gbar[kry_dim+1];

  double tmp;
  double norm_w;
  complex c_tmp1,c_tmp2;
  double inorm = spinor_field_sqnorm_f_flt(in);

#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  v = alloc_spinor_field_f_flt(kry_dim+1,in->type);
  w=v+kry_dim;
  alloc_mem_t=std_mem_t; /* set the allocation memory type back */

  spinor_field_copy_f_flt(&v[0], in);  // p2 now has the rhs-spinorfield "b"
    
  if(trial!=NULL) {
    M.flt(w,trial);
    ++cgiter;
    spinor_field_sub_assign_f_flt(&v[0],w);
    
    if(out!=trial){
      spinor_field_copy_f_flt(out,trial);
    }
    
  } else {
    spinor_field_zero_f_flt(out);
  }

  beta=sqrt(spinor_field_sqnorm_f_flt(&v[0]));

  lprintf("INVERTER",100,"GMRES_core_flt Initial norm: %e\n",beta);

  spinor_field_mul_f_flt(&v[0],1./beta,&v[0]);
  gbar[0].re=beta; 
  for (j=0;j<kry_dim;++j){
    ++cgiter;
    M.flt(w,&v[j]);

    /* Gram-Schmidth orthogonalization (modified)  */
    tmp = sqrt(spinor_field_sqnorm_f_flt(w));
    for (i=0;i<j+1;++i){	
      h[Hind(i,j)]=spinor_field_prod_f_flt(w,&v[i]);
      _complex_minus(c_tmp1,h[Hind(i,j)]);
      spinor_field_mulc_add_assign_f_flt(w,c_tmp1,&v[i]);
    }
    norm_w = sqrt(spinor_field_sqnorm_f_flt(w));

    if ( norm_w/tmp < 10e-1){
      for (i=0;i<j+1;++i){	
	c_tmp1=spinor_field_prod_f_flt(w,&v[i]);
	//	_complex_add_assign(h[i*kry_dim+j],c_tmp1);
	_complex_minus(c_tmp1,c_tmp1);
	spinor_field_mulc_add_assign_f_flt(w,c_tmp1,&v[i]);
      }
      norm_w = sqrt(spinor_field_sqnorm_f_flt(w));
    }

    if (norm_w*norm_w < par->err2){
      lprintf("INVERTER",500,"GMRES_core_flt: Lucky termination!:\n");
      break;
    }

    /* Do the givens rotations */
    //First rotate new column with old rot matrixes
    for (i=0;i<j;++i){
      c_tmp1=h[Hind(i,j)];
      _complex_mul_star(h[Hind(i,j)],c_tmp1,cs[i]);
      _complex_mulr_assign(h[Hind(i,j)],ss[i],h[Hind((i+1),j)]);
      c_tmp2=h[Hind((i+1),j)];
      _complex_mulr(h[Hind((i+1),j)],-ss[i],c_tmp1);
      _complex_mul_assign(h[Hind((i+1),j)],cs[i],c_tmp2);
    }

    //Calculate new sin and cos term 
    //    tmp=sqrt(_complex_prod_re(h[Hind(j,j)],h[Hind(j,j)])+norm_w*norm_w);
    tmp = pythag3(h[Hind(j,j)].re,h[Hind(j,j)].im,norm_w);
    ss[j]= norm_w/ tmp;
    _complex_mulr(cs[j],1./tmp,h[Hind(j,j)]);

    //Rotate gbar
    tmp=gbar[j].re;
    _complex_mulr_star(gbar[j],tmp,cs[j]);
    gbar[j+1].re = -ss[j]*tmp;

    //Rotate new column with the new rot matrix
    c_tmp1=h[Hind(j,j)];
    _complex_mul_star(h[Hind(j,j)],c_tmp1,cs[j]);
    h[Hind(j,j)].re += ss[j]*norm_w;

    lprintf("INVERTER",100,"GMRES_core_flt Error in iteration %d: %e\n",j,gbar[j+1].re*gbar[j+1].re/inorm);
    /*Check Convergence*/
    *valid = (fabs(gbar[j+1].re*gbar[j+1].re/inorm) < par->err2_flt );
    if (*valid){
      j++;
      break;
    }

    /* If not converged pick new vector */
    if (j<kry_dim-1){
      spinor_field_mul_f_flt(&v[j+1],1./norm_w,w);
    }
  }

  // h is now upper triangular... 
  // now the vector that minimizes |\beta e_1 - h y| is found by y = h^-1 g
  // reusing the variable cs (not needed anymore)
  for (i=j-1;i>=0;--i){
    cs[i]=gbar[i];
    for (k=j-1;k>i;--k){
      _complex_mul_sub_assign(cs[i],h[Hind(i,k)],cs[k]);
    }
    _complex_mulr(cs[i],1./h[Hind(i,i)].re,cs[i]);
  }

  //The solutions is (out_m = out_0 + V y)
  for (i=0;i<j;++i){
    spinor_field_mulc_add_assign_f_flt(out,cs[i],&v[i]);
  }
  
  lprintf("INVERTER",50,"GMRES Error after GMRES_core_flt: %e\n",gbar[j].re*gbar[j].re/inorm);
  
  //  Debug the error
  /*  M.flt(w,out);
  spinor_field_sub_assign_f_flt(w,in);
  double tau=spinor_field_sqnorm_f_flt(w)/inorm;
  lprintf("INVERTER",100,"test  = %e (req. %e)\n",tau,par->err2);*/

  /* free memory */
  free_spinor_field_f_flt(v);
  /* return number of cg iter */
  return cgiter;
}

int GMRES(inverter_par *par, spinor_operator M, spinor_field *in, spinor_field *out, spinor_field *trial){
  int iter, rep=0;
  short int valid=0;

  iter = GMRES_core(&valid,par,M,in,out,trial);
  while (!valid && (par->max_iter==0 || iter < par->max_iter)){
    iter += GMRES_core(&valid,par,M,in,out,out);
    if ((++rep %10 == 0 )){
      lprintf("INVERTER",-10,"GMRES recursion = %d (precision too high?)\n",rep);
    }
  }
  lprintf("INVERTER",10,"GMRES: MVM = %d\n",iter);
  return iter;
}




int GMRES_flt(inverter_par *par, spinor_operator M, spinor_field_flt *in, spinor_field_flt *out, spinor_field_flt *trial){
  int iter, rep=0;
  short int valid=0;

  iter = GMRES_core_flt(&valid,par,M,in,out,NULL);
  while (!valid && (par->max_iter==0 || iter < par->max_iter)){
    iter += GMRES_core_flt(&valid,par,M,in,out,out);
    if ((++rep %10 == 0 )){
      lprintf("INVERTER",-10,"GMRES recursion = %d (precision too high?)\n",rep);
    }
  }
  lprintf("INVERTER",10,"GMRES_flt: MVM_flt = %d\n",iter);
  return iter;
}

