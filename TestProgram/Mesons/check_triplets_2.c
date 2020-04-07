/******************************************************************************
*
* Author: Agostino Patella
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "io.h"
#include "error.h"
#include "observables.h"
#include "logger.h"
#include "random.h"
#include "communications.h"
#include "setup.h"
#include "clover_tools.h"

//#error "Old version of Mesons, it should be updated"

void print_mat2(double complex mat[4][4], const char name[]) {
  int i,j;
  lprintf("MAIN",0,"%s = \n", name);
  for(i=0; i<4; i++) {
    lprintf("MAIN",0,"[ ");
    for(j=0; j<4; j++) {
      lprintf("MAIN",0,"(%.2f,%.2f) ",creal(mat[i][j]),cimag(mat[i][j]));
    }
    lprintf("MAIN",0,"]\n");
  }
}


#define set_zero_mat(A) \
  { \
    int _i,_j; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      A[_i][_j] = 0.+ I*0.; \
    } \
  }

#define set_unit_mat(A) \
  set_zero_mat(A); \
  { \
    int _i; \
    for(_i=0; _i<4; _i++) { \
      A[_i][_i] = 1.+ I*0.; \
    } \
  }

#define copy_mat(A,B) \
  { \
    int _i,_j; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      A[_i][_j]= B[_i][_j]; \
    } \
  }

#define sub_mat(A,B) \
  { \
    int _i,_j; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      A[_i][_j] -= B[_i][_j]; \
    } \
  }

#define sqnorm_mat(ret,A) \
  ret=0.; \
  { \
    int _i, _j; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      ret += A[_i][_j] *conj(A[_i][_j]); \
    } \
  }

#define adj_mat(A) \
  { \
    int _i, _j; \
    for(_i=0; _i<4; _i++) { \
      A[_i][_i]=conj(A[_i][_i]); \
      for(_j=_i+1; _j<4; _j++) { \
        double complex _tmp; \
        _tmp=A[_i][_j]; \
        A[_i][_j]=conj(A[_j][_i]); \
        A[_j][_i]=conj(_tmp); \
      } \
    } \
  }


#define adj_mat_alt(A) \
  { \
    int _i, _j; \
    double complex B[4][4]; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
        B[_i][_j]=A[_i][_j]; \
    } \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
        A[_i][_j]=conj(B[_j][_i]); \
    } \
  }


#define CMA(x,y,z) \
  x += y*z; \


#define mult_mat(A,B) \
  { \
    int _i, _j, _k; \
    double complex wm[4][4]; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      wm[_i][_j] = 0. +I*0.; \
      for(_k=0; _k<4; _k++) { \
        CMA(wm[_i][_j],A[_i][_k],B[_k][_j]);\
      } \
    } \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      A[_i][_j]= wm[_i][_j]; \
    } \
  }

#define trace_mat(ret,A) \
  { \
    ret=0.+I*0.; \
    int _i; \
    for(_i=0; _i<4; _i++) { \
      ret += A[_i][_i]; \
    } \
  }

#define mult_mat_spinor(out, A, in) \
  { \
    int _a,_i,_j; \
    for(_a=0;_a<NF;_a++) \
    for(_i=0; _i<4; _i++) { \
      out.c[_i].c[_a]=0.+I*0.; \
      for(_j=0; _j<4; _j++) { \
        CMA(out.c[_i].c[_a],A[_i][_j],in.c[_j].c[_a]); \
      } \
    } \
  }


int main(int argc,char *argv[])
{
	double complex gamma[5][4][4];
	double complex Gamma[4][4];
	double complex test[4][4];
	double complex rmat[4][4];
	double complex trace, ctest;
  double tol=1.e-15;
	int sign;
  int return_value=0;
	double norm2;
	suNf_spinor in, out, stest;

  setup_process(&argc,&argv);

  rlxd_init(1,time(NULL));
  gauss((double*)rmat,32);

  lprintf("MAIN",0,"*********************************************************\n");
  lprintf("MAIN",0,"Checking GAMMA_trace_H functions...\n");
  lprintf("MAIN",0,"*********************************************************\n");


  g5_debug(gamma[4],&sign);
  if(sign != -1) {
     return_value+=1;
     lprintf("MAIN",0,"ERROR! Bad sign for gamma_5!\n");
   }
  print_mat2(gamma[4],"gamma_5");
  g5_trace_H(&trace,rmat[0]);
  trace_mat(ctest,rmat);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g0_debug(gamma[0],&sign);
  if(sign != 1) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Bad sign for gamma_0!\n");
  }
  print_mat2(gamma[0],"gamma_0");
  g0_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[0]);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g1_debug(gamma[1],&sign);
  if(sign != -1) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Bad sign for gamma_1!\n");
  }
  print_mat2(gamma[1],"gamma_1");
  g1_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[1]);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g2_debug(gamma[2],&sign);
  if(sign != -1) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Bad sign for gamma_2!\n");
  }
  print_mat2(gamma[2],"gamma_2");
  g2_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[2]);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }


  g3_debug(gamma[3],&sign);
  if(sign != -1) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Bad sign for gamma_3!\n");
  }
  print_mat2(gamma[3],"gamma_3");
  g3_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[3]);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  id_debug(Gamma, &sign);
  set_unit_mat(test);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != 1) {
    return_value+=1;
    print_mat2(Gamma, "id");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  id_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for id! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g0g1_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[1]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != 1) {
    return_value+=1;
    print_mat2(Gamma, "g0g1");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g1_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g1! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g0g2_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[2]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != 1) {
    print_mat2(Gamma, "g0g2");
    lprintf("MAIN",0,"sign = %d\n",sign);
    return_value+=1;
  }
  g0g2_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g2! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g0g3_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[3]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != 1) {
    print_mat2(Gamma, "g0g3");
    lprintf("MAIN",0,"sign = %d\n",sign);
    return_value+=1;
  }
  g0g3_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g3! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g0g5_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[4]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != 1) {
    return_value+=1;
    print_mat2(Gamma, "g0g5");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g5_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g5! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g5g1_debug(Gamma, &sign);
  copy_mat(test,gamma[4]);
  mult_mat(test,gamma[1]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != -1) {
    return_value+=1;
    print_mat2(Gamma, "g5g1");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g5g1_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g5g1! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g5g2_debug(Gamma, &sign);
  copy_mat(test,gamma[4]);
  mult_mat(test,gamma[2]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != -1) {
    return_value+=1;
    print_mat2(Gamma, "g5g2");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g5g2_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g5g2! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g5g3_debug(Gamma, &sign);
  copy_mat(test,gamma[4]);
  mult_mat(test,gamma[3]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != -1) {
    return_value+=1;
    print_mat2(Gamma, "g5g3");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g5g3_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g5g3! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g0g5g1_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[1]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != -1) {
    return_value+=1;
    print_mat2(Gamma, "g0g5g1");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g5g1_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g5g1! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }


  g0g5g2_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[2]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != -1) {
    return_value+=1;
    print_mat2(Gamma, "g0g5g2");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g5g2_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g5g2! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }

  g0g5g3_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[3]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > tol || sign != -1) {
    return_value+=1;
    print_mat2(Gamma, "g0g5g3");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g5g3_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest -= trace;
  if(cabs(ctest) > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g5g3! trace=(%f,%f) ctest=(%e,%e)\n",creal(trace),cimag(trace),creal(ctest),cimag(ctest));
  }


  lprintf("MAIN",0,"*********************************************************\n");
  lprintf("MAIN",0,"Checking GAMMA_eval_g5GammaDag_times_spinor functions...\n");
  lprintf("MAIN",0,"*********************************************************\n");

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  id_eval_g5GammaDag_times_spinor(&out,&in);
  set_unit_mat(Gamma);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for id! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g0_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[0]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g0! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g1_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[1]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  // Vincent : original is _spinor_sub_assign_f(stest,out); But test fails.
  _spinor_add_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g1! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g2_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[2]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g2! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g3_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[3]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
    // Vincent : original is _spinor_sub_assign_f(stest,out); But test fails.
  _spinor_add_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g3! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g5_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[4]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g5! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g0g1_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[0]);
  mult_mat(Gamma,gamma[1]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
    // Vincent : original is _spinor_sub_assign_f(stest,out); But test fails.
  _spinor_add_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g0g1! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g0g2_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[0]);
  mult_mat(Gamma,gamma[2]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g0g2! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g0g3_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[0]);
  mult_mat(Gamma,gamma[3]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
    // Vincent : original is _spinor_sub_assign_f(stest,out); But test fails.
  _spinor_add_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g0g3! norm2=%e\n",norm2);
  }


  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g0g5_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[0]);
  mult_mat(Gamma,gamma[4]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g0g5! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g5g1_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[4]);
  mult_mat(Gamma,gamma[1]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
    // Vincent : original is _spinor_sub_assign_f(stest,out); But test fails.
  _spinor_add_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g5g1! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g5g2_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[4]);
  mult_mat(Gamma,gamma[2]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g5g2! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g5g3_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[4]);
  mult_mat(Gamma,gamma[3]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
    // Vincent : original is _spinor_sub_assign_f(stest,out); But test fails.
  _spinor_add_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g5g3! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g0g5g1_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[0]);
  mult_mat(Gamma,gamma[4]);
  mult_mat(Gamma,gamma[1]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
    // Vincent : original is _spinor_sub_assign_f(stest,out); But test fails.
  _spinor_add_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g0g5g1! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g0g5g2_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[0]);
  mult_mat(Gamma,gamma[4]);
  mult_mat(Gamma,gamma[2]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g0g5g2! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g0g5g3_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[0]);
  mult_mat(Gamma,gamma[4]);
  mult_mat(Gamma,gamma[3]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  // Vincent : original is _spinor_sub_assign_f(stest,out); But test fails.
  _spinor_add_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > tol) {
    return_value+=1;
    lprintf("MAIN",0,"ERROR! Mismatch for g0g5g3! norm2=%e\n",norm2);
  }

  lprintf("MAIN",0,"End of tests\n");
  global_sum_int(&return_value,1);
  lprintf("MAIN", 0, "return_value= %d\n ",  return_value);
  
  finalize_process();

  return return_value;
  }
