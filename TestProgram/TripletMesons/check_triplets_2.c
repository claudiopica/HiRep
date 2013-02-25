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


#error "Old version of Mesons, it should be updated"

void print_mat(complex mat[4][4], const char name[]) {
  int i,j;
  lprintf("MAIN",0,"%s = \n", name);
  for(i=0; i<4; i++) {
    lprintf("MAIN",0,"[ ");
    for(j=0; j<4; j++) {
      lprintf("MAIN",0,"(%.2f,%.2f) ",mat[i][j].re,mat[i][j].im);
    }
    lprintf("MAIN",0,"]\n");
  }
}


#define set_zero_mat(A) \
  { \
    int _i,_j; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      A[_i][_j].re = A[_i][_j].im = 0.; \
    } \
  }

#define set_unit_mat(A) \
  set_zero_mat(A); \
  { \
    int _i; \
    for(_i=0; _i<4; _i++) { \
      A[_i][_i].re = 1.; \
    } \
  }

#define copy_mat(A,B) \
  { \
    int _i,_j; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      A[_i][_j].re = B[_i][_j].re; \
      A[_i][_j].im = B[_i][_j].im; \
    } \
  }

#define sub_mat(A,B) \
  { \
    int _i,_j; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      A[_i][_j].re -= B[_i][_j].re; \
      A[_i][_j].im -= B[_i][_j].im; \
    } \
  }

#define sqnorm_mat(ret,A) \
  ret=0.; \
  { \
    int _i, _j; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      ret += A[_i][_j].re*A[_i][_j].re + A[_i][_j].im*A[_i][_j].im; \
    } \
  }

#define adj_mat(A) \
  { \
    int _i, _j; \
    for(_i=0; _i<4; _i++) { \
      A[_i][_i].im=-A[_i][_i].im; \
      for(_j=_i+1; _j<4; _j++) { \
        double _tmp; \
        _tmp=A[_i][_j].re; \
        A[_i][_j].re=A[_j][_i].re; \
        A[_j][_i].re=_tmp; \
        _tmp=A[_i][_j].im; \
        A[_i][_j].im=-A[_j][_i].im; \
        A[_j][_i].im=-_tmp; \
      } \
    } \
  }

#define CMA(x,y,z) \
  (x).re += (y).re*(z).re-(y).im*(z).im; \
  (x).im += (y).re*(z).im+(y).im*(z).re;
  

#define mult_mat(A,B) \
  { \
    int _i, _j, _k; \
    complex wm[4][4]; \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      wm[_i][_j].re = wm[_i][_j].im = 0.; \
      for(_k=0; _k<4; _k++) { \
        CMA(wm[_i][_j],A[_i][_k],B[_k][_j]);\
      } \
    } \
    for(_i=0; _i<4; _i++) \
    for(_j=0; _j<4; _j++) { \
      A[_i][_j].re = wm[_i][_j].re; \
      A[_i][_j].im = wm[_i][_j].im; \
    } \
  }

#define trace_mat(ret,A) \
  { \
    ret.re=ret.im=0.; \
    int _i; \
    for(_i=0; _i<4; _i++) { \
      ret.re += A[_i][_i].re; \
      ret.im += A[_i][_i].im; \
    } \
  }

#define mult_mat_spinor(out, A, in) \
  { \
    int _a,_i,_j; \
    for(_a=0;_a<NF;_a++) \
    for(_i=0; _i<4; _i++) { \
      out.c[_i].c[_a].re=out.c[_i].c[_a].im=0.; \
      for(_j=0; _j<4; _j++) { \
        CMA(out.c[_i].c[_a],A[_i][_j],in.c[_j].c[_a]); \
      } \
    } \
  }


int main(int argc,char *argv[])
{
	complex gamma[5][4][4];
	complex Gamma[4][4];
	complex test[4][4];
	complex rmat[4][4];
	complex trace, ctest;
	int sign;
	double norm2;
	suNf_spinor in, out, stest;

  rlxd_init(1,time(NULL));
  gauss((double*)rmat,32);

  lprintf("MAIN",0,"*********************************************************\n");
  lprintf("MAIN",0,"Checking GAMMA_trace_H functions...\n");
  lprintf("MAIN",0,"*********************************************************\n");


  g5_debug(gamma[4],&sign);
  if(sign != -1) lprintf("MAIN",0,"ERROR! Bad sign for gamma_5!\n");
  print_mat(gamma[4],"gamma_5");
  g5_trace_H(&trace,rmat[0]);
  trace_mat(ctest,rmat);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g0_debug(gamma[0],&sign);
  if(sign != 1) lprintf("MAIN",0,"ERROR! Bad sign for gamma_0!\n");
  print_mat(gamma[0],"gamma_0");
  g0_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[0]);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }
  
  g1_debug(gamma[1],&sign);
  if(sign != -1) lprintf("MAIN",0,"ERROR! Bad sign for gamma_1!\n");
  print_mat(gamma[1],"gamma_1");
  g1_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[1]);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g2_debug(gamma[2],&sign);
  if(sign != -1) lprintf("MAIN",0,"ERROR! Bad sign for gamma_2!\n");
  print_mat(gamma[2],"gamma_2");
  g2_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[2]);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g3_debug(gamma[3],&sign);
  if(sign != -1) lprintf("MAIN",0,"ERROR! Bad sign for gamma_3!\n");
  print_mat(gamma[3],"gamma_3");
  g3_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[3]);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }



  id_debug(Gamma, &sign);
  set_unit_mat(test);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != 1) {
    print_mat(Gamma, "id");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  id_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for id! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g0g1_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[1]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != 1) {
    print_mat(Gamma, "g0g1");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g1_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g1! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g0g2_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[2]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != 1) {
    print_mat(Gamma, "g0g2");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g2_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g2! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g0g3_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[3]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != 1) {
    print_mat(Gamma, "g0g3");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g3_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g3! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g0g5_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[4]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != 1) {
    print_mat(Gamma, "g0g5");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g5_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g5! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g5g1_debug(Gamma, &sign);
  copy_mat(test,gamma[4]);
  mult_mat(test,gamma[1]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != -1) {
    print_mat(Gamma, "g5g1");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g5g1_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g5g1! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g5g2_debug(Gamma, &sign);
  copy_mat(test,gamma[4]);
  mult_mat(test,gamma[2]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != -1) {
    print_mat(Gamma, "g5g2");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g5g2_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g5g2! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g5g3_debug(Gamma, &sign);
  copy_mat(test,gamma[4]);
  mult_mat(test,gamma[3]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != -1) {
    print_mat(Gamma, "g5g3");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g5g3_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g5g3! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g0g5g1_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[1]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != -1) {
    print_mat(Gamma, "g0g5g1");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g5g1_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g5g1! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g0g5g2_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[2]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != -1) {
    print_mat(Gamma, "g0g5g2");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g5g2_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g5g2! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
  }

  g0g5g3_debug(Gamma, &sign);
  copy_mat(test,gamma[0]);
  mult_mat(test,gamma[4]);
  mult_mat(test,gamma[3]);
  sub_mat(test,Gamma);
  sqnorm_mat(norm2,test);
  if(norm2 > 1e-15 || sign != -1) {
    print_mat(Gamma, "g0g5g3");
    lprintf("MAIN",0,"sign = %d\n",sign);
  }
  g0g5g3_trace_H(&trace,rmat[0]);
  copy_mat(test,rmat);
  mult_mat(test,gamma[4]);
  mult_mat(test,Gamma);
  trace_mat(ctest,test);
  ctest.re-=trace.re;ctest.im-=trace.im;
  if(ctest.re > 1e-15 || ctest.im > 1e-15) {
    lprintf("MAIN",0,"ERROR! Trace mismatch for g0g5g3! trace=(%f,%f) ctest=(%e,%e)\n",trace.re,trace.im,ctest.re,ctest.im);
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
  if(norm2 > 1e-15) {
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
  if(norm2 > 1e-15) {
    lprintf("MAIN",0,"ERROR! Mismatch for g0! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g1_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[1]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > 1e-15) {
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
  if(norm2 > 1e-15) {
    lprintf("MAIN",0,"ERROR! Mismatch for g2! norm2=%e\n",norm2);
  }

  ranlxd((double*)(&in),sizeof(suNf_spinor)/sizeof(double));
  _spinor_zero_f(out);
  g3_eval_g5GammaDag_times_spinor(&out,&in);
  copy_mat(Gamma,gamma[3]);
  mult_mat(Gamma,gamma[4]);
  adj_mat(Gamma);
  mult_mat_spinor(stest,Gamma,in);
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > 1e-15) {
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
  if(norm2 > 1e-15) {
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
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > 1e-15) {
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
  if(norm2 > 1e-15) {
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
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > 1e-15) {
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
  if(norm2 > 1e-15) {
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
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > 1e-15) {
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
  if(norm2 > 1e-15) {
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
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > 1e-15) {
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
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > 1e-15) {
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
  if(norm2 > 1e-15) {
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
  _spinor_sub_assign_f(stest,out);
  _spinor_prod_re_f(norm2,stest,stest);
  if(norm2 > 1e-15) {
    lprintf("MAIN",0,"ERROR! Mismatch for g0g5g3! norm2=%e\n",norm2);
  }

  lprintf("MAIN",0,"End of tests\n");



  exit(0);
}

