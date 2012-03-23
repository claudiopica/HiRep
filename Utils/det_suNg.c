/*******************************************************************************
*
* File  det_suNg.c
*
* Function to calculate determinant of a complex suNg-matrix using 
* LU-decomposition. See NR.
* 
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "suN.h"
#include "representation.h"

/* The incoming matrix will be destroyed */

void det_suNg(complex* res, suNg *a){
  double vv[NG];
  complex csum,ctmp,cdum;
  double big,dum,tmp;
  int i,j,k,imax;
  int d = 1;
  for (j=0;j<NG;j++){
    big = 0.0;
    for (i=0;i<NG;i++){
      tmp = _complex_abs_sqr(a->c[NG*j+i]);
      if (tmp>big) big = tmp;
    }
    if (big == 0.0) {
      _complex_0(*res);
      return;
    }
    vv[j] = 1.0/sqrt(big);
  }
  imax = 0;
  for (j=0;j<NG;j++){
    for (i=0;i<j;i++){
      csum = a->c[NG*i+j];
      for (k=0;k<i;k++){
	_complex_mul(ctmp,a->c[NG*i+k],a->c[NG*k+j]);
	_complex_sub_assign(csum,ctmp);
      }
      a->c[NG*i+j] = csum;
    }
    big = 0.0;
    for (i=j;i<NG;i++){
      csum = a->c[NG*i+j];
      for (k=0;k<j;k++){
	_complex_mul(ctmp,a->c[NG*i+k],a->c[NG*k+j]);
	_complex_sub_assign(csum,ctmp);	  
      }
      a->c[NG*i+j] = csum;
      dum=vv[i]*_complex_abs(csum);
      if (dum >= big){
	big=dum;
	imax = i;
      }
    }

    if (j != imax){
      for (k=0;k<NG;k++){
	cdum = a->c[NG*imax+k];
	a->c[NG*imax+k] = a->c[NG*j+k];
	a->c[NG*j+k]=cdum;
      }
      d = -d;
      vv[imax] = vv[j];
    }

    if (_complex_abs(a->c[NG*j+j]) == 0.0){
      a->c[NG*j+j].re = 1e-20;
      a->c[NG*j+j].im = 0;	
    }

    if (j != NG -1 ){
      _complex_inv(cdum,a->c[NG*j+j]);
      for (i=j+1;i<NG;i++){
	_complex_mul(ctmp,a->c[NG*i+j],cdum);
	a->c[NG*i+j] = ctmp;
      }
    }
  }
  csum.re = (double) d;
  csum.im = 0.0;
  for (j=0;j<NG;j++){
    _complex_mul(ctmp,a->c[NG*j+j],csum);
    csum = ctmp;
  }
  *res = csum;
}
