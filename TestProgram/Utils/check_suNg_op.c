/*
* NOCOMPILE = SO
*/
#define MAIN_PROGRAM

#include <stdio.h>
#include <math.h>

#include "io.h"
#include "random.h"
#include "error.h"
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"

void print_cmplx(double complex c){
  printf("%g",creal(c));
  if (cimag(c)>0){
    printf("+%gi",cimag(c));
  }
  else{
    printf("-%gi",fabs(cimag(c)));
  }
}

double  norm_suNg_minus_id(suNg* a){
  int i,j;
  double r=0.;
  for (i=0;i<NG;i++){
      for (j=1;j<NG;j++){
        #ifdef WITH_QUATERNIONS
        r +=fabs(a->c[i*NG+j]);
        #else
        r +=cabs(a->c[i*NG+j]);
        #endif
        if (i==j) r-= 1.0;
    }

  }
return r;
}
void printML(suNg* a){
  int i,j;
  printf("[");
  for (i=0;i<NG;i++){
    print_cmplx(a->c[i*NG]);
    for (j=1;j<NG;j++){
      printf(",");
      print_cmplx(a->c[i*NG+j]);
    }
    printf(";");
  }
printf("]\n");
}

int main(){
  int return_value=0;
  double complex det;
  double test;
  suNg a,b,c;
  random_suNg(&a);
  random_suNg(&b);
  random_suNg(&c); // Vincent:   NOT USED !
  det_hermNg(&det,&a);
  printf("Det: ");
  print_cmplx(det);
  if (cabs(det - 1.0) >1e-14) return_value +=1;


  _suNg_add_assign(a,b); /* a = a+ b */
  det_hermNg(&det,&a);
  printf("Det: ");
  print_cmplx(det);
  printf("\n");
  if (cabs(det) < 1e-14) return_value +=1;


  printML(&a);
  c = a;  /* c = a+b */
  inv_hermNg(&a); /* a= (a+b)^-1 does not work for SO ? */
  printML(&a);
  _suNg_times_suNg(b,a,c); /* b = (a+b)^-1 *(a+b) so it should be 1 */
  printf("Should be the idendity matrix:\n");
  printML(&b);
  test= norm_suNg_minus_id(&b);
  if (test > 1e-14) return_value +=1;

  c=b; /* c =  (a+b)^-1 *(a+b) */
  printML(&b); /* Vincent useless */
  inv_hermNg(&b); /* b = ((a+b)^-1 *(a+b))^-1 */
  det_hermNg(&det,&b);
  printf("Should be the determinant of the idendity matrix:\n");
  printf("Det: ");
  print_cmplx(det);
  if (cabs(det - 1.0) >1e-14) return_value +=1;
  printf("\n");
  printML(&b);
  _suNg_times_suNg(a,b,c); /* a = ((a+b)^-1 *(a+b))^-1  * (a+b)^-1 *(a+b) */
  printf("Should be the idendity matrix:\n");
  printML(&a);
  test= norm_suNg_minus_id(&a);
  if (test > 1e-14) return_value +=1;


  return return_value;

}
