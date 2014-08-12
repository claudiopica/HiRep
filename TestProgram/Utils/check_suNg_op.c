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

void print_cmplx(complex c){
  printf("%g",c.re);
  if (c.im>0){
    printf("+%gi",c.im);
  }
  else{
    printf("%gi",c.im);
  }
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
  complex det;
  suNg a,b,c;
  random_suNg(&a);
  random_suNg(&b);
  random_suNg(&c);
  _suNg_add_assign(a,b);
  det_suNg(&det,&a);
  printf("Det: ");
  print_cmplx(det);
  printf("\n");
  printML(&a);
  c = a;
  inv_suNg(&a);
  printML(&a);
  _suNg_times_suNg(b,a,c);
  printML(&b);
  b.c[8].re=1e-15;
  c=b;
  printML(&b);
  inv_suNg(&b);
  det_suNg(&det,&b);
  printf("Det: ");
  print_cmplx(det);
  printf("\n");
  printML(&b);
  _suNg_times_suNg(a,b,c);
  printML(&a);
  
}
