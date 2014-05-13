#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "error.h"
#include "ranlux.h"
#include "random.h"
#include "global.h"
#include "communications.h"


//Random array of +/- 1/\sqrt(2)
void ranz2(double r[],int n)
{
  double plus = sqrt(.5);
  double minus = -sqrt(.5);
  int i;
  ranlxd(r,n);
  for(i = 0; i < n; i++)
    r[i] = (r[i]<.5) ? plus : minus;
}

void generate_random_point(int *pr) {
  double ran;
  ranlxd(&ran,1);
  pr[1] = (int)(ran*GLB_X);
  ranlxd(&ran,1);
  pr[2] = (int)(ran*GLB_Y);
  ranlxd(&ran,1);
  pr[3] = (int)(ran*GLB_Z);
  ranlxd(&ran,1);
  pr[0] = (int)(ran*GLB_T);
  bcast_int(pr,4);
}
