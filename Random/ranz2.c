#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "error.h"
#include "ranlux.h"
#include "random.h"


void ranz2(double r[],int n)
{
  double plus = sqrt(.5);
  double minus = -sqrt(.5);
  int i;
  ranlxd(r,n);
  for(i = 0; i < n; i++)
    r[i] = (r[i]<.5) ? plus : minus;
}

