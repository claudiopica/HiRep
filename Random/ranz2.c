#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "error.h"


void rz2_init(int seed)
{
  error(seed<=0,1,"rz2_init [ranz2.c]",
	"Bad choice of seed (should be between 1 and 2^31-1)");

  srand(seed);
}


void ranz2(double r[],int n)
{
  double plus = sqrt(.5);
  double minus = -sqrt(.5);
  unsigned int i;
  for(i = 0; i < n; i++)
    r[i] = ((rand() & 1) == 0) ? plus : minus;
}

