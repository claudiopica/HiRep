#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "representation.h"
#include "random.h"
#include <math.h>

void gaussian_momenta(suNg_algebra_vector *momenta) {
   int i;

   const double c3=1./sqrt(_FUND_NORM2);
   const int ngen=NG*NG-1;
   
   gauss((double*)momenta,ngen*4*VOLUME);
   for (i=0; i<ngen*4*VOLUME; ++i) {
	*(((double*)momenta)+i)*=c3;
   }
}
