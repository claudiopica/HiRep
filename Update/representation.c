#include "global.h"
#include "representation.h"
#include "utils.h"

void represent_gauge_field() {
#ifndef REPR_FUNDAMENTAL
   int i;
   suNf *Ru=u_gauge_f;
   suNg *u=u_gauge;

   for (i=0; i<VOLUME*4; ++i){
      _group_represent(*Ru,*u);
      ++Ru;
      ++u;
   }
	 apply_bc();
#else
	static short int first_time=1;
	 if(first_time) {
		 first_time=0;
	   u_gauge_f=(suNf *)((void*)u_gauge);
		 apply_bc();
	 }
#endif
}
