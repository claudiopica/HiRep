#include <stdlib.h>

#include "geometry.h" 
#include "global.h" 
#include "error.h"
#include "logger.h"
#include "inverters.h" 
#include "linear_algebra.h"
#include "memory.h"

void SAP_prec(int nu, inverter_ptr inv, mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out)
{

   // lprintf("SAP_prec",40,"spinor_field_sqnorm_f(in)=%e, spinor_field_sqnorm_f(out)=%e\n",spinor_field_sqnorm_f(in),spinor_field_sqnorm_f(out));
	spinor_field *res, *tmp_spinor, *tmp_spinor2;
	res=alloc_spinor_field_f(3,in->type);
	tmp_spinor=res+1;
	tmp_spinor2=res+2;
	
	for (;nu>0;nu--)
	{
		// compute/update res (This is a bit waste of time in the first round)
		M(tmp_spinor,out);
		spinor_field_sub_f(res,in,tmp_spinor);
		
	    spinor_field_zero_f(tmp_spinor);	// Temporary 
		empty_buffers(tmp_spinor);


		// Invert black
		if (out->type==&glattice){
			res->type=&glat_black;			
  		  	tmp_spinor->type=&glat_black;
  		  	(void) inv(par, M, res, tmp_spinor);
          	res->type=&glattice;
          	tmp_spinor->type=&glattice;
		} else if (out->type==&glat_even){
			res->type=&glat_even_black;			
  		  	tmp_spinor->type=&glat_even_black;
  		  	(void) inv(par, M, res, tmp_spinor);
          	res->type=&glat_even;
          	tmp_spinor->type=&glat_even;
		} else {
			res->type=&glat_odd_black;			
  		  	tmp_spinor->type=&glat_odd_black;
  		  	(void) inv(par, M, res, tmp_spinor);
          	res->type=&glat_odd;
          	tmp_spinor->type=&glat_odd;
		}
		
		
		// Update res
		M(tmp_spinor2,tmp_spinor);
		spinor_field_sub_assign_f(res,tmp_spinor2);
		// Update solution
		spinor_field_add_assign_f(out,tmp_spinor);
		
	    spinor_field_zero_f(tmp_spinor);	// Temporary 
		empty_buffers(tmp_spinor);
		
		// Invert red
			if (out->type==&glattice){
			res->type=&glat_red;			
  		  	tmp_spinor->type=&glat_red;
  		  	(void) inv(par, M, res, tmp_spinor);
          	res->type=&glattice;
          	tmp_spinor->type=&glattice;
		} else if (out->type==&glat_even){
			res->type=&glat_even_red;			
  		  	tmp_spinor->type=&glat_even_red;
  		  	(void) inv(par, M, res, tmp_spinor);
          	res->type=&glat_even;
          	tmp_spinor->type=&glat_even;
		} else {
			res->type=&glat_odd_red;			
  		  	tmp_spinor->type=&glat_odd_red;
  		  	(void) inv(par, M, res, tmp_spinor);
          	res->type=&glat_odd;
          	tmp_spinor->type=&glat_odd;
		}

		// Update solution
		spinor_field_add_assign_f(out,tmp_spinor);
	}

// remove temporary spinors	
	free_spinor_field_f(res);
}