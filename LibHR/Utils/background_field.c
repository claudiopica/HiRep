/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "utils.h"
#include "suN.h"
#include "update.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "logger.h"
#include <math.h>
#include <stdlib.h>



#define PI 3.141592653589793238462643383279502884197

#define DEBUG_BACKGROUND

void apply_background_field_zdir(suNg_field* V,double Q,int n) {
int mu;
int index;
double A=0;
complex phase_bulk,phase_boundary;
complex phase;
suNg utmp[1];
static suNg_field* V_old=NULL;
int c[4];


double E=2.*PI*(double)n/(Q*(double)GLB_T*(double)GLB_Z);
int x3,x4;

// copy gauge field
V_old=alloc_gfield(&glattice);  	
suNg_field_copy(V_old,V);

start_gf_sendrecv(V_old);
complete_gf_sendrecv(V_old);

//debug
#ifdef DEBUG_BACKGROUND
static suNg_field* Vtest_old=NULL;
static suNg_field* Vtest_new=NULL;
complex tmp,tmp2;
double diff_re,diff_im;

Vtest_old=alloc_gfield(&glattice);  	
Vtest_new=alloc_gfield(&glattice);  	
suNg_field_copy(Vtest_old,V);

start_gf_sendrecv(Vtest_old);
complete_gf_sendrecv(Vtest_old);

#endif // DEBUG_BACKGROUND


for(c[0]=0; c[0]<T; c[0]++) for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	index=ipt(c[0],c[1],c[2],c[3]);
	x4 = c[0]+zerocoord[0]; 
//	lprintf("DEBUG",0," x4 %d T %d GLB_T %d c0 %d zerocoord %d %d %d %d   \n",x4,T,GLB_T,c[0],zerocoord[0],zerocoord[1],zerocoord[2],zerocoord[3]);
	for (mu=0;mu<4;++mu){


		if (mu != 3) A =0 ;
		else A=-E*(double)x4; // -E x_4 =  - 2pi n x_4/(Q*T*L)

		phase_bulk.re = cos(Q*A);
		phase_bulk.im = sin(Q*A);

		// if on the boundary 
		if (mu==0 && x4 == GLB_T-1 )
		{
			x3 = c[3]+zerocoord[3]; 
//			lprintf("DEBUG",0,"x3 %d x4 %d mu %d \n",x3,x4,mu);

			phase_boundary.re = cos(Q*E*GLB_T*x3);
			phase_boundary.im = sin(Q*E*GLB_T*x3);
			_complex_mul(phase,phase_bulk,phase_boundary);

		}
		else
		{
			phase=phase_bulk;
		}	

//		_suNg_mul(utmp[0],1.,*_4FIELD_AT(V_old,index,mu));
	utmp[0] = *_4FIELD_AT(V_old,index,mu);
		_suNg_mulc(*_4FIELD_AT(V,index,mu),phase,utmp[0]); 

	}// end loop mu

}// end loop local volume 

start_gf_sendrecv(V);
complete_gf_sendrecv(V);

// test sum of the difference between old and new plaquette for each site. Real and imaginary parts are treated independantly 
#ifdef DEBUG_BACKGROUND

diff_re = 0.;
diff_im = 0.;

suNg_field_copy(Vtest_new,V);
complete_gf_sendrecv(Vtest_new);

for(c[0]=0; c[0]<T; c[0]++) for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){

				tmp.re=0.;
				tmp.im=0.;
				tmp2.re=0.;
				tmp2.im=0.;

				index=ipt(c[0],c[1],c[2],c[3]);

				// calculate new plaquettes in direction 3,4
				cplaq(&tmp,index,3,0);

				// now divide by e^{iQE}
				phase.re = cos(Q*E);
				phase.im = sin(Q*E);
				_complex_div(tmp2,tmp,phase); 

				// restore old gauge field
				suNg_field_copy(V,Vtest_old);

				start_gf_sendrecv(V);
				complete_gf_sendrecv(V);
				
				cplaq(&tmp,index,3,0);

				// compute absolute value of the difference between old and new for the real and imaginary parts independantly
				diff_re += fabs(tmp2.re - tmp.re);
				diff_im += fabs(tmp2.im - tmp.im);


				// restore new gauge field
				suNg_field_copy(V,Vtest_new);
				
				start_gf_sendrecv(V);
				complete_gf_sendrecv(V);



}
global_sum(&diff_re,1);
global_sum(&diff_im,1);

free_gfield(Vtest_old);
free_gfield(Vtest_new);

lprintf("DEBUG",0," difference of the real and imaginary part of the old and new plaquette summed on the volume: %e %e \n",diff_re,diff_im);

#endif // DEBUG_BACKGROUND

free_gfield(V_old);
}

