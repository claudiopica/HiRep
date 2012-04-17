/***************************************************************************\
* Copyright (c) 2012, Ulrik Ishøj Søndergaard                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/
#ifdef WITH_GPU

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"
#include "memory.h"

#include <stdio.h>
#include <math.h>

extern rhmc_par _update_par;

template<unsigned int is_ix_odd> // is_ix_odd = 0 if ix<vol4h,      is_ix_even = 1 if ix>vol4h
__device__ void staples_device(int ix,int mu,suNg *v, const int *iup, const int *idn, const suNg* gauge,const int vol4h)
{
	suNg u1up,u2up,u3up;
	suNg u1dn,u2dn,u3dn;
	suNg staple, tr1, tr2;
	
	int i,nu,ixpmu,ixpnu,ixmnu,ixpmumnu;
	unsigned int is_ixpmu_odd = (is_ix_odd)? 0 : 1; 
	unsigned int is_ixpnu_odd = (is_ix_odd)? 0 : 1; 
	unsigned int is_ixmnu_odd = (is_ix_odd)? 0 : 1; 
	unsigned int is_ixpmumnu_odd = (is_ix_odd)? 1 : 0; 
	
	_suNg_zero(*v);
	ixpmu=iup(ix,mu);

	
	for (i=1;i<4;i++)
	{
		nu=(mu+i)&0x3;
		ixpnu=iup(ix,nu);
		ixmnu=idn(ix,nu);
		ixpmumnu=idn(ixpmu,nu);  // if ix is even the this is also even. -"- odd

		ixpnu	-=	is_ixpnu_odd*vol4h;
		ixmnu	-=	is_ixmnu_odd*vol4h;
		ixpmumnu -=	is_ixpmumnu_odd*vol4h; 
		
		
		
		_suNg_read_gpu(vol4h,u1up,gauge+is_ix_odd*4*vol4h,ix-is_ix_odd*vol4h,nu);
		_suNg_read_gpu(vol4h,u2up,gauge+is_ixpnu_odd*4*vol4h,ixpnu,mu);
		_suNg_read_gpu(vol4h,u3up,gauge+is_ixpmu_odd*4*vol4h,ixpmu-is_ixpmu_odd*vol4h,nu);

		_suNg_read_gpu(vol4h,u1dn,gauge+is_ixmnu_odd*4*vol4h,ixmnu,nu);
		_suNg_read_gpu(vol4h,u2dn,gauge+is_ixmnu_odd*4*vol4h,ixmnu,mu);
		_suNg_read_gpu(vol4h,u3dn,gauge+is_ixpmumnu_odd*4*vol4h,ixpmumnu,nu);

		
		//up_staple();
		_suNg_times_suNg(tr2,u1up,u2up);
		_suNg_dagger(tr1,u3up);
		_suNg_times_suNg(staple,tr2,tr1);
		
		//add_to_v(v);
		  _suNg_add_assign(*v,staple);
		
		//dn_staple();
		_suNg_times_suNg(tr2,u2dn,u3dn);
		_suNg_dagger(tr1,u1dn);
		_suNg_times_suNg(staple,tr1,tr2);
		
		//add_to_v(v);
		_suNg_add_assign(*v,staple);
	}
}




__global__ void gauge_force_kernel(const suNg* gauge, suNg_algebra_vector* force, const int *iup, const int *idn, int N, double dt_beta_over_NG){
	
	suNg s1,s2,s_tmp;
	suNg_algebra_vector f;
//	double nsq;
	int mu;
	int vol4h=N/2;
	
	int ix = blockIdx.x*blockDim.x+ threadIdx.x;
	ix = min(ix,N-1);
	unsigned int is_ix_odd = (ix>=vol4h)? 1 : 0; 
	
	for (mu=0; mu<4; ++mu) {
		
		// staples
		if (is_ix_odd){	staples_device<1>(ix,mu,&s1, iup, idn, gauge,vol4h);}
		else {		staples_device<0>(ix,mu,&s1, iup, idn, gauge,vol4h);}

		_suNg_read_gpu(vol4h,s_tmp,gauge+is_ix_odd*4*vol4h,ix-is_ix_odd*vol4h,mu); 
		_suNg_times_suNg_dagger(s2,s_tmp,s1);
		
		_fund_algebra_project(f,s2);

		_algebra_vector_mul_add_assign_gpu_g(vol4h,force+is_ix_odd*4*vol4h,ix-is_ix_odd*vol4h,mu, -dt_beta_over_NG, f); 

	}
}

void force0(double dt, suNg_av_field *force, void *vpar){
	
  //	gfield_copy_to_gpu(u_gauge);
  //	suNg_av_field_copy_to_gpu(force);
	/* check input types */
	_TWO_SPINORS_MATCHING(u_gauge,force);
	
	int N = T*X*Y*Z;//u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
	int grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);

	gauge_force_kernel<<<grid,BLOCK_SIZE>>>(u_gauge->gpu_ptr, force->gpu_ptr, iup_gpu, idn_gpu, N,dt*_update_par.beta/((double)(NG)));

	
	//	suNg_av_field_copy_from_gpu(force);
  }


#endif //WITH_GPU


