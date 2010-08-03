/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include <math.h>

#include "global.h"
#include "utils.h"
#include "suN.h"

#ifdef ROTATED_SF
#include "update.h"
extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif /* ROTATED_SF */


void apply_bc(){
#if defined(ANTIPERIODIC_BC_T) && !defined(ROTATED_SF) && !defined(BASIC_SF)
	if(COORD[0]==0) {
		int index;
		int ix,iy,iz;
		suNf *u;
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(2*T_BORDER,ix,iy,iz);
			if(index!=-1) {
				u=pu_gauge_f(index,0);
				_suNf_minus(*u,*u);
			}
		}
	}
#elif defined(ROTATED_SF)
#warning SCRIVERE UN TEST PER VERIFICARE CHE QUESTO EE CONSISTENTE CON I BORDI
	if(COORD[0] == 0) {
		int index;
		int ix,iy,iz;
		suNf *u;
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T_BORDER+1,ix,iy,iz);
			if(index!=-1) {
				if(ix!=X_EXT-1){
					u=pu_gauge_f(index,1);
					_suNf_mul(*u,_update_par.SF_ds,*u);
				}
				if(iy!=Y_EXT-1){
					u=pu_gauge_f(index,2);
					_suNf_mul(*u,_update_par.SF_ds,*u);
				}
				if(iz!=Z_EXT-1){
					u=pu_gauge_f(index,3);
					_suNf_mul(*u,_update_par.SF_ds,*u);
				}
			}
		}
	}
	if(COORD[0] == NP_T-1) {
		int index;
		int ix,iy,iz;
		suNf *u;
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T+T_BORDER-1,ix,iy,iz);
			if(index!=-1) {
				if(ix!=X_EXT-1){
					u=pu_gauge_f(index,1);
					_suNf_mul(*u,_update_par.SF_ds,*u);
				}
				if(iy!=Y_EXT-1){
					u=pu_gauge_f(index,2);
					_suNf_mul(*u,_update_par.SF_ds,*u);
				}
				if(iz!=Z_EXT-1){
					u=pu_gauge_f(index,3);
					_suNf_mul(*u,_update_par.SF_ds,*u);
				}
			}
		}
	}
#endif

#ifdef ANTIPERIODIC_BC_X
	if(COORD[1]==0) {
		int index;
		int it,iy,iz;
		suNf *u;
		for (it=0;it<T_EXT;++it)
		for (iy=0;iy<Y_EXT;++iy)
		for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(it,2*X_BORDER,iy,iz);
			if(index!=-1) {
				suNf *u=pu_gauge_f(index,1);
				_suNf_minus(*u,*u);
			}
		}
	}
#endif

#ifdef ANTIPERIODIC_BC_Y
	if(COORD[2]==0) {
		int index;
		int ix,it,iz;
		suNf *u;
		for (it=0;it<T_EXT;++it)
		for (ix=0;ix<X_EXT;++ix)
		for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(it,ix,2*Y_BORDER,iz);
			if(index!=-1) {
				suNf *u=pu_gauge_f(index,2);
				_suNf_minus(*u,*u);
			}
		}
	}
#endif

#ifdef ANTIPERIODIC_BC_Z
	if(COORD[3]==0) {
		int index;
		int ix,iy,it;
		suNf *u;
		for (it=0;it<T_EXT;++it)
		for (ix=0;ix<X_EXT;++ix)
		for (iy=0;iy<Y_EXT;++iy){
			index=ipt_ext(it,ix,iy,2*Z_BORDER);
			if(index!=-1) {
				suNf *u=pu_gauge_f(index,3);
				_suNf_minus(*u,*u);
			}
		}
	}
#endif
}


void apply_bc_flt(){
#if defined(ANTIPERIODIC_BC_T) && !defined(ROTATED_SF) && !defined(BASIC_SF)
	if(COORD[0]==0) {
		int index;
		int ix,iy,iz;
		suNf_flt *u;
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(2*T_BORDER,ix,iy,iz);
			if(index!=-1) {
				u=pu_gauge_f_flt(index,0);
				_suNf_minus(*u,*u);
			}
		}
	}
#elif defined(ROTATED_SF)
#warning SCRIVERE UN TEST PER VERIFICARE CHE QUESTO EE CONSISTENTE CON I BORDI
	if(COORD[0] == 0) {
		int index;
		int ix,iy,iz;
		suNf_flt *u;
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T_BORDER+1,ix,iy,iz);
			if(index!=-1) {
				if(ix!=X_EXT-1){
					u=pu_gauge_f_flt(index,1);
					_suNf_mul(*u,((float)_update_par.SF_ds),*u);
				}
				if(iy!=Y_EXT-1){
					u=pu_gauge_f_flt(index,2);
					_suNf_mul(*u,((float)_update_par.SF_ds),*u);
				}
				if(iz!=Z_EXT-1){
					u=pu_gauge_f_flt(index,3);
					_suNf_mul(*u,((float)_update_par.SF_ds),*u);
				}
			}
		}
	}
	if(COORD[0] == NP_T-1) {
		int index;
		int ix,iy,iz;
		suNf_flt *u;
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T+T_BORDER-1,ix,iy,iz);
			if(index!=-1) {
				if(ix!=X_EXT-1){
					u=pu_gauge_f_flt(index,1);
					_suNf_mul(*u,((float)_update_par.SF_ds),*u);
				}
				if(iy!=Y_EXT-1){
					u=pu_gauge_f_flt(index,2);
					_suNf_mul(*u,((float)_update_par.SF_ds),*u);
				}
				if(iz!=Z_EXT-1){
					u=pu_gauge_f_flt(index,3);
					_suNf_mul(*u,((float)_update_par.SF_ds),*u);
				}
			}
		}
	}
#endif

#ifdef ANTIPERIODIC_BC_X
	if(COORD[1]==0) {
		int index;
		int it,iy,iz;
		suNf_flt *u;
		for (it=0;it<T_EXT;++it) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(it,2*X_BORDER,iy,iz);
			if(index!=-1) {
				u=pu_gauge_f_flt(index,1);
				_suNf_minus(*u,*u);
			}
		}
	}
#endif

#ifdef ANTIPERIODIC_BC_Y
	if(COORD[2]==0) {
		int index;
		int ix,it,iz;
		suNf_flt *u;
		for (it=0;it<T_EXT;++it) for (ix=0;ix<X_EXT;++ix) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(it,ix,2*Y_BORDER,iz);
			if(index!=-1) {
				u=pu_gauge_f_flt(index,2);
				_suNf_minus(*u,*u);
			}
		}
	}
#endif

#ifdef ANTIPERIODIC_BC_Z
	if(COORD[3]==0) {
		int index;
		int ix,iy,it;
		suNf_flt *u;
		for (it=0;it<T_EXT;++it) for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy){
			index=ipt_ext(it,ix,iy,2*Z_BORDER);
			if(index!=-1) {
				u=pu_gauge_f_flt(index,3);
				_suNf_minus(*u,*u);
			}
		}
	}
#endif
}





#ifdef BASIC_SF

void SF_spinor_bcs(spinor_field *sp)
{
	int ix,iy,iz,index;

#warning SCRIVERE UN TEST PER VERIFICARE CHE QUESTO EE CONSISTENTE CON I BORDI
	if(COORD[0] == 0) {
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T_BORDER,ix,iy,iz);
			if(index!=-1 && index<sp->type->gsize) {
				_spinor_zero_g(*_FIELD_AT(sp,index));
			}
			index=ipt_ext(T_BORDER+1,ix,iy,iz);
			if(index!=-1 && index<sp->type->gsize) {
				_spinor_zero_g(*_FIELD_AT(sp,index));
			}
		}
	}
	if(COORD[0] == NP_T-1) {
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T+T_BORDER-1,ix,iy,iz);
			if(index!=-1 && index<sp->type->gsize) {
				_spinor_zero_g(*_FIELD_AT(sp,index));
			}
		}
	}

}


double SF_test_spinor_bcs(spinor_field *sp)
{
	_DECLARE_INT_ITERATOR(i);
	int ix, iy, iz;
	double temp = 0;
	double total = 0;
	suNf_spinor *sp_temp;
	start_sf_sendrecv(sp);

	if(COORD[0]==0)
	{
	_PIECE_FOR(sp->type,i)
	{
		_SITE_FOR(sp->type,i)
		{
		       	for (ix=0; ix<GLB_X/NP_X; ++ix)
		        for (iy=0; iy<GLB_Y/NP_Y; ++iy)
		        for (iz=0; iz<GLB_Z/NP_Z; ++iz)
		        {
			{
			{
			if (ipt(0,ix,iy,iz)==i||ipt(1,ix,iy,iz)==i)
			{
			sp_temp=_FIELD_AT(sp,i);
		_vector_prod_re_f(temp,(*sp_temp).c[0],(*sp_temp).c[0]);
		total+=temp;
		_vector_prod_re_f(temp,(*sp_temp).c[1],(*sp_temp).c[1]);
		total+=temp;
		_vector_prod_re_f(temp,(*sp_temp).c[2],(*sp_temp).c[2]);
		total+=temp;
		_vector_prod_re_f(temp,(*sp_temp).c[3],(*sp_temp).c[3]);
		total+=temp;

		_vector_prod_im_f(temp,(*sp_temp).c[0],(*sp_temp).c[0]);
		total+=temp;
		_vector_prod_im_f(temp,(*sp_temp).c[1],(*sp_temp).c[1]);
		total+=temp;
		_vector_prod_im_f(temp,(*sp_temp).c[2],(*sp_temp).c[2]);
		total+=temp;
		_vector_prod_im_f(temp,(*sp_temp).c[3],(*sp_temp).c[3]);
		total+=temp;
			}
			}
			}
			}
		}
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(sp);
		}
	}
	}

	if(COORD[0]==NP_T-1)
	{
	_PIECE_FOR(sp->type,i)
	{
		_SITE_FOR(sp->type,i)
		{
		       	for (ix=0; ix<GLB_X/NP_X; ++ix)
		        for (iy=0; iy<GLB_Y/NP_Y; ++iy)
		        for (iz=0; iz<GLB_Z/NP_Z; ++iz)
		        {
			{
			{
			if (ipt((GLB_T/NP_T)-1,ix,iy,iz)==i)
			{
			sp_temp=_FIELD_AT(sp,i);
		_vector_prod_re_f(temp,(*sp_temp).c[0],(*sp_temp).c[0]);
		total+=temp;
		_vector_prod_re_f(temp,(*sp_temp).c[1],(*sp_temp).c[1]);
		total+=temp;
		_vector_prod_re_f(temp,(*sp_temp).c[2],(*sp_temp).c[2]);
		total+=temp;
		_vector_prod_re_f(temp,(*sp_temp).c[3],(*sp_temp).c[3]);
		total+=temp;

		_vector_prod_im_f(temp,(*sp_temp).c[0],(*sp_temp).c[0]);
		total+=temp;
		_vector_prod_im_f(temp,(*sp_temp).c[1],(*sp_temp).c[1]);
		total+=temp;
		_vector_prod_im_f(temp,(*sp_temp).c[2],(*sp_temp).c[2]);
		total+=temp;
		_vector_prod_im_f(temp,(*sp_temp).c[3],(*sp_temp).c[3]);
		total+=temp;
			}
			}
			}
			}
		}
		if(_PIECE_INDEX(i)==0)
		{
			complete_sf_sendrecv(sp);
		}
	}
	}

	return total;

}

#endif /* BASIC_SF */





#if defined(BASIC_SF) || defined(ROTATED_SF)

void SF_force_bcs(suNg_av_field *force) {
	int ix,iy,iz,index;

#warning SCRIVERE UN TEST PER VERIFICARE CHE QUESTO EE CONSISTENTE CON I BORDI
	if(COORD[0] == 0) {
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T_BORDER,ix,iy,iz);
			if(index!=-1 && index<force->type->gsize) {
				_algebra_vector_zero_g(*_4FIELD_AT(force,index,0));
				if(ix!=X_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,1));
				}
				if(iy!=Y_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,2));
				}
				if(iz!=Z_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,3));
				}
			}
			index=ipt_ext(T_BORDER+1,ix,iy,iz);
			if(index!=-1 && index<force->type->gsize) {
				if(ix!=X_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,1));
				}
				if(iy!=Y_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,2));
				}
				if(iz!=Z_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,3));
				}
			}
		}
	}
	if(COORD[0] == NP_T-1) {
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T+T_BORDER-1,ix,iy,iz);
			if(index!=-1 && index<force->type->gsize) {
				if(ix!=X_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,1));
				}
				if(iy!=Y_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,2));
				}
				if(iz!=Z_EXT-1){
					_algebra_vector_zero_g(*_4FIELD_AT(force,index,3));
				}
			}
		}
	}

}


#endif /* BASIC_SF || ROTATED_SF */



#define PI 3.141592653589793238462643383279502884197
#define ST 1.414213562373095048801688724209698078570


#if NG==2

double SF_eta = PI/4.0;
double SF_phi0_dn[NG] = {0., 0.};
double SF_phi1_dn[NG] = {-1., 1.};
double SF_phi0_up[NG] = {-PI, PI};
double SF_phi1_up[NG] = {1., -1.};

#elif NG==3

double SF_eta = 0.;
double SF_phi0_dn[NG] = {-PI/3., 0., PI/3.};
double SF_phi1_dn[NG] = {1., -.5, -.5};
double SF_phi0_up[NG] = {-PI, PI/3., 2.*PI/3.};
double SF_phi1_up[NG] = {-1., .5, .5};

#elif NG==4

double SF_eta = 0.;
double SF_phi0_dn[NG] = {-ST*PI/4., ST*PI/4.-PI/2., PI/2.-ST*PI/4., ST*PI/4.};
double SF_phi1_dn[NG] = {-.5, -.5, .5, .5};
double SF_phi0_up[NG] = {-ST*PI/4.-PI/2., -PI+ST*PI/4., PI-ST*PI/4., PI/2.+ST*PI/4.};
double SF_phi1_up[NG] = {.5, .5, -.5, -.5};

#else

#error SF boundary conditions not defined at this NG

#endif




void SF_gauge_bcs(suNg_field *gf, int strength)
{
  int index;
  int ix, iy, iz;
  int k;
  suNg *u;

  error(gf==NULL,1,"SF_gauge_bcs [random_fields.c]",
	"Attempt to access unallocated memory space");   

  /*Boundary gauge fields*/
	suNg Bound0, BoundT;
	if(strength==1) {  /*SF bcs*/
		_suNg_zero(Bound0);
		for(k=0; k<NG; k++) {
			Bound0.c[(1+NG)*k].re = cos((SF_phi0_dn[k]+SF_phi1_dn[k]*SF_eta)/(GLB_T-2));
			Bound0.c[(1+NG)*k].im = sin((SF_phi0_dn[k]+SF_phi1_dn[k]*SF_eta)/(GLB_T-2));
		}
		_suNg_zero(BoundT);
		for(k=0; k<NG; k++) {
			BoundT.c[(1+NG)*k].re = cos((SF_phi0_up[k]+SF_phi1_up[k]*SF_eta)/(GLB_T-2));
			BoundT.c[(1+NG)*k].im = sin((SF_phi0_up[k]+SF_phi1_up[k]*SF_eta)/(GLB_T-2));
		}
	} else { /*UNIT bcs*/
		_suNg_unit(Bound0);
		_suNg_unit(BoundT);
	}	
  
	
#warning SCRIVERE UN TEST PER VERIFICARE CHE QUESTO EE CONSISTENTE CON I BORDI
	if(COORD[0] == 0) {
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			if(T_BORDER>0) {
				index=ipt_ext(T_BORDER-1,ix,iy,iz);
				if(index!=-1) {
					u=pu_gauge(index,0);
					_suNg_unit(*u);
				}
			}
			index=ipt_ext(T_BORDER,ix,iy,iz);
			if(index!=-1) {
				u=pu_gauge(index,0);
				_suNg_unit(*u);
				if(ix!=X_EXT-1){
					u=pu_gauge(index,1);
					_suNg_unit(*u);
				}
				if(iy!=Y_EXT-1){
					u=pu_gauge(index,2);
					_suNg_unit(*u);
				}
				if(iz!=Z_EXT-1){
					u=pu_gauge(index,3);
					_suNg_unit(*u);
				}
			}
			index=ipt_ext(T_BORDER+1,ix,iy,iz);
			if(index!=-1) {
				if(ix!=X_EXT-1){
					u=pu_gauge(index,1);
					*u = Bound0;
				}
				if(iy!=Y_EXT-1){
					u=pu_gauge(index,2);
					*u = Bound0;
				}
				if(iz!=Z_EXT-1){
					u=pu_gauge(index,3);
					*u = Bound0;
				}
			}
		}
	}
	if(COORD[0] == NP_T-1) {
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			index=ipt_ext(T+T_BORDER-1,ix,iy,iz);
			if(index!=-1) {
				u=pu_gauge(index,0);
				_suNg_unit(*u);
				if(ix!=X_EXT-1){
					u=pu_gauge(index,1);
					*u = BoundT;
				}
				if(iy!=Y_EXT-1){
					u=pu_gauge(index,2);
					*u = BoundT;
				}
				if(iz!=Z_EXT-1){
					u=pu_gauge(index,3);
					*u = BoundT;
				}
			}
		}
	}
	
}

