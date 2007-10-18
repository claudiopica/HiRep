#include "global.h"
#include "utils.h"
#include "suN.h"

void apply_bc(){
#ifdef ANTIPERIODIC_BC_T
	{
		int ix,iy,iz;
		for (ix=0;ix<X;++ix){
			for (iy=0;iy<Y;++iy){
				for (iz=0;iz<Z;++iz){
					suNf *u=pu_gauge_f(ipt(T-1,ix,iy,iz),0);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_X
	{
		int it,iy,iz;
		for (it=0;it<T;++it){
			for (iy=0;iy<Y;++iy){
				for (iz=0;iz<Z;++iz){
					suNf *u=pu_gauge_f(ipt(it,X-1,iy,iz),1);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_Y
	{
		int it,ix,iz;
		for (it=0;it<T;++it){
			for (ix=0;ix<X;++ix){
				for (iz=0;iz<Z;++iz){
					suNf *u=pu_gauge_f(ipt(it,ix,Y-1,iz),2);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_Z
	{
		int it,ix,iy;
		for (it=0;it<T;++it){
			for (ix=0;ix<X;++ix){
				for (iy=0;iy<Y;++iy){
					suNf *u=pu_gauge_f(ipt(it,ix,iy,Z-1),3);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
}

void apply_bc_flt(){
#ifdef ANTIPERIODIC_BC_T
	{
		int ix,iy,iz;
		for (ix=0;ix<X;++ix){
			for (iy=0;iy<Y;++iy){
				for (iz=0;iz<Z;++iz){
					suNf_flt *u=pu_gauge_f_flt(ipt(T-1,ix,iy,iz),0);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_X
	{
		int it,iy,iz;
		for (it=0;it<T;++it){
			for (iy=0;iy<Y;++iy){
				for (iz=0;iz<Z;++iz){
					suNf_flt *u=pu_gauge_f_flt(ipt(it,X-1,iy,iz),1);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_Y
	{
		int it,ix,iz;
		for (it=0;it<T;++it){
			for (ix=0;ix<X;++ix){
				for (iz=0;iz<Z;++iz){
					suNf_flt *u=pu_gauge_f_flt(ipt(it,ix,Y-1,iz),2);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_Z
	{
		int it,ix,iy;
		for (it=0;it<T;++it){
			for (ix=0;ix<X;++ix){
				for (iy=0;iy<Y;++iy){
					suNf_flt *u=pu_gauge_f_flt(ipt(it,ix,iy,Z-1),3);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
}
