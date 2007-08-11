#include "global.h"
#include "utils.h"
#include "suN.h"

void apply_bc(){
#ifdef ANTIPERIODIC_BC_T
	{
		int ix,iy,iz;
		for (ix=0;ix<L;++ix){
			for (iy=0;iy<L;++iy){
				for (iz=0;iz<L;++iz){
					suNf *u=pu_gauge_f(ipt[T-1][ix][iy][iz],0);
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
			for (iy=0;iy<L;++iy){
				for (iz=0;iz<L;++iz){
					suNf *u=pu_gauge_f(ipt[it][L-1][iy][iz],1);
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
			for (ix=0;ix<L;++ix){
				for (iz=0;iz<L;++iz){
					suNf *u=pu_gauge_f(ipt[it][ix][L-1][iz],2);
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
			for (ix=0;ix<L;++ix){
				for (iy=0;iy<L;++iy){
					suNf *u=pu_gauge_f(ipt[it][ix][iy][L-1],3);
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
		for (ix=0;ix<L;++ix){
			for (iy=0;iy<L;++iy){
				for (iz=0;iz<L;++iz){
					suNf_flt *u=pu_gauge_f_flt(ipt[T-1][ix][iy][iz],0);
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
			for (iy=0;iy<L;++iy){
				for (iz=0;iz<L;++iz){
					suNf_flt *u=pu_gauge_f_flt(ipt[it][L-1][iy][iz],1);
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
			for (ix=0;ix<L;++ix){
				for (iz=0;iz<L;++iz){
					suNf_flt *u=pu_gauge_f_flt(ipt[it][ix][L-1][iz],2);
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
			for (ix=0;ix<L;++ix){
				for (iy=0;iy<L;++iy){
					suNf_flt *u=pu_gauge_f_flt(ipt[it][ix][iy][L-1],3);
					_suNf_minus(*u,*u);
				}
			}
		}
	}
#endif
}
