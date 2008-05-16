#include "global.h"
#include "utils.h"
#include "suN.h"

void apply_bc(){
#ifdef ANTIPERIODIC_BC_T
	if(proc_t(myid)==0) {
		int index;
		int ix,iy,iz;
		for (ix=0;ix<X_EXT;++ix){
			for (iy=0;iy<Y_EXT;++iy){
				for (iz=0;iz<Z_EXT;++iz){
					index=ipt_ext(2*T_BORDER,ix,iy,iz);
					if(index!=-1) {
						suNf *u=pu_gauge_f(index,0);
						_suNf_minus(*u,*u);
					}
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_X
	if(proc_x(myid)==0) {
		int index;
		int it,iy,iz;
		for (it=0;it<T_EXT;++it){
			for (iy=0;iy<Y_EXT;++iy){
				for (iz=0;iz<Z_EXT;++iz){
					index=ipt_ext(it,2*X_BORDER,iy,iz);
					if(index!=-1) {
						suNf *u=pu_gauge_f(index,1);
						_suNf_minus(*u,*u);
					}
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_Y
	if(proc_y(myid)==0) {
		int index;
		int it,ix,iz;
		for (it=0;it<T_EXT;++it){
			for (ix=0;ix<X_EXT;++ix){
				for (iz=0;iz<Z_EXT;++iz){
					index=ipt_ext(it,ix,2*Y_BORDER,iz);
					if(index!=-1) {
						suNf *u=pu_gauge_f(index,2);
						_suNf_minus(*u,*u);
					}
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_Z
	if(proc_z(myid)==0) {
		int index;
		int it,ix,iy;
		for (it=0;it<T_EXT;++it){
			for (ix=0;ix<X_EXT;++ix){
				for (iy=0;iy<Y_EXT;++iy){
					index=ipt_ext(it,ix,iy,2*Z_BORDER);
					if(index!=-1) {
						suNf *u=pu_gauge_f(index,3);
						_suNf_minus(*u,*u);
					}
				}
			}
		}
	}
#endif
}

void apply_bc_flt(){
#ifdef ANTIPERIODIC_BC_T
	if(proc_t(myid)==0) {
		int index;
		int ix,iy,iz;
		for (ix=0;ix<X_EXT;++ix){
			for (iy=0;iy<Y_EXT;++iy){
				for (iz=0;iz<Z_EXT;++iz){
					index=ipt_ext(2*T_BORDER,ix,iy,iz);
					if(index!=-1) {
						suNf_flt *u=pu_gauge_f_flt(index,0);
						_suNf_minus(*u,*u);
					}
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_X
	if(proc_x(myid)==0) {
		int index;
		int it,iy,iz;
		for (it=0;it<T_EXT;++it){
			for (iy=0;iy<Y_EXT;++iy){
				for (iz=0;iz<Z_EXT;++iz){
					index=ipt_ext(it,2*X_BORDER,iy,iz);
					if(index!=-1) {
						suNf_flt *u=pu_gauge_f_flt(index,1);
						_suNf_minus(*u,*u);
					}
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_Y
	if(proc_y(myid)==0) {
		int index;
		int it,ix,iz;
		for (it=0;it<T_EXT;++it){
			for (ix=0;ix<X_EXT;++ix){
				for (iz=0;iz<Z_EXT;++iz){
					index=ipt_ext(it,ix,2*Y_BORDER,iz);
					if(index!=-1) {
						suNf_flt *u=pu_gauge_f_flt(index,2);
						_suNf_minus(*u,*u);
					}
				}
			}
		}
	}
#endif
#ifdef ANTIPERIODIC_BC_Z
	if(proc_z(myid)==0) {
		int index;
		int it,ix,iy;
		for (it=0;it<T_EXT;++it){
			for (ix=0;ix<X_EXT;++ix){
				for (iy=0;iy<Y_EXT;++iy){
					index=ipt_ext(it,ix,iy,2*Z_BORDER);
					if(index!=-1) {
						suNf_flt *u=pu_gauge_f_flt(index,3);
						_suNf_minus(*u,*u);
					}
				}
			}
		}
	}
#endif
}
