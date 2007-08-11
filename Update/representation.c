#include "global.h"
#include "representation.h"
#include "utils.h"
#include <math.h>


void _group_represent2(suNf* v, suNg *u) {

#define XG(m,a,b) ((m)+(a)*NG+(b))
#define XF(m,a,b) ((m)+(a)*NF+(b))


#ifdef REPR_ADJOINT

	int A, C;
	int a, b, i, j, k, c, d;
	float* vf = (float*)v;
	complex* uf = (complex*)u;
	
	suNg m;
	complex* mf = (complex*)(&m);
	
	
	A = 0;
	for(a = 0; a < NG; a++) for(b = (a==0)?1:0; b < NG; b++) {
		if(a > b)
		{
		  for(i = 0; i < NG; i++) for(j = i; j < NG; j++) {
		    XG(mf,i,j)->re = XG(uf,i,a)->re*XG(uf,j,b)->re+XG(uf,i,a)->im*XG(uf,j,b)->im+XG(uf,i,b)->re*XG(uf,j,a)->re+XG(uf,i,b)->im*XG(uf,j,a)->im;
		    XG(mf,i,j)->im = -XG(uf,i,a)->re*XG(uf,j,b)->im+XG(uf,i,a)->im*XG(uf,j,b)->re-XG(uf,i,b)->re*XG(uf,j,a)->im+XG(uf,i,b)->im*XG(uf,j,a)->re;
		  }
		}
		else if(a < b)
		{
		  for(i = 0; i < NG; i++) for(j = i; j < NG; j++) {
		    XG(mf,i,j)->im = XG(uf,i,a)->re*XG(uf,j,b)->re+XG(uf,i,a)->im*XG(uf,j,b)->im-XG(uf,i,b)->re*XG(uf,j,a)->re-XG(uf,i,b)->im*XG(uf,j,a)->im;
		    XG(mf,i,j)->re = +XG(uf,i,a)->re*XG(uf,j,b)->im-XG(uf,i,a)->im*XG(uf,j,b)->re-XG(uf,i,b)->re*XG(uf,j,a)->im+XG(uf,i,b)->im*XG(uf,j,a)->re;
		  }
		}
		else if(a == b)
		{
		  for(i = 0; i < NG; i++) for(j = i; j < NG; j++) {
		  	XG(mf,i,j)->re = -a*(XG(uf,i,a)->re*XG(uf,j,a)->re+XG(uf,i,a)->im*XG(uf,j,a)->im);
		  	XG(mf,i,j)->im = -a*(XG(uf,i,a)->im*XG(uf,j,a)->re-XG(uf,i,a)->re*XG(uf,j,a)->im);
		  	for(k = 0; k < a; k++) {
			  	XG(mf,i,j)->re += XG(uf,i,k)->re*XG(uf,j,k)->re+XG(uf,i,k)->im*XG(uf,j,k)->im;
			  	XG(mf,i,j)->im += XG(uf,i,k)->im*XG(uf,j,k)->re-XG(uf,i,k)->re*XG(uf,j,k)->im;
			  }
			  XG(mf,i,j)->re *= sqrt(2./(a*(a+1.)));
			  XG(mf,i,j)->im *= sqrt(2./(a*(a+1.)));
		  }
		}

		C = 0;
		for(c = 0; c < NG; c++) for(d = (c==0)?1:0; d < NG; d++) {
			if(c > d)
			{
				*(XF(vf,C,A)) = XG(mf,d,c)->re;
			}
			else if(c < d)
			{
				*(XF(vf,C,A)) = XG(mf,c,d)->im;
			}
			else if(c == d)
			{
				*(XF(vf,C,A)) = -c*XG(mf,c,c)->re;
		  	for(k = 0; k < c; k++) {
					*(XF(vf,C,A)) += XG(mf,k,k)->re;
		  	}
		  	*(XF(vf,C,A)) *= sqrt(.5/(c*(c+1.)));
			}
			
			C++;
		}
			
		A++;
	}

#elif defined REPR_SYMMETRIC

	const float st = sqrt(2.);
	int A, C;
	int a, b, i, j, c, d;
	complex* vf = (complex*)v;
	complex* uf = (complex*)u;
	
	suNg m;
	complex* mf = (complex*)(&m);

	A = 0;
	for(a = 0; a < NG; a++) {
		for(b = 0; b < a; b++) {
		  for(i = 0; i < NG; i++) for(j = i; j < NG; j++) {
		    XG(mf,i,j)->re = XG(uf,i,a)->re*XG(uf,j,b)->re-XG(uf,i,a)->im*XG(uf,j,b)->im+XG(uf,i,b)->re*XG(uf,j,a)->re-XG(uf,i,b)->im*XG(uf,j,a)->im;
		    XG(mf,i,j)->im = XG(uf,i,a)->re*XG(uf,j,b)->im+XG(uf,i,a)->im*XG(uf,j,b)->re+XG(uf,i,b)->re*XG(uf,j,a)->im+XG(uf,i,b)->im*XG(uf,j,a)->re;
		  }

			C = 0;
			for(c = 0; c < NG; c++) {
				for(d = 0; d < c; d++) {
					XF(vf,C,A)->re = XG(mf,d,c)->re;
					XF(vf,C,A)->im = XG(mf,d,c)->im;
					C++;
				}
				XF(vf,C,A)->re = XG(mf,c,c)->re/st;
				XF(vf,C,A)->im = XG(mf,c,c)->im/st;
				C++;
			}

			A++;
		}
			
	  for(i = 0; i < NG; i++) for(j = i; j < NG; j++) {
	    XG(mf,i,j)->re = XG(uf,i,a)->re*XG(uf,j,a)->re-XG(uf,i,a)->im*XG(uf,j,a)->im;
	    XG(mf,i,j)->im = XG(uf,i,a)->re*XG(uf,j,a)->im+XG(uf,i,a)->im*XG(uf,j,a)->re;
	  }

		C = 0;
		for(c = 0; c < NG; c++) {
			for(d = 0; d < c; d++) {
				XF(vf,C,A)->re = XG(mf,d,c)->re*st;
				XF(vf,C,A)->im = XG(mf,d,c)->im*st;
				C++;
			}
			XF(vf,C,A)->re = XG(mf,c,c)->re;
			XF(vf,C,A)->im = XG(mf,c,c)->im;
			C++;
		}

		A++;
	}

#elif defined REPR_ANTISYMMETRIC

	int A, C;
	int a, b, i, j, c, d;
	complex* vf = (complex*)v;
	complex* uf = (complex*)u;
	
	suNg m;
	complex* mf = (complex*)(&m);

	A = 0;
	for(a = 1; a < NG; a++) for(b = 0; b < a; b++)
	{
	  for(i = 0; i < NG; i++) for(j = i; j < NG; j++) {
	    XG(mf,i,j)->re = XG(uf,i,a)->re*XG(uf,j,b)->re-XG(uf,i,a)->im*XG(uf,j,b)->im-XG(uf,i,b)->re*XG(uf,j,a)->re+XG(uf,i,b)->im*XG(uf,j,a)->im;
	    XG(mf,i,j)->im = XG(uf,i,a)->re*XG(uf,j,b)->im+XG(uf,i,a)->im*XG(uf,j,b)->re-XG(uf,i,b)->re*XG(uf,j,a)->im-XG(uf,i,b)->im*XG(uf,j,a)->re;
	  }

		C = 0;
		for(c = 1; c < NG; c++) for(d = 0; d < c; d++) {
			XF(vf,C,A)->re = -XG(mf,d,c)->re;
			XF(vf,C,A)->im = -XG(mf,d,c)->im;
			C++;
		}

		A++;
	}

#elif defined REPR_FUNDAMENTAL

	*v = *((suNf *)u); 
#endif

#undef XG
#undef XF

}




void represent_gauge_field() {
#ifndef REPR_FUNDAMENTAL
   int i;
   suNf *Ru=u_gauge_f;
   suNg *u=u_gauge;

   for (i=0; i<VOLUME*4; ++i){
      /*_group_represent2(Ru,u); */
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

void represent_gauge_field_dble() {
#ifndef REPR_FUNDAMENTAL
   int i;
   suNf_dble *Ru=u_gauge_f_dble;
   suNg_dble *u=u_gauge_dble;

   for (i=0; i<VOLUME*4; ++i){
      /*_group_represent2(Ru,u); */
      _group_represent(*Ru,*u);
      ++Ru;
      ++u;
   }
	 apply_bc_dble();
#else
	static short int first_time=1;
	 if(first_time) {
		 first_time=0;
	   u_gauge_f_dble=(suNf_dble *)((void*)u_gauge_dble);
		 apply_bc_dble();
	 }
#endif
}



