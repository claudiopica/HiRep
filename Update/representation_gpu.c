/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Agostino Patella, Ari Hietanen          *   
* All rights reserved.                                                      * 
\***************************************************************************/
#ifdef WITH_GPU
#include "global.h"
#include "representation.h"
#include "utils.h"
#include <math.h>
#include "gpu.h"
#include "suN.h"

#define XG(m,a,b) ((m)+(a)*NG+(b))
#define XF(m,a,b) ((m)+(a)*NF+(b))

__global__ void group_represent_kernel(suNf* gauge_f, suNg *gauge, int direction, int N) {
  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int vol4h = N/2;
  if (ix>=N/2){
    ix -= vol4h;
    gauge += 4*vol4h;
    gauge_f += 4*vol4h;
  }
  suNg u;  
  suNf v;  
  _suNg_read_gpu(vol4h,u,gauge,ix,direction);
  do {
#ifdef REPR_ADJOINT
    int A, C;
    int a, b, i, j, k, c, d;
    double* vf = (double*)(&v);
    complex* uf = (complex*)(&u);

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

    const double st = sqrt(2.);
    int A, C;
    int a, b, i, j, c, d;
    complex* vf = (complex*)(&v);
    complex* uf = (complex*)(&u);

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
    complex* vf = (complex*)(&v);
    complex* uf = (complex*)(&u);

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

    v = *((suNf*) (&u));
#endif
  } while(0);
  _suNf_write_gpu(vol4h,v,gauge_f,ix,direction);
}


__global__ void group_represent_kernel_flt(suNf_flt* gauge_f, suNg_flt *gauge,int direction,int N) {

  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int vol4h = N/2;
  if (ix>=N/2){
    ix -= vol4h;
    gauge += 4*vol4h;
    gauge_f += 4*vol4h;
  }
  suNg_flt u;  
  suNf_flt v;  
  _suNg_read_gpu(vol4h,u,gauge,ix,direction);
  do {
#ifdef REPR_ADJOINT

    int A, C;
    int a, b, i, j, k, c, d;
    float* vf = (float*)(&v);
    complex_flt* uf = (complex_flt*)(&u);

    suNg_flt m;
    complex_flt* mf = (complex_flt*)(&m);


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

    const double st = sqrt(2.);
    int A, C;
    int a, b, i, j, c, d;
    complex_flt* vf = (complex_flt*)(&v);
    complex_flt* uf = (complex_flt*)(&u);

    suNg_flt m;
    complex_flt* mf = (complex_flt*)(&m);

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
    complex_flt* vf = (complex_flt*)(&v);
    complex_flt* uf = (complex_flt*)(&u);

    suNg_flt m;
    complex_flt* mf = (complex_flt*)(&m);

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

#endif
  }while (0);
  _suNf_flt_write_gpu(vol4h,v,gauge_f,ix,direction);
}

#undef XG
#undef XF


#include "communications.h"

void represent_gauge_field() {
#if (defined(ALLOCATE_REPR_GAUGE_FIELD))
  int N = T*X*Y*Z;
  int grid = (N-1)/BLOCK_SIZE + 1;
  int d;
  for (d=0;d<4;++d){
    group_represent_kernel<<<grid,BLOCK_SIZE>>>(u_gauge_f->gpu_ptr,u_gauge->gpu_ptr,d,N);
  }
  /* loop on local lattice first */
  apply_bc();
#else
  static int first_time=1;
   /* wait gauge field transfer */
  complete_gf_sendrecv(u_gauge);

  if(first_time) {
    first_time=0;
    u_gauge_f=(suNf_field *)((void*)u_gauge);
    apply_bc();
  }
#endif
}

void represent_gauge_field_flt() {
#if (defined(ALLOCATE_REPR_GAUGE_FIELD))
  int N = T*X*Y*Z;
  int grid = (N-1)/BLOCK_SIZE + 1;
  int d;
  for (d=0;d<4;++d){
    group_represent_kernel_flt<<<grid,BLOCK_SIZE>>>(u_gauge_f_flt->gpu_ptr,u_gauge_flt->gpu_ptr,d,N);
  }  
  apply_bc_flt();
#else
  static int first_time=1;
  /* wait gauge field transfer */
  complete_gf_sendrecv_flt(u_gauge_flt);
  if(first_time) {
    first_time=0;
    u_gauge_f_flt=(suNf_field_flt *)((void*)u_gauge_flt);
    apply_bc_flt();
  }
#endif
}
#endif


