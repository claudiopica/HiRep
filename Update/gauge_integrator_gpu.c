/***************************************************************************\
 * Copyright (c) 2012, Claudio Pica                                        *   
 * All rights reserved.                                                    * 
 \**************************************************************************/

#ifdef WITH_GPU

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "gpu.h"
#include "utils.h"
#include "suN.h"
#include <math.h>
#include <assert.h>

/* These YtoU and su2_rotate are copied from Utils/suN_exp.c.tmpl */


/*
 *  U = (1+i*y.sigma/4)(1-i*y.sigma/4)^{-1}
 *  U = u[0] + i*u.sigma/4
 */
__device__ void YtoU(double* u, double *y)
{
    double y2 = y[0]*y[0] + y[1]*y[1] +y[2]*y[2];
    double detY = 1.0 + y2/16.;
    u[0] = (1.0 - y2/16.)/detY;
    u[1] = y[0]/(2.*detY);
    u[2] = y[1]/(2.*detY);
    u[3] = y[2]/(2.*detY);
}


/*
 *  Applica la rotazione di SU(2) definita da
 *  U = s[0] + i*s.sigma/4
 *  al vettore (v1,v2)
 */
__device__ void su2_rotate(double *s,complex *v1,complex *v2)
{
    complex z1, z2;
    z1.re=s[0]*(*v1).re-s[1]*(*v2).im+s[2]*(*v2).re-s[3]*(*v1).im;
    z1.im=s[0]*(*v1).im+s[1]*(*v2).re+s[2]*(*v2).im+s[3]*(*v1).re;
    z2.re=s[0]*(*v2).re-s[1]*(*v1).im-s[2]*(*v1).re+s[3]*(*v2).im;
    z2.im=s[0]*(*v2).im+s[1]*(*v1).re-s[2]*(*v1).im-s[3]*(*v2).re;
    (*v1)=z1;
    (*v2)=z2;
}


#define XG(m,a,b) ((m)+(a)*NG+(b))
#define XA(h,a) ((h)+(a))
#define TI(a,b) ((a)*NG+(b)-1)


__device__ void ExpX_gpu(double dt, suNg_algebra_vector *h, suNg *u)
{
#ifdef WITH_QUATERNIONS
    
    suNg v_tmp, u_tmp;
    
    u_tmp=*u;
    _suNg_exp(dt,*h,v_tmp);
    _suNg_times_suNg(*u,v_tmp,u_tmp);
    
#else //WITH_QUATERNIONS 
    
    int i, j, k, n;
    double y[3];
    double d[NG];
    double s[(NG*(NG-1))/2][4];
    double *hf = (double*)h;
    complex *uf = (complex*)u;
    
    d[0] = 0.0; k = NG;
    for(i = 1; i < NG; i++) {
        double tmp = sqrt( 2./(i*(i+1)) ) * (*XA(hf,k)) / (double)NG;
        d[i] = -i *tmp;
        for(j = 0; j < i; j++)
            d[j] += tmp;
        k+=NG+1;
    }
    
    k = 0; dt*=.5;
    for(j = 1; j<NG; j++)
        for(i = 0; i < j; i++) {
            y[0] =  dt * (*XA(hf,TI(j,i)));
            y[1] = -dt * (*XA(hf,TI(i,j)));
            y[2] =  dt * (d[i]-d[j]);
            YtoU(s[k],y);
            for(n = 0; n < NG; n++)
                su2_rotate(s[k],XG(uf,i,n),XG(uf,j,n));
            k++;
        }
    
    k = NG*(NG-1)/2 - 1;
    for(j = NG-1; j >= 1; j--)
        for(i = j-1; i >= 0; i--) {
            for(n = 0; n < NG; n++)
                su2_rotate(s[k],XG(uf,i,n),XG(uf,j,n));
            k--;
        }
    
#endif //WITH_QUATERNIONS
}

#undef XA
#undef XG
#undef TI

__global__ void gauge_integrator_gpu_kernel(suNg* gauge, suNg_algebra_vector* momenta, double dt, int N)
{   
  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int i;
  int vol4h = N/2;
  suNg u;
  suNg_algebra_vector h;
  ix = min(ix,N-1);
  if (ix>=N/2){
    ix -= vol4h;
    gauge += 4*vol4h;
    momenta += 4*vol4h;
  }
  for (i=0;i<4;++i){
    _suNg_av_read_gpu(vol4h,h,momenta,ix,i);
    _suNg_read_gpu(vol4h,u,gauge,ix,i);
      ExpX_gpu(dt,&h,&u);
    _suNg_write_gpu(vol4h,u,gauge,ix,i);
  }
}


void gauge_integrator_gpu(suNg_av_field *momenta, double dt){
  int N = T*X*Y*Z;//u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
  int grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
  gauge_integrator_gpu_kernel<<<grid,BLOCK_SIZE>>>(u_gauge->gpu_ptr,momenta->gpu_ptr,dt,N);
}

#endif
