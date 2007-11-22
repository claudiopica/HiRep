/*******************************************************************************
*
* Gauge covariance of the Dirac operator
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

static double hmass=0.1;
static suNg *g;


static void D(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

static void random_g(void)
{
   int ix;

   for (ix=0;ix<VOLUME;ix++)
      random_suNg(g+ix);
}

static void transform_u(void)
{
   int ix,iy,mu;
   suNg *u,v;

   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         iy=iup(ix,mu);
         u=pu_gauge(ix,mu);
         _suNg_times_suNg_dagger(v,*u,g[iy]);
         _suNg_times_suNg(*u,g[ix],v);
      }
   }
   
   represent_gauge_field();
}

static void transform_s(spinor_field *out, spinor_field *in)
{
   int ix;
   suNf gfx;
   suNf_spinor r,s;

   for (ix=0;ix<VOLUME;ix++)
   {
      s = *_SPINOR_AT_SITE(in,ix);
      
      _group_represent2(&gfx,&(g[ix]));

      _suNf_multiply(r.c[0],gfx,s.c[0]);
      _suNf_multiply(r.c[1],gfx,s.c[1]);
      _suNf_multiply(r.c[2],gfx,s.c[2]);
      _suNf_multiply(r.c[3],gfx,s.c[3]);

      *_SPINOR_AT_SITE(out,ix) = r;
   }   
}


int main(int argc,char *argv[])
{
   double sig,tau;
   spinor_field *s0,*s1,*s2,*s3;

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   geometry_eo_lexi();
   printf("The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
   printf("\n");
   
   level=0;
   seed=123;
   rlxs_init(level,seed);
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
 
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   printf("allocating gfield_f\n");
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();

   /* allocate memory */
	 g=malloc(sizeof(suNg)*VOLUME);
	 set_spinor_len(VOLUME);
	 s0=alloc_spinor_field_f(4);
	 s1=s0+1;
	 s2=s1+1;
	 s3=s2+1;

   printf("Generating a random gauge field... ");
   fflush(stdout);
   random_u();
   represent_gauge_field();
   printf("done.\n");

   set_spinor_len(VOLUME);   
   gaussian_spinor_field(&(s0[0]));
   tau = 1./sqrt(spinor_field_sqnorm_f(s0));
   spinor_field_mul_f(s0,tau,s0);
   
   printf("Generating a random gauge transf... ");
   fflush(stdout);
   random_g();
   printf("done.\n");

   printf("Gauge covariance of the Dirac operator:\n");
   
   D(s1,s0);
   transform_s(s2,s1);
   
   transform_s(s3,s0);
   transform_u();

   spinor_field_zero_f(s1);
   D(s1,s3);
   
   spinor_field_mul_add_assign_f(s1,-1.0,s2);
   sig=spinor_field_sqnorm_f(s1);

   printf("Maximal normalized difference = %.2e\n",sqrt(sig));
   printf("(should be around 1*10^(-15) or so)\n\n");
   
   free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif
	 free_field(s0);

	 free(g);
   exit(0);
}
