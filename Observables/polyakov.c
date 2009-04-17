#include "global.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "logger.h"

static void polyakov_1P(complex* poly, int tx[3], int mu) {
   int i, x[4];
   int glb[4]={GLB_T,GLB_X,GLB_Y,GLB_Z};
   suNg omega;
   suNg omega2;

   error(mu<0 || mu>3,1,"polyakov_1P [polyakov.c]","Bad value for direction parameter.");
   
   x[(mu+1)%4]=tx[0];
   x[(mu+2)%4]=tx[1];
   x[(mu+3)%4]=tx[2];
   
   _suNg_unit(omega);
   for(x[mu]=0; x[mu]<glb[mu]; ) {
      i=ipt(x[0],x[1],x[2],x[3]);
      _suNg_times_suNg(omega2,omega,*pu_gauge(i,mu));
      x[mu]++;
      i=ipt(x[0],x[1],x[2],x[3]);
      _suNg_times_suNg(omega,omega2,*pu_gauge(i,mu));
      x[mu]++;
   }
   
   _suNg_trace_re(poly->re,omega);
   _suNg_trace_im(poly->im,omega);
}

void polyakov(int mu) {
   int tx[3];
   int glb[4]={GLB_T,GLB_X,GLB_Y,GLB_Z};
   complex poly;
   poly.re=poly.im=0.;

   error(mu<0 || mu>3,1,"polyakov [polyakov.c]","Bad value for direction parameter.");

   for(tx[0]=0; tx[0]<glb[(mu+1)%4]; tx[0]++)
   for(tx[1]=0; tx[1]<glb[(mu+2)%4]; tx[1]++)
   for(tx[2]=0; tx[2]<glb[(mu+3)%4]; tx[2]++) {
      complex tmp;
      polyakov_1P(&tmp, tx, mu);
      _complex_add_assign(poly,tmp);
   }
   
   poly.re /= GLB_X*GLB_Y*GLB_Z*GLB_T/glb[mu];
   poly.im /= GLB_X*GLB_Y*GLB_Z*GLB_T/glb[mu];

   lprintf("POLYAKOV",0,"Polyakov direction %d = %1.8e %1.8e\n",mu,poly.re,poly.im);
}

