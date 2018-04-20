/*******************************************************************************
*
* Action of the Dirac operator on plane waves
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
#include "logger.h"
#include "communications.h"

double hmass=0.1;

#ifndef SAFE_MOD_H
#define SAFE_MOD_H

static int safe_mod(int x,int y)
{
   if (x>=0)
      return(x%y);
   else
      return((y-(abs(x)%y))%y);
}

#endif

static suNf_spinor mul_gamma(int mu,suNf_spinor s)
{
   suNf_spinor r;

   if (mu==0)
   {
      _vector_mulc_f(r.c[0],-1.0,s.c[2]);
      _vector_mulc_f(r.c[1],-1.0,s.c[3]);
      _vector_mulc_f(r.c[2],-1.0,s.c[0]);
      _vector_mulc_f(r.c[3],-1.0,s.c[1]);
   }
   else if (mu==1)
   {
      _vector_mulc_f(r.c[0],-I,s.c[3]);
      _vector_mulc_f(r.c[1],-I,s.c[2]);
      _vector_mulc_f(r.c[2],I,s.c[1]);
      _vector_mulc_f(r.c[3],I,s.c[0]);
   }
   else if (mu==2)
   {
      _vector_mulc_f(r.c[0],-1,s.c[3]);
      r.c[1]=s.c[2];
      r.c[2]=s.c[1];
      _vector_mulc_f(r.c[3],-1,s.c[0]);
   }
   else if (mu==3)
   {
      _vector_mulc_f(r.c[0],-I,s.c[2]);
      _vector_mulc_f(r.c[1],I,s.c[3]);
      _vector_mulc_f(r.c[2],I,s.c[0]);
      _vector_mulc_f(r.c[3],-I,s.c[1]);
   }
   else
   {
      r.c[0]=s.c[0];
      r.c[1]=s.c[1];
      _vector_mulc_f(r.c[2],-1.0,s.c[2]);
      _vector_mulc_f(r.c[3],-1.0,s.c[3]);
   }

   return r;
}


int main(int argc,char *argv[])
{
   int i,j,n,ix,mu;
   int bc_t,bc_x,bc_y,bc_z;
   int x0,x1,x2,x3;
   int np[4];
   double ran[4],sp[4];
   double pi,p[4];
   double *rs,r,mp,sig,px;
   double complex z;
   suNf_spinor s,s0,s1,s2,s3;
   spinor_field *ps0,*ps1,*ps2;
   char tmp[256];
   
   /* setup process id and communications */
   setup_process(&argc,&argv);
      
   /* setup lattice geometry */
   geometry_mpi_eo();
   /* test_geometry_mpi_eo(); */
   
    /* setup random numbers */
    read_input(rlx_var.read,"input_file");
    lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
    rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */

   
   lprintf("MAIN",0,"Action of Qhat on plane waves\n");
   lprintf("MAIN",0,"-----------------------------\n\n");
   
   
   u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f(&glattice);
#endif

	
   unit_u(u_gauge);
   start_gf_sendrecv(u_gauge);

   represent_gauge_field();
   
   ps0=alloc_spinor_field_f(3,&glattice);
   ps1=ps0+1;
   ps2=ps1+1;
   
   pi=4.0*atan(1.0);
   n=10;
   bc_t=bc_x=bc_y=bc_z=0;
#ifdef BC_T_ANTIPERIODIC
   bc_t=1;
#endif
#ifdef BC_X_ANTIPERIODIC
   bc_x=1;
#endif
#ifdef BC_Y_ANTIPERIODIC
   bc_y=1;
#endif
#ifdef BC_Z_ANTIPERIODIC
   bc_z=1;
#endif
   lprintf("MAIN",0,"Anti periodic boundary conditions: %d %d %d %d\n",bc_t,bc_x,bc_y,bc_z);
   
   for (i=0;i<n;i++)
   {
      ranlxd(ran,4);
      
      np[0]=(int)(ran[0]*(double)(GLB_T));
      np[1]=(int)(ran[1]*(double)(GLB_X));
      np[2]=(int)(ran[2]*(double)(GLB_Y));
      np[3]=(int)(ran[3]*(double)(GLB_Z));

      p[0]=((double)(np[0])*2.0+bc_t)*pi/(double)(GLB_T);
      p[1]=((double)(np[1])*2.0+bc_x)*pi/(double)(GLB_X);
      p[2]=((double)(np[2])*2.0+bc_y)*pi/(double)(GLB_Y);
      p[3]=((double)(np[3])*2.0+bc_z)*pi/(double)(GLB_Z);

      mp=(double)(hmass);
      mp+=(double)(1.0-cos(p[0]));
      mp+=(double)(1.0-cos(p[1]));
      mp+=(double)(1.0-cos(p[2]));
      mp+=(double)(1.0-cos(p[3]));      
      
      sp[0]=(double)(sin(p[0]));
      sp[1]=(double)(sin(p[1]));
      sp[2]=(double)(sin(p[2]));
      sp[3]=(double)(sin(p[3]));
      
      rs=(double*)(&s);
      r=0.0f;
      while ((1.0f+r)==1.0f)
	{
	  gauss(rs,8*NF);
	  r=0.0f;
	  
	  for (j=0;j<8*NF;j++)
	    r+=rs[j]*rs[j];
	  
	  r=(double)(sqrt((double)(r)));
	}
      
      
      /*We define directly the spinor also on the buffer so that the test can work with or without mpi*/
      for (x0=-T_BORDER;x0<T+T_BORDER;x0++) 
	for (x1=-X_BORDER;x1<X+X_BORDER;x1++)
	  for (x2=-Y_BORDER;x2<Y+Y_BORDER;x2++)
	    for (x3=-Z_BORDER;x3<Z+Z_BORDER;x3++)
	      {
		ix=ipt(x0,x1,x2,x3);
		if(ix==-1 || ix >= glattice.gsize_spinor) continue;
		
		/* Attention, the definition of the plane wave depends on the slice used for the BC*/
		px=p[0]*(double)(safe_mod(x0+zerocoord[0]-T_BORDER-1,GLB_T)+T_BORDER+1)
		  +p[1]*(double)(safe_mod(x1+zerocoord[1]-X_BORDER-1,GLB_X)+X_BORDER+1)
		  +p[2]*(double)(safe_mod(x2+zerocoord[2]-Y_BORDER-1,GLB_Y)+Y_BORDER+1)
		  +p[3]*(double)(safe_mod(x3+zerocoord[3]-Z_BORDER-1,GLB_Z)+Z_BORDER+1);
		
		z=(cos(px))+I*(sin(px));
		
		_vector_mulc_f(s0.c[0],z,s.c[0]);
		_vector_mulc_f(s0.c[1],z,s.c[1]);
		_vector_mulc_f(s0.c[2],z,s.c[2]);
		_vector_mulc_f(s0.c[3],z,s.c[3]);
		
		*_FIELD_AT(ps0,ix) = s0;
		
		z=mp;
		
		_vector_mulc_f(s1.c[0],z,s0.c[0]);
		_vector_mulc_f(s1.c[1],z,s0.c[1]);
		_vector_mulc_f(s1.c[2],z,s0.c[2]);
		_vector_mulc_f(s1.c[3],z,s0.c[3]);
		
		for (mu=0;mu<4;mu++) {
		  s2=mul_gamma(mu,s0);
		  
		  z=I*sp[mu];
		  
		  _vector_mulc_f(s3.c[0],z,s2.c[0]);
		  _vector_mulc_f(s3.c[1],z,s2.c[1]);
		  _vector_mulc_f(s3.c[2],z,s2.c[2]);
		  _vector_mulc_f(s3.c[3],z,s2.c[3]);
		  
		  _vector_add_assign_f(s1.c[0],s3.c[0]);
		  _vector_add_assign_f(s1.c[1],s3.c[1]);
		  _vector_add_assign_f(s1.c[2],s3.c[2]);
		  _vector_add_assign_f(s1.c[3],s3.c[3]);
		}
		*_FIELD_AT(ps1,ix) = s1;
	      }
#pragma omp parallel default(shared)
{      
      Dphi(hmass,ps2,ps0);
}
     
      spinor_field_mul_add_assign_f(ps1,-1.0,ps2);
      sig=spinor_field_sqnorm_f(ps1)/spinor_field_sqnorm_f(ps0);
      
      lprintf("MAIN",0,"Maximal normalized difference = %.2e at p=(%d,%d,%d,%d)\n",sqrt(sig),
             np[0],np[1],np[2],np[3]);
      lprintf("MAIN",0,"should be around 1*10^(-15) or so)\n\n");
   }
   
   finalize_process();
   
   exit(0);
   
   
}
