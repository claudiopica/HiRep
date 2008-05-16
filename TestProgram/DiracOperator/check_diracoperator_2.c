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

double hmass=0.1;

static suNf_spinor mul_gamma(int mu,suNf_spinor s)
{
   suNf_spinor r;
   complex i,m_i,m_1;

   i.re=0.0f;
   i.im=1.0f;

   m_i.re=0.0f;
   m_i.im=-1.0f;

   m_1.re=-1.0f;
   m_1.im=0.0f;

   if (mu==0)
   {
      _vector_mulc_f(r.c[0],m_1,s.c[2]);
      _vector_mulc_f(r.c[1],m_1,s.c[3]);
      _vector_mulc_f(r.c[2],m_1,s.c[0]);
      _vector_mulc_f(r.c[3],m_1,s.c[1]);
   }
   else if (mu==1)
   {
      _vector_mulc_f(r.c[0],m_i,s.c[3]);
      _vector_mulc_f(r.c[1],m_i,s.c[2]);
      _vector_mulc_f(r.c[2],i,s.c[1]);
      _vector_mulc_f(r.c[3],i,s.c[0]);
   }
   else if (mu==2)
   {
      _vector_mulc_f(r.c[0],m_1,s.c[3]);
      r.c[1]=s.c[2];
      r.c[2]=s.c[1];
      _vector_mulc_f(r.c[3],m_1,s.c[0]);
   }
   else if (mu==3)
   {
      _vector_mulc_f(r.c[0],m_i,s.c[2]);
      _vector_mulc_f(r.c[1],i,s.c[3]);
      _vector_mulc_f(r.c[2],i,s.c[0]);
      _vector_mulc_f(r.c[3],m_i,s.c[1]);
   }
   else
   {
      r.c[0]=s.c[0];
      r.c[1]=s.c[1];
      _vector_mulc_f(r.c[2],m_1,s.c[2]);
      _vector_mulc_f(r.c[3],m_1,s.c[3]);
   }

   return r;
}


int main(int argc,char *argv[])
{
   int i,j,n,ix,mu,level,seed;
   int bc_t,bc_x,bc_y,bc_z;
   int x0,x1,x2,x3;
   int np[4];
   double ran[4],sp[4];
   double pi,p[4];
   double *rs,r,mp,sig,px;
   complex z;
   suNf_spinor s,s0,s1,s2,s3;
   spinor_field *ps0,*ps1,*ps2;

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   geometry_mpi_eo();
   printf("The coordinates of the node are (%d,%d,%d,%d)\n",proc_t(myid),proc_x(myid),proc_y(myid),proc_z(myid));
   printf("The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
   printf("The lattice global size is %dx%dx%dx%d\n",GLOBAL_T,GLOBAL_X,GLOBAL_Y,GLOBAL_Z);
   printf("The lattice borders are (%d,%d,%d,%d)\n",T_BORDER,X_BORDER,Y_BORDER,Z_BORDER);
   printf("\n");
   printf("\n");
   printf("Action of Qhat on plane waves\n");
   printf("-----------------------------\n\n");

   level=1;
   seed=123;
   rlxd_init(level,seed);
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
 
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   printf("allocating gfield_f\n");
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();
	 
	 ps0=alloc_spinor_field_f(3,&glattice);
	 ps1=ps0+1;
	 ps2=ps1+1;

   pi=4.0*atan(1.0);
   n=10;
   bc_t=bc_x=bc_y=bc_z=0;
#ifdef ANTIPERIODIC_BC_T
   bc_t=1;
#endif
#ifdef ANTIPERIODIC_BC_X
   bc_x=1;
#endif
#ifdef ANTIPERIODIC_BC_Y
   bc_y=1;
#endif
#ifdef ANTIPERIODIC_BC_Z
   bc_z=1;
#endif
   printf("bc=%d %d %d %d\n",bc_t,bc_x,bc_y,bc_z);
   
   for (i=0;i<n;i++)
   {
      ranlxd(ran,4);
      
      np[0]=(int)(ran[0]*(double)(GLOBAL_T));
      np[1]=(int)(ran[1]*(double)(GLOBAL_X));
      np[2]=(int)(ran[2]*(double)(GLOBAL_Y));
      np[3]=(int)(ran[3]*(double)(GLOBAL_Z));

      p[0]=((double)(np[0])*2.0+bc_t)*pi/(double)(GLOBAL_T);
      p[1]=((double)(np[1])*2.0+bc_x)*pi/(double)(GLOBAL_X);
      p[2]=((double)(np[2])*2.0+bc_y)*pi/(double)(GLOBAL_Y);
      p[3]=((double)(np[3])*2.0+bc_z)*pi/(double)(GLOBAL_Z);

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
      
      for (x0=-T_BORDER;x0<T+T_BORDER;x0++) /* DEFINIAMO LO SPINORE ANCHE SUI BUFFERS, IN MODO CHE IL TEST FUNZIONI CON O SENZA MPI */
      for (x1=-X_BORDER;x1<X+X_BORDER;x1++)
      for (x2=-Y_BORDER;x2<Y+Y_BORDER;x2++)
      for (x3=-Z_BORDER;x3<Z+Z_BORDER;x3++)
      {
         ix=ipt(x0,x1,x2,x3);
         if(ix==-1) continue;
   
/* LA DEFINIZIONE DELL'ONDA PIANA DIPENDE DA QUALI TIME SLICES HO USATO PER LE BC */    
      
         px=p[0]*(double)(SAFE_MOD(global_t(myid,x0)-T_BORDER-1,GLOBAL_T)+T_BORDER+1)
            +p[1]*(double)(SAFE_MOD(global_x(myid,x1)-X_BORDER-1,GLOBAL_X)+X_BORDER+1)
            +p[2]*(double)(SAFE_MOD(global_y(myid,x2)-Y_BORDER-1,GLOBAL_Y)+Y_BORDER+1)
            +p[3]*(double)(SAFE_MOD(global_z(myid,x3)-Z_BORDER-1,GLOBAL_Z)+Z_BORDER+1);

         z.re=(double)(cos(px));
         z.im=(double)(sin(px));
         
         _vector_mulc_f(s0.c[0],z,s.c[0]);
         _vector_mulc_f(s0.c[1],z,s.c[1]);
         _vector_mulc_f(s0.c[2],z,s.c[2]);
         _vector_mulc_f(s0.c[3],z,s.c[3]);
         
         *_SPINOR_AT(ps0,ix) = s0;
         
         z.re=mp;
         z.im=0.0f;
         
         _vector_mulc_f(s1.c[0],z,s0.c[0]);
         _vector_mulc_f(s1.c[1],z,s0.c[1]);
         _vector_mulc_f(s1.c[2],z,s0.c[2]);
         _vector_mulc_f(s1.c[3],z,s0.c[3]);
         
         for (mu=0;mu<4;mu++) {
            s2=mul_gamma(mu,s0);
            
            z.re=0.0f;
            z.im=sp[mu];

            _vector_mulc_f(s3.c[0],z,s2.c[0]);
            _vector_mulc_f(s3.c[1],z,s2.c[1]);
            _vector_mulc_f(s3.c[2],z,s2.c[2]);
            _vector_mulc_f(s3.c[3],z,s2.c[3]);
            
            _vector_add_assign_f(s1.c[0],s3.c[0]);
            _vector_add_assign_f(s1.c[1],s3.c[1]);
            _vector_add_assign_f(s1.c[2],s3.c[2]);
            _vector_add_assign_f(s1.c[3],s3.c[3]);
         }
         *_SPINOR_AT(ps1,ix) = s1;
      }
      
      Dphi(hmass,ps2,ps0);
      
      spinor_field_mul_add_assign_f(ps1,-1.0,ps2);
      sig=spinor_field_sqnorm_f(ps1)/spinor_field_sqnorm_f(ps0);
         
      printf("Maximal normalized difference = %.2e at p=(%d,%d,%d,%d)\n",sqrt(sig),
             np[0],np[1],np[2],np[3]);
      printf("(should be around 1*10^(-15) or so)\n\n");
   }
   
   free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif
	 free_field(ps0);

   exit(0);
}
