/*! \file
 * \brief Routines for the SF gauge observable
 *
 * Routines for the Schrodinger Functional gauge observable.
 *
 */

/*! \defgroup sfgauge SF Gauge Observable
 * \ingroup obs
 */

#include "global.h"
#include "suN.h"
#include "error.h"
#include "communications.h"
#include "geometry.h"
#include "observables.h"
#include <stdio.h>

void print_matrix(suNg v, char label[])
{
   int i;
   printf("\n\n SuNg Matrix:\n");
   printf(label);
   for(i=0;i<NG*NG;i++)
   {
      if((i%NG)==0)
      {
         printf("\n");
      }
      printf("[%+.5f,%+.5f] ",((v).c[i]).re,((v).c[i]).im);
   }
}

/*! \ingroup sfgauge
 * Bottom E_8 component
 *
 * \f$ E_k^8({\bf x})= {\frac{1}{a^2}}\mathrm{Re}\, \mathrm{tr}\, \left\{ i \lambda_8 U(x,k) U(x+a\hat k,0) U(x+a\hat 0,k)^{-1} U(x,0)^{-1} \right\}_{x_0=0}\f$
 *
 * The matrix \f$\lambda_8\f$ depends on the size of the gauge group, for SU(3) it is given by \f$\lambda_8={\rm diag}(1,-1/2,-1/2)\f$
 *
 */

double E_8(int ix,int k)
{
   double p;
   suNg *v1,*v2,*v3,*v4,w1,w2,w3,w4,ilambda8;
   int iy, iz;
/*   printf("\n\n E_8: ix=%d, k=%d \n\n",ix,k);*/

if(NG==4)
{
   _suNg_zero(ilambda8);
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 0.5;
   ((ilambda8).c[5]).re = 0.0;
   ((ilambda8).c[5]).im = 0.5;
   ((ilambda8).c[10]).re = 0.0;
   ((ilambda8).c[10]).im = -0.5;
   ((ilambda8).c[15]).re = 0.0;
   ((ilambda8).c[15]).im = -0.5;
}
else if(NG==3)
{
   _suNg_zero(ilambda8);
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 1.0;
   ((ilambda8).c[4]).re = 0.0;
   ((ilambda8).c[4]).im = -0.5;
   ((ilambda8).c[8]).re = 0.0;
   ((ilambda8).c[8]).im = -0.5;
}
else if (NG==2)
{
   _suNg_zero(ilambda8);
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 1.0;
   ((ilambda8).c[3]).re = 0.0;
   ((ilambda8).c[3]).im = -1.0;
}
else
{
   error(1,1,"E [SF_action.c]",
         "No explicit form for observable for given NG");
}  

   iy=iup(ix,k);
   iz=iup(ix,0);

   v1=pu_gauge(ix,k);
/*   print_matrix((*v1),"v1\n");*/
   v2=pu_gauge(iy,0);
/*   print_matrix((*v2),"v2\n");*/
   v3=pu_gauge(iz,k);
/*   print_matrix((*v3),"v3\n");*/
   v4=pu_gauge(ix,0);
/*   print_matrix((*v4),"v4\n");*/

   _suNg_times_suNg(w1,(*v1),(*v2));
/*   print_matrix(w1,"w1=v1v2\n");*/
   _suNg_times_suNg(w2,(*v4),(*v3));
/*   print_matrix(w2,"w2=v4v3\n");*/

   _suNg_times_suNg_dagger(w3,w1,w2);      
/*   print_matrix(w3,"w3=w1(w2)'=v1v2(v3)'(v4)'\n");*/
      
   _suNg_times_suNg(w4,ilambda8,w3);
/*   print_matrix((*v4),"i lambda_8 w3\n");*/

   _suNg_trace_re(p,w4);
/*   printf("Re Tr = %f \n",p);*/

   return p;
}

/*! \ingroup sfgauge
 * Top E_8 component
 *
 * \f$ (E_k^8)'({\bf x})= -{\frac{1}{a^2}}\mathrm{Re}\, \mathrm{tr}\, \left\{ i \lambda_8 U(x+a\hat 0,k)^{-1} U(x,0)^{-1} U(x,k) U(x+a\hat k,0) \right\}_{x_0=T-2a}\f$
 *
 * The matrix \f$\lambda_8\f$ depends on the size of the gauge group, for SU(3) it is given by \f$\lambda_8={\rm diag}(1,-1/2,-1/2)\f$

 */

double E_8_top(int ix,int k)
{
   double p;
   suNg *v1,*v2,*v3,*v4,w1,w2,w3,w4,ilambda8;
   int iy, iz;
/*   printf("\n\n E_8: ix=%d, k=%d \n\n",ix,k);*/

if(NG==4)
{
   _suNg_zero(ilambda8);
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 0.5;
   ((ilambda8).c[5]).re = 0.0;
   ((ilambda8).c[5]).im = 0.5;
   ((ilambda8).c[10]).re = 0.0;
   ((ilambda8).c[10]).im = -0.5;
   ((ilambda8).c[15]).re = 0.0;
   ((ilambda8).c[15]).im = -0.5;
}
else if(NG==3)
{
   _suNg_zero(ilambda8);
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 1.0;
   ((ilambda8).c[4]).re = 0.0;
   ((ilambda8).c[4]).im = -0.5;
   ((ilambda8).c[8]).re = 0.0;
   ((ilambda8).c[8]).im = -0.5;
}
else if (NG==2)
{
   _suNg_zero(ilambda8);
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 1.0;
   ((ilambda8).c[3]).re = 0.0;
   ((ilambda8).c[3]).im = -1.0;
}
else
{
   error(1,1,"E [SF_action.c]",
         "No explicit form for observable for given NG");
}  
   iy=iup(ix,k);
   iz=iup(ix,0);

   v1=pu_gauge(ix,k);
/*   print_matrix((*v1),"v1\n");*/
   v2=pu_gauge(iy,0);
/*   print_matrix((*v2),"v2\n");*/
   v3=pu_gauge(iz,k);
/*   print_matrix((*v3),"v3\n");*/
   v4=pu_gauge(ix,0);
/*   print_matrix((*v4),"v4\n");*/

   _suNg_times_suNg(w1,(*v1),(*v2));
/*   print_matrix(w1,"w1=v1v2\n");*/
   _suNg_times_suNg(w2,ilambda8,(*v3));
/*   print_matrix(w2,"w2=v4v3\n");*/

   _suNg_times_suNg_dagger(w3,w2,w1);      
/*   print_matrix(w3,"w3=w1(w2)'=v1v2(v3)'(v4)'\n");*/
      
   _suNg_times_suNg(w4,(*v4),w3);
/*   print_matrix((*v4),"i lambda_8 w3\n");*/

   _suNg_trace_re(p,w4);
/*   printf("Re Tr = %f \n",p);*/
   return -p;
}

/*! \ingroup sfgauge
 *  SF Gauge Observable
 *
 * \f$\frac{\partial S_g}{\partial\eta}= - {2 \over g_0^2 \, L} a^3 \sum_{\bf x}\left\{E_k^8({\bf x}) - (E_k^8)'({\bf x}) \right\}\f$
 *
 * where
 *
 * \f$ E_k^8({\bf x})= {\frac{1}{a^2}}\mathrm{Re}\, \mathrm{tr}\, \left\{ i \lambda_8 U(x,k) U(x+a\hat k,0) U(x+a\hat 0,k)^{-1} U(x,0)^{-1} \right\}_{x_0=0}\f$
 *
 * \f$ (E_k^8)'({\bf x})= -{\frac{1}{a^2}}\mathrm{Re}\, \mathrm{tr}\, \left\{ i \lambda_8 U(x+a\hat 0,k)^{-1} U(x,0)^{-1} U(x,k) U(x+a\hat k,0) \right\}_{x_0=T-2a}\f$
 *
 * The matrix \f$\lambda_8\f$ depends on the size of the gauge group, for SU(3) it is given by \f$\lambda_8={\rm diag}(1,-1/2,-1/2)\f$

 */

double sf_action(double beta)
{
/*single processor method*/
/*
  int x, y, z, k;
  double pa=0.0;
  
	for (x=0;x<GLB_X;x++)
	for (y=0;y<GLB_Y;y++)
	for (z=0;z<GLB_Z;z++)
		for (k=1;k<4;k++)
			pa+=(double)(E_8(ipt(0,x,y,z),k) + E_8_top(ipt(GLB_T-2,x,y,z),k));

  return pa*(double)(beta/(NG*GLB_X));
*/
/*mpi method*/
  _DECLARE_INT_ITERATOR(i);
  double pa=0.;
  int ix, iy, iz, k;

if(COORD[0]==0)
{
  _PIECE_FOR(&glattice,i) {
    _SITE_FOR(&glattice,i) {
    
	for (ix=0; ix<GLB_X/NP_X; ++ix)
        for (iy=0; iy<GLB_Y/NP_Y; ++iy)
        for (iz=0; iz<GLB_Z/NP_Z; ++iz)
	{
	{
	{
		if (ipt(1,ix,iy,iz)==i)/*TEMP ix,iy,iz->0,0,0 etc*/
		{
			for (k=1;k<4;k++)
			{
				pa+=(double)(E_8(i,k));
			}
    		}
	}
	}
	}

    }
    if(_PIECE_INDEX(i)==0) {
      complete_gf_sendrecv(u_gauge);
    }
  }
}
if(COORD[0]==NP_T-1)
{
  _PIECE_FOR(&glattice,i) {
    _SITE_FOR(&glattice,i) {
    
	for (ix=0; ix<GLB_X/NP_X; ++ix)
        for (iy=0; iy<GLB_Y/NP_Y; ++iy)
        for (iz=0; iz<GLB_Z/NP_Z; ++iz)
	{
	{
	{
		if (ipt((GLB_T/NP_T)-2,ix,iy,iz)==i)
		{
			for (k=1;k<4;k++)
			{
				pa+=(double)(E_8_top(i,k));
			}
    		}
	}
	}
	}

    }
    if(_PIECE_INDEX(i)==0) {
      complete_gf_sendrecv(u_gauge);
    }
  }
}
  global_sum(&pa, 1);

  return pa*(double)(beta/(NG*GLB_X));

}

double sf_test_gauge_bcs()
{
/*calculates average of all plaquettes that should remain fixed for SF*/
  _DECLARE_INT_ITERATOR(i);
  double pa=0.;
  int ix, iy, iz;

if(COORD[0]==0)
{
  _PIECE_FOR(&glattice,i) {
    _SITE_FOR(&glattice,i) {
    
	for (ix=0; ix<GLB_X/NP_X; ++ix)
        for (iy=0; iy<GLB_Y/NP_Y; ++iy)
        for (iz=0; iz<GLB_Z/NP_Z; ++iz)
	{
	{
	{
		if (ipt(1,ix,iy,iz)==i)
		{
			pa+=(double)(plaq(i,2,1));
			pa+=(double)(plaq(i,3,1));
			pa+=(double)(plaq(i,3,2));
    		}
/*		if (ipt(0,ix,iy,iz)==i)
		{
			pa+=(double)(plaq(i,1,0));
			pa+=(double)(plaq(i,2,0));
			pa+=(double)(plaq(i,2,1));
			pa+=(double)(plaq(i,3,0));
			pa+=(double)(plaq(i,3,1));
			pa+=(double)(plaq(i,3,2));
    		}
*/		
	}
	}
	}

    }
    if(_PIECE_INDEX(i)==0) {
      complete_gf_sendrecv(u_gauge);
    }
  }
}
if(COORD[0]==NP_T-1)
{
  _PIECE_FOR(&glattice,i) {
    _SITE_FOR(&glattice,i) {
    
	for (ix=0; ix<GLB_X/NP_X; ++ix)
        for (iy=0; iy<GLB_Y/NP_Y; ++iy)
        for (iz=0; iz<GLB_Z/NP_Z; ++iz)
	{
	{
	{
		if (ipt((GLB_T/NP_T)-1,ix,iy,iz)==i)
		{
/*			pa+=(double)(plaq(i,1,0));
			pa+=(double)(plaq(i,2,0));
			pa+=(double)(plaq(i,2,1));
*/			pa+=(double)(plaq(i,3,0));
			pa+=(double)(plaq(i,3,1));
			pa+=(double)(plaq(i,3,2));
    		}
	}
	}
	}

    }
    if(_PIECE_INDEX(i)==0) {
      complete_gf_sendrecv(u_gauge);
    }
  }
}
  global_sum(&pa, 1);

  return pa/(double)(GLB_X*GLB_Y*GLB_Z*NG*(6+6+3));

}
