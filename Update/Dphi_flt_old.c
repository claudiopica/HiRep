/*******************************************************************************
*
* File Dphi_dble.c
*
* Action of the Wilson-Dirac operator D and hermitian g5D on a given 
* double-precision spinor field
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "error.h"
#include "global.h"
#include "dirac.h"

/* p = out ; q = in */
void Dphi_flt_old(double m0, suNf_spinor_flt *p, suNf_spinor_flt *q)
{
   int ix,iy;
   double rho;
   suNf_flt *up,*um;
   suNf_vector_flt psi,chi;
   suNf_spinor_flt *r, *s,*sp,*sm;

   error((q==NULL)||(p==NULL),1,"Qphi [Qphi.c]",
         "Attempt to access unallocated memory space");
   
   error(q==p,1,"Qphi [Qphi.c]",
         "Input and output fields must be different");
   
   rho=-8.0-2.0*m0;
   r=p-1;
   s=q-1;

/************************ loop over all lattice sites *************************/

   for (ix=0;ix<VOLUME;ix++) 
   {
      ++r;
      ++s;

      _vector_mul_f((*r).c[0],rho,(*s).c[0]);
      _vector_mul_f((*r).c[1],rho,(*s).c[1]);
      _vector_mul_f((*r).c[2],rho,(*s).c[2]);
      _vector_mul_f((*r).c[3],rho,(*s).c[3]);

/******************************* direction +0 *********************************/

      iy=iup(ix,0);
      sp=q+iy;
      up=pu_gauge_f_flt(ix,0);
      
      _vector_add_f(psi,(*sp).c[0],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);
      
      _vector_add_assign_f((*r).c[0],chi);
      _vector_add_assign_f((*r).c[2],chi);

      _vector_add_f(psi,(*sp).c[1],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);
            
      _vector_add_assign_f((*r).c[1],chi);
      _vector_add_assign_f((*r).c[3],chi);

/******************************* direction -0 *********************************/

      iy=idn(ix,0);
      sm=q+iy;
      um=pu_gauge_f_flt(iy,0);
      
      _vector_sub_f(psi,(*sm).c[0],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_sub_assign_f((*r).c[2],chi);

      _vector_sub_f(psi,(*sm).c[1],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c[1],chi);
      _vector_sub_assign_f((*r).c[3],chi);

/******************************* direction +1 *********************************/

      iy=iup(ix,1);
      sp=q+iy;
      up=pu_gauge_f_flt(ix,1);      
      
      _vector_i_add_f(psi,(*sp).c[0],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_sub_assign_f((*r).c[3],chi);

      _vector_i_add_f(psi,(*sp).c[1],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_sub_assign_f((*r).c[2],chi);

/******************************* direction -1 *********************************/

      iy=idn(ix,1);
      sm=q+iy;
      um=pu_gauge_f_flt(iy,1);
      
      _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_add_assign_f((*r).c[3],chi);

      _vector_i_sub_f(psi,(*sm).c[1],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_add_assign_f((*r).c[2],chi);

/******************************* direction +2 *********************************/

      iy=iup(ix,2);
      sp=q+iy;
      up=pu_gauge_f_flt(ix,2);
      
      _vector_add_f(psi,(*sp).c[0],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_add_assign_f((*r).c[3],chi);

      _vector_sub_f(psi,(*sp).c[1],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);
      
      _vector_add_assign_f((*r).c[1],chi);
      _vector_sub_assign_f((*r).c[2],chi);

/******************************* direction -2 *********************************/

      iy=idn(ix,2);
      sm=q+iy;
      um=pu_gauge_f_flt(iy,2);
      
      _vector_sub_f(psi,(*sm).c[0],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_sub_assign_f((*r).c[3],chi);

      _vector_add_f(psi,(*sm).c[1],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c[1],chi);
      _vector_add_assign_f((*r).c[2],chi);

/******************************* direction +3 *********************************/

      iy=iup(ix,3);
      sp=q+iy;
      up=pu_gauge_f_flt(ix,3);
      
      _vector_i_add_f(psi,(*sp).c[0],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_sub_assign_f((*r).c[2],chi);

      _vector_i_sub_f(psi,(*sp).c[1],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_add_assign_f((*r).c[3],chi);

/******************************* direction -3 *********************************/

      iy=idn(ix,3);
      sm=q+iy;
      um=pu_gauge_f_flt(iy,3);
      
      _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_add_assign_f((*r).c[2],chi);

      _vector_i_add_f(psi,(*sm).c[1],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_sub_assign_f((*r).c[3],chi);
      
/******************************** end of loop *********************************/

      _vector_mul_f((*r).c[0],-0.5,(*r).c[0]);
      _vector_mul_f((*r).c[1],-0.5,(*r).c[1]);
      _vector_mul_f((*r).c[2],-0.5,(*r).c[2]);
      _vector_mul_f((*r).c[3],-0.5,(*r).c[3]);
   }
}

void g5Dphi_flt_old(double m0, suNf_spinor_flt *p, suNf_spinor_flt *q)
{
   int ix,iy;
   double rho;
   suNf_flt *up,*um;
   suNf_vector_flt psi,chi;
   suNf_spinor_flt *r, *s,*sp,*sm;

   error((q==NULL)||(p==NULL),1,"Qphi [Qphi.c]",
         "Attempt to access unallocated memory space");
   
   error(q==p,1,"Qphi [Qphi.c]",
         "Input and output fields must be different");
   
   rho=-8.0-2.0*m0;
   r=p-1;
   s=q-1;

/************************ loop over all lattice sites *************************/

   for (ix=0;ix<VOLUME;ix++) 
   {
      ++r;
      ++s;

      _vector_mul_f((*r).c[0],rho,(*s).c[0]);
      _vector_mul_f((*r).c[1],rho,(*s).c[1]);
      _vector_mul_f((*r).c[2],rho,(*s).c[2]);
      _vector_mul_f((*r).c[3],rho,(*s).c[3]);

/******************************* direction +0 *********************************/

      iy=iup(ix,0);
      sp=q+iy;
      up=pu_gauge_f_flt(ix,0);
      
      _vector_add_f(psi,(*sp).c[0],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);
      
      _vector_add_assign_f((*r).c[0],chi);
      _vector_add_assign_f((*r).c[2],chi);

      _vector_add_f(psi,(*sp).c[1],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);
            
      _vector_add_assign_f((*r).c[1],chi);
      _vector_add_assign_f((*r).c[3],chi);

/******************************* direction -0 *********************************/

      iy=idn(ix,0);
      sm=q+iy;
      um=pu_gauge_f_flt(iy,0);
      
      _vector_sub_f(psi,(*sm).c[0],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_sub_assign_f((*r).c[2],chi);

      _vector_sub_f(psi,(*sm).c[1],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c[1],chi);
      _vector_sub_assign_f((*r).c[3],chi);

/******************************* direction +1 *********************************/

      iy=iup(ix,1);
      sp=q+iy;
      up=pu_gauge_f_flt(ix,1);      
      
      _vector_i_add_f(psi,(*sp).c[0],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_sub_assign_f((*r).c[3],chi);

      _vector_i_add_f(psi,(*sp).c[1],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_sub_assign_f((*r).c[2],chi);

/******************************* direction -1 *********************************/

      iy=idn(ix,1);
      sm=q+iy;
      um=pu_gauge_f_flt(iy,1);
      
      _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_add_assign_f((*r).c[3],chi);

      _vector_i_sub_f(psi,(*sm).c[1],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_add_assign_f((*r).c[2],chi);

/******************************* direction +2 *********************************/

      iy=iup(ix,2);
      sp=q+iy;
      up=pu_gauge_f_flt(ix,2);
      
      _vector_add_f(psi,(*sp).c[0],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_add_assign_f((*r).c[3],chi);

      _vector_sub_f(psi,(*sp).c[1],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);
      
      _vector_add_assign_f((*r).c[1],chi);
      _vector_sub_assign_f((*r).c[2],chi);

/******************************* direction -2 *********************************/

      iy=idn(ix,2);
      sm=q+iy;
      um=pu_gauge_f_flt(iy,2);
      
      _vector_sub_f(psi,(*sm).c[0],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_sub_assign_f((*r).c[3],chi);

      _vector_add_f(psi,(*sm).c[1],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c[1],chi);
      _vector_add_assign_f((*r).c[2],chi);

/******************************* direction +3 *********************************/

      iy=iup(ix,3);
      sp=q+iy;
      up=pu_gauge_f_flt(ix,3);
      
      _vector_i_add_f(psi,(*sp).c[0],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_sub_assign_f((*r).c[2],chi);

      _vector_i_sub_f(psi,(*sp).c[1],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_add_assign_f((*r).c[3],chi);

/******************************* direction -3 *********************************/

      iy=idn(ix,3);
      sm=q+iy;
      um=pu_gauge_f_flt(iy,3);
      
      _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_add_assign_f((*r).c[2],chi);

      _vector_i_add_f(psi,(*sm).c[1],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_sub_assign_f((*r).c[3],chi);
      
/******************************** end of loop *********************************/

      _vector_mul_f((*r).c[0],-0.5,(*r).c[0]);
      _vector_mul_f((*r).c[1],-0.5,(*r).c[1]);
      _vector_mul_f((*r).c[2],0.5,(*r).c[2]);
      _vector_mul_f((*r).c[3],0.5,(*r).c[3]);
   }
}

