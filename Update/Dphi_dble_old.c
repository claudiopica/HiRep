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
void Dphi_dble_old(double m0, suNf_spinor_dble *p, suNf_spinor_dble *q)
{
   int ix,iy;
   double rho;
   suNf_dble *up,*um;
   suNf_vector_dble psi,chi;
   suNf_spinor_dble *r, *s,*sp,*sm;

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

      _vector_mul_f((*r).c1,rho,(*s).c1);
      _vector_mul_f((*r).c2,rho,(*s).c2);
      _vector_mul_f((*r).c3,rho,(*s).c3);
      _vector_mul_f((*r).c4,rho,(*s).c4);

/******************************* direction +0 *********************************/

      iy=iup[ix][0];
      sp=q+iy;
      up=pu_gauge_f_dble(ix,0);
      
      _vector_add_f(psi,(*sp).c1,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);
      
      _vector_add_assign_f((*r).c1,chi);
      _vector_add_assign_f((*r).c3,chi);

      _vector_add_f(psi,(*sp).c2,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);
            
      _vector_add_assign_f((*r).c2,chi);
      _vector_add_assign_f((*r).c4,chi);

/******************************* direction -0 *********************************/

      iy=idn[ix][0];
      sm=q+iy;
      um=pu_gauge_f_dble(iy,0);
      
      _vector_sub_f(psi,(*sm).c1,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_sub_assign_f((*r).c3,chi);

      _vector_sub_f(psi,(*sm).c2,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c2,chi);
      _vector_sub_assign_f((*r).c4,chi);

/******************************* direction +1 *********************************/

      iy=iup[ix][1];
      sp=q+iy;
      up=pu_gauge_f_dble(ix,1);      
      
      _vector_i_add_f(psi,(*sp).c1,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_i_sub_assign_f((*r).c4,chi);

      _vector_i_add_f(psi,(*sp).c2,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_sub_assign_f((*r).c3,chi);

/******************************* direction -1 *********************************/

      iy=idn[ix][1];
      sm=q+iy;
      um=pu_gauge_f_dble(iy,1);
      
      _vector_i_sub_f(psi,(*sm).c1,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_i_add_assign_f((*r).c4,chi);

      _vector_i_sub_f(psi,(*sm).c2,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_add_assign_f((*r).c3,chi);

/******************************* direction +2 *********************************/

      iy=iup[ix][2];
      sp=q+iy;
      up=pu_gauge_f_dble(ix,2);
      
      _vector_add_f(psi,(*sp).c1,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_add_assign_f((*r).c4,chi);

      _vector_sub_f(psi,(*sp).c2,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);
      
      _vector_add_assign_f((*r).c2,chi);
      _vector_sub_assign_f((*r).c3,chi);

/******************************* direction -2 *********************************/

      iy=idn[ix][2];
      sm=q+iy;
      um=pu_gauge_f_dble(iy,2);
      
      _vector_sub_f(psi,(*sm).c1,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_sub_assign_f((*r).c4,chi);

      _vector_add_f(psi,(*sm).c2,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c2,chi);
      _vector_add_assign_f((*r).c3,chi);

/******************************* direction +3 *********************************/

      iy=iup[ix][3];
      sp=q+iy;
      up=pu_gauge_f_dble(ix,3);
      
      _vector_i_add_f(psi,(*sp).c1,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_i_sub_assign_f((*r).c3,chi);

      _vector_i_sub_f(psi,(*sp).c2,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_add_assign_f((*r).c4,chi);

/******************************* direction -3 *********************************/

      iy=idn[ix][3];
      sm=q+iy;
      um=pu_gauge_f_dble(iy,3);
      
      _vector_i_sub_f(psi,(*sm).c1,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c1,chi);
      _vector_i_add_assign_f((*r).c3,chi);

      _vector_i_add_f(psi,(*sm).c2,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_sub_assign_f((*r).c4,chi);
      
/******************************** end of loop *********************************/

      _vector_mul_f((*r).c1,-0.5,(*r).c1);
      _vector_mul_f((*r).c2,-0.5,(*r).c2);
      _vector_mul_f((*r).c3,-0.5,(*r).c3);
      _vector_mul_f((*r).c4,-0.5,(*r).c4);
   }
}

void g5Dphi_dble_old(double m0, suNf_spinor_dble *p, suNf_spinor_dble *q)
{
   int ix,iy;
   double rho;
   suNf_dble *up,*um;
   suNf_vector_dble psi,chi;
   suNf_spinor_dble *r, *s,*sp,*sm;

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

      _vector_mul_f((*r).c1,rho,(*s).c1);
      _vector_mul_f((*r).c2,rho,(*s).c2);
      _vector_mul_f((*r).c3,rho,(*s).c3);
      _vector_mul_f((*r).c4,rho,(*s).c4);

/******************************* direction +0 *********************************/

      iy=iup[ix][0];
      sp=q+iy;
      up=pu_gauge_f_dble(ix,0);
      
      _vector_add_f(psi,(*sp).c1,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);
      
      _vector_add_assign_f((*r).c1,chi);
      _vector_add_assign_f((*r).c3,chi);

      _vector_add_f(psi,(*sp).c2,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);
            
      _vector_add_assign_f((*r).c2,chi);
      _vector_add_assign_f((*r).c4,chi);

/******************************* direction -0 *********************************/

      iy=idn[ix][0];
      sm=q+iy;
      um=pu_gauge_f_dble(iy,0);
      
      _vector_sub_f(psi,(*sm).c1,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_sub_assign_f((*r).c3,chi);

      _vector_sub_f(psi,(*sm).c2,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c2,chi);
      _vector_sub_assign_f((*r).c4,chi);

/******************************* direction +1 *********************************/

      iy=iup[ix][1];
      sp=q+iy;
      up=pu_gauge_f_dble(ix,1);      
      
      _vector_i_add_f(psi,(*sp).c1,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_i_sub_assign_f((*r).c4,chi);

      _vector_i_add_f(psi,(*sp).c2,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_sub_assign_f((*r).c3,chi);

/******************************* direction -1 *********************************/

      iy=idn[ix][1];
      sm=q+iy;
      um=pu_gauge_f_dble(iy,1);
      
      _vector_i_sub_f(psi,(*sm).c1,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_i_add_assign_f((*r).c4,chi);

      _vector_i_sub_f(psi,(*sm).c2,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_add_assign_f((*r).c3,chi);

/******************************* direction +2 *********************************/

      iy=iup[ix][2];
      sp=q+iy;
      up=pu_gauge_f_dble(ix,2);
      
      _vector_add_f(psi,(*sp).c1,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_add_assign_f((*r).c4,chi);

      _vector_sub_f(psi,(*sp).c2,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);
      
      _vector_add_assign_f((*r).c2,chi);
      _vector_sub_assign_f((*r).c3,chi);

/******************************* direction -2 *********************************/

      iy=idn[ix][2];
      sm=q+iy;
      um=pu_gauge_f_dble(iy,2);
      
      _vector_sub_f(psi,(*sm).c1,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_sub_assign_f((*r).c4,chi);

      _vector_add_f(psi,(*sm).c2,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c2,chi);
      _vector_add_assign_f((*r).c3,chi);

/******************************* direction +3 *********************************/

      iy=iup[ix][3];
      sp=q+iy;
      up=pu_gauge_f_dble(ix,3);
      
      _vector_i_add_f(psi,(*sp).c1,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c1,chi);
      _vector_i_sub_assign_f((*r).c3,chi);

      _vector_i_sub_f(psi,(*sp).c2,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_add_assign_f((*r).c4,chi);

/******************************* direction -3 *********************************/

      iy=idn[ix][3];
      sm=q+iy;
      um=pu_gauge_f_dble(iy,3);
      
      _vector_i_sub_f(psi,(*sm).c1,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c1,chi);
      _vector_i_add_assign_f((*r).c3,chi);

      _vector_i_add_f(psi,(*sm).c2,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_sub_assign_f((*r).c4,chi);
      
/******************************** end of loop *********************************/

      _vector_mul_f((*r).c1,-0.5,(*r).c1);
      _vector_mul_f((*r).c2,-0.5,(*r).c2);
      _vector_mul_f((*r).c3,0.5,(*r).c3);
      _vector_mul_f((*r).c4,0.5,(*r).c4);
   }
}

