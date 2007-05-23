/*******************************************************************************
*
* File Dphi.c
*
* Action of the Wilson-Dirac operator D and hermitian g5D on a given 
* single-precision spinor field
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "global.h"
#include "error.h"
#include "dirac.h"

/*
 * NOTE :
 * here we are making the assumption that the geometry is such that
 * all even sites are in the range [0,VOLUME/2[ and all odd sites are
 * in the range [VOLUME/2,VOLUME[
 */

/* prende 2 spinor lunghi VOLUME/2 definiti solo su siti con la stessa parita' */
void Dphi_(block_selector B, suNf_spinor *out, suNf_spinor *in)
{
   int ix,iy, smin=0, smax=0;
   suNf *up,*um;
   suNf_vector psi,chi;
   suNf_spinor *r=0,*sp,*sm;

   error((in==NULL)||(out==NULL),1,"Dphi_ [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi_ [Dphi.c]",
         "Input and output fields must be different");

   switch(B) {
      case EO:
         in-=VOLUME/2;
         smin=0;
         smax=VOLUME/2;
         break;
      case OE:
         /* in=in; */
         smin=VOLUME/2;
         smax=VOLUME;
         break;
      default:
         error(1,1,"Dphi_ [Dphi.c]",
               "Invalid block parity selection");
   }

   r=out;
  

/************************ loop over all lattice sites *************************/

   for (ix=smin;ix<smax;++ix) 
   {

/******************************* direction +0 *********************************/

      iy=iup[ix][0];
      sp=in+iy;
      up=pu_gauge_f(ix,0);
      
      _vector_add_f(psi,(*sp).c1,(*sp).c3);
      _suNf_multiply(chi,(*up),psi);
      
      (*r).c1=chi;
      (*r).c3=chi;

      _vector_add_f(psi,(*sp).c2,(*sp).c4);
      _suNf_multiply(chi,(*up),psi);
            
      (*r).c2=chi;
      (*r).c4=chi;

/******************************* direction -0 *********************************/

      iy=idn[ix][0];
      sm=in+iy;
      um=pu_gauge_f(iy,0);
      
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
      sp=in+iy;
      up=pu_gauge_f(ix,1);
      
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
      sm=in+iy;
      um=pu_gauge_f(iy,1);
      
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
      sp=in+iy;
      up=pu_gauge_f(ix,2);
      
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
      sm=in+iy;
      um=pu_gauge_f(iy,2);
      
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
      sp=in+iy;
      up=pu_gauge_f(ix,3);
      
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
      sm=in+iy;
      um=pu_gauge_f(iy,3);
      
      _vector_i_sub_f(psi,(*sm).c1,(*sm).c3);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c1,chi);
      _vector_i_add_assign_f((*r).c3,chi);

      _vector_i_add_f(psi,(*sm).c2,(*sm).c4);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c2,chi);
      _vector_i_sub_assign_f((*r).c4,chi);
      
/******************************** end of loop *********************************/

      _spinor_mul_f(*r,-0.5f,*r);
      /*
      _vector_mul_f((*r).c1,-0.5f,(*r).c1);
      _vector_mul_f((*r).c2,-0.5f,(*r).c2);
      _vector_mul_f((*r).c3,-0.5f,(*r).c3);
      _vector_mul_f((*r).c4,-0.5f,(*r).c4);
      */
      r+=1;
   }
}

/*
 * this function takes 2 spinors defined on the whole lattice
 * of size VOLUME
 */
void Dphi(float m0, suNf_spinor *out, suNf_spinor *in)
{
   int ix;
   float rho;
   suNf_spinor *r, *s;

   error((in==NULL)||(out==NULL),1,"Qphi [Qphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Qphi [Qphi.c]",
         "Input and output fields must be different");

   Dphi_(OE, out+(VOLUME/2), in);
   Dphi_(EO, out, in+(VOLUME/2));

   rho=+4.0f+m0;
   r=out;
   s=in;

/************************ loop over all lattice sites *************************/

   for (ix=0;ix<VOLUME;ix++) 
   {
      _spinor_mul_add_assign_f(*r,rho,*s);
      ++r;
      ++s;
   }

}

void g5Dphi(float m0, suNf_spinor *out, suNf_spinor *in)
{
   int ix;
   float rho;
   suNf_spinor *r, *s;

   error((in==NULL)||(out==NULL),1,"Qphi [Qphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Qphi [Qphi.c]",
         "Input and output fields must be different");

   Dphi_(OE, out+(VOLUME/2), in);
   Dphi_(EO, out, in+(VOLUME/2));
   
   rho=4.0f+m0;
   r=out;
   s=in;

/************************ loop over all lattice sites *************************/

   for (ix=0;ix<VOLUME;ix++) 
   {
      _spinor_mul_add_assign_f(*r,rho,*s);
      _spinor_g5_assign_f(*r);
      ++r;
      ++s;
   }
}


