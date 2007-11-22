/*******************************************************************************
*
* File Dphi.c
*
* Action of the Wilson-Dirac operator D and hermitian g5D on a given 
* double-precision spinor field
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
 * the following variable is used to keep trace of
 * matrix-vector multoplication.
 * we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter=0;

unsigned long int getMVM() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	MVMcounter=0; /* reset counter */

	return res;
}

/*
 * NOTE :
 * here we are making the assumption that the geometry is such that
 * all even sites are in the range [0,VOLUME/2[ and all odd sites are
 * in the range [VOLUME/2,VOLUME[
 */

/* prende 2 spinor lunghi VOLUME/2 definiti solo su siti con la stessa parita' */
void Dphi_(block_selector B, spinor_field *out, spinor_field *in)
{
   int ix,iy,sx,block=0;
   suNf *up,*um;
   suNf_vector psi,chi;
   suNf_spinor *r=0,*sp,*sm;

   error((in==NULL)||(out==NULL),1,"Dphi_ [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi_ [Dphi.c]",
         "Input and output fields must be different");

   switch(B) {
      case EO:
         block=EVENSITES;
         break;
      case OE:
         block=ODDSITES;
         break;
      default:
         error(1,1,"Dphi_ [Dphi.c]",
               "Invalid block parity selection");
   }

	 ++MVMcounter; /* count matrix call */

 
/************************ loop over all lattice sites *************************/

   FOR_SOME_SC(block,ix,sx)
   {
      r=_SPINOR_AT(out,sx);
 
/******************************* direction +0 *********************************/

      iy=iup(ix,0);
      sp=_SPINOR_AT_SITE(in,iy);
      up=pu_gauge_f(ix,0);
      
      _vector_add_f(psi,(*sp).c[0],(*sp).c[2]);
      _suNf_multiply(chi,(*up),psi);
      
      (*r).c[0]=chi;
      (*r).c[2]=chi;

      _vector_add_f(psi,(*sp).c[1],(*sp).c[3]);
      _suNf_multiply(chi,(*up),psi);
            
      (*r).c[1]=chi;
      (*r).c[3]=chi;

/******************************* direction -0 *********************************/

      iy=idn(ix,0);
      sm=_SPINOR_AT_SITE(in,iy);
      um=pu_gauge_f(iy,0);
      
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
      sp=_SPINOR_AT_SITE(in,iy);
      up=pu_gauge_f(ix,1);
      
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
      sm=_SPINOR_AT_SITE(in,iy);
      um=pu_gauge_f(iy,1);
      
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
      sp=_SPINOR_AT_SITE(in,iy);
      up=pu_gauge_f(ix,2);
      
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
      sm=_SPINOR_AT_SITE(in,iy);
      um=pu_gauge_f(iy,2);
      
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
      sp=_SPINOR_AT_SITE(in,iy);
      up=pu_gauge_f(ix,3);
      
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
      sm=_SPINOR_AT_SITE(in,iy);
      um=pu_gauge_f(iy,3);
      
      _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[2]);
      _suNf_inverse_multiply(chi,(*um),psi);
      
      _vector_add_assign_f((*r).c[0],chi);
      _vector_i_add_assign_f((*r).c[2],chi);

      _vector_i_add_f(psi,(*sm).c[1],(*sm).c[3]);
      _suNf_inverse_multiply(chi,(*um),psi);

      _vector_add_assign_f((*r).c[1],chi);
      _vector_i_sub_assign_f((*r).c[3],chi);
      
/******************************** end of loop *********************************/

      _spinor_mul_f(*r,-0.5f,*r);
   }
}



/*
 * this function takes 2 spinors defined on the whole lattice
 * of size VOLUME
 */
void Dphi(double m0, spinor_field *out, spinor_field *in)
{
   int ix, sx;
   double rho;
   suNf_spinor *r, *s;

   error((in==NULL)||(out==NULL),1,"Dphi [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi [Dphi.c]",
         "Input and output fields must be different");

   Dphi_(OE, out, in);
   Dphi_(EO, out, in);

   rho=+4.0f+m0;

/************************ loop over all lattice sites *************************/

   FOR_LOCAL_SC(ix,sx) {
      r=_SPINOR_AT(out,sx);
      s=_SPINOR_AT(in,sx);
      _spinor_mul_add_assign_f(*r,rho,*s);
   }

}

void g5Dphi(double m0, spinor_field *out, spinor_field *in)
{
   int ix,sx;
   double rho;
   suNf_spinor *r, *s;

   error((in==NULL)||(out==NULL),1,"g5Dphi [Dphi.c]",
         "Attempt to access unallocated memory space");

   error(in==out,1,"g5Dphi [Dphi.c]",
         "Input and output fields must be different");

   Dphi_(OE, out, in);
   Dphi_(EO, out, in);
   
   rho=4.0f+m0;

/************************ loop over all lattice sites *************************/

   FOR_LOCAL_SC(ix,sx) {
      r=_SPINOR_AT(out,sx);
      s=_SPINOR_AT(in,sx);
      _spinor_mul_add_assign_f(*r,rho,*s);
      _spinor_g5_assign_f(*r);
   }
}

