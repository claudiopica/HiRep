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
#include "linear_algebra.h"
#include "memory.h"

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
 *
 * D_[EO/OE] in = -0.5 *[ (1-g\mu) U in + (1+g\mu) U^+ in ]
 */

/* prende 2 spinor lunghi VOLUME/2 definiti solo su siti di parita' definita */
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

	 ++MVMcounter; /* count matrix call */

   r=out;
  

/************************ loop over all lattice sites *************************/

   for (ix=smin;ix<smax;++ix) 
   {

/******************************* direction +0 *********************************/

      iy=iup(ix,0);
      sp=in+iy;
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
      sm=in+iy;
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
      sp=in+iy;
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
      sm=in+iy;
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
      sp=in+iy;
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
      sm=in+iy;
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
      sp=in+iy;
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
      sm=in+iy;
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

      _spinor_mul_f(*r,-0.5,*r);
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
 * Dphi in = (4+m0)*in + D_ in
 */
void Dphi(double m0, suNf_spinor *out, suNf_spinor *in)
{
   int ix;
   double rho;
   suNf_spinor *r, *s;

   error((in==NULL)||(out==NULL),1,"Dphi [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi [Dphi.c]",
         "Input and output fields must be different");

   Dphi_(OE, out+(VOLUME/2), in);
   Dphi_(EO, out, in+(VOLUME/2));

   rho=+4.0+m0;
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

void g5Dphi(double m0, suNf_spinor *out, suNf_spinor *in)
{
   int ix;
   double rho;
   suNf_spinor *r, *s;

   error((in==NULL)||(out==NULL),1,"g5Dphi [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"g5Dphi [Dphi.c]",
         "Input and output fields must be different");

   Dphi_(OE, out+(VOLUME/2), in);
   Dphi_(EO, out, in+(VOLUME/2));
   
   rho=4.0+m0;
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

static suNf_spinor *stmp_pre=NULL;

void free_tmp() {
	if(stmp_pre!=NULL) {
		free_field(stmp_pre);
		stmp_pre=NULL;
	}
}

/*
 * this function takes 2 spinors defined on the even lattice
 * of size VOLUME/2
 * Dphi in = (4+m0)^2*in - D_EO D_OE in
 *
 */
void Dphi_eopre(double m0, suNf_spinor *out, suNf_spinor *in)
{
   int ix;
   double rho;
   suNf_spinor *r, *s;

   error((in==NULL)||(out==NULL),1,"Dphi_eopre [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi_eopre [Dphi.c]",
         "Input and output fields must be different");

	 /* alloc memory for temporary spinor field */
	 if (stmp_pre==NULL) {
		 unsigned int slen;
		 get_spinor_len(&slen);
		 set_spinor_len(VOLUME/2);
		 stmp_pre=alloc_spinor_field_f(1);
		 error(stmp_pre==NULL,1,"Dphi_eopre [Dphi.c]",
				 "cannot allocate memory space for temporary storage");
		 set_spinor_len(slen);
		 atexit(&free_tmp);
	 }

   Dphi_(OE, stmp_pre, in);
   Dphi_(EO, out, stmp_pre);

   rho=+4.0+m0;
	 rho*=-rho; /* this minus sign is taken into account below */

/************************ loop over all lattice sites *************************/

	 r=out;
	 s=in;
   for (ix=0;ix<(VOLUME/2);ix++) 
   {
      _spinor_mul_add_assign_f(*r,rho,*s);
			_spinor_minus_f(*r,*r);
      ++r;
      ++s;
   }

}

/*
 * this function takes 2 spinors defined on the even lattice
 * of size VOLUME/2
 * Dphi in = (4+m0)^2*in - D_EO D_OE in
 *
 */
void g5Dphi_eopre(double m0, suNf_spinor *out, suNf_spinor *in)
{
   int ix;
   double rho;
   suNf_spinor *r, *s;

   error((in==NULL)||(out==NULL),1,"Dphi_eopre [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi_eopre [Dphi.c]",
         "Input and output fields must be different");

	 /* alloc memory for temporary spinor field */
	 if (stmp_pre==NULL) {
		 unsigned int slen;
		 get_spinor_len(&slen);
		 set_spinor_len(VOLUME/2);
		 stmp_pre=alloc_spinor_field_f(1);
		 error(stmp_pre==NULL,1,"Dphi_eopre [Dphi.c]",
				 "cannot allocate memory space for temporary storage");
		 set_spinor_len(slen);
		 atexit(&free_tmp);
	 }

   Dphi_(OE, stmp_pre, in);
   Dphi_(EO, out, stmp_pre);

   rho=+4.0+m0;
	 rho*=-rho; /* this minus sign is taken into account below */

/************************ loop over all lattice sites *************************/

	 r=out;
	 s=in;
   for (ix=0;ix<(VOLUME/2);ix++) 
   {
		 _spinor_mul_add_assign_f(*r,rho,*s);
			_spinor_minus_f(*r,*r);
      _spinor_g5_assign_f(*r);
      ++r;
      ++s;
   }

}
