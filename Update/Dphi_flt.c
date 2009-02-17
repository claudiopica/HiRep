/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File Dphi_flt.c
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
#include "communications_flt.h"
#include "linear_algebra.h"

/*
 * the following variable is used to keep trace of
 * matrix-vector multoplication.
 * we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter=0;

/*
 * NOTE :
 * here we are making the assumption that the geometry is such that
 * all even sites are in the range [0,VOLUME/2[ and all odd sites are
 * in the range [VOLUME/2,VOLUME[
 */

void Dphi_flt_(spinor_field_flt *out, spinor_field_flt *in)
{
   int iy;
   _DECLARE_INT_ITERATOR(ix);
   suNf_flt *up,*um;
   suNf_vector_flt psi,chi;
   suNf_spinor_flt *r=0,*sp,*sm;

   error((in==NULL)||(out==NULL),1,"Dphi_flt_ [Dphi_flt.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi_flt_ [Dphi_flt.c]",
         "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
   error(out->type==&glat_even && in->type!=&glat_odd,1,"Dphi_ [Dphi.c]", "Spinors don't match!");
   error(out->type==&glat_odd && in->type!=&glat_even,1,"Dphi_ [Dphi.c]", "Spinors don't match!");
   error(out->type==&glattice && in->type!=&glattice,1,"Dphi_ [Dphi.c]", "Spinors don't match!");
#endif /* CHECK_SPINOR_MATCHING */

   ++MVMcounter; /* count matrix call */
   if(out->type==&glattice) ++MVMcounter;
 
/************************ loop over all lattice sites *************************/
   /* start communication of input spinor field */
   start_sf_sendrecv_flt(in);
   
   _PIECE_FOR(out->type,ix) {
     _SITE_FOR(out->type,ix) {
       r=_FIELD_AT(out,ix);
 
/******************************* direction +0 *********************************/
       
       iy=iup(ix,0);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f_flt(ix,0);
       
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
       sm=_FIELD_AT(in,iy);
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
       sp=_FIELD_AT(in,iy);
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
       sm=_FIELD_AT(in,iy);
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
       sp=_FIELD_AT(in,iy);
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
       sm=_FIELD_AT(in,iy);
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
       sp=_FIELD_AT(in,iy);
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
       sm=_FIELD_AT(in,iy);
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

       _spinor_mul_f(*r,-0.5f,*r);
     } /* SITE_FOR */
     if(_PIECE_INDEX(ix)==0) {
       /* wait for spinor to be transfered */
       complete_sf_sendrecv_flt(in);
     }
   } /* PIECE FOR */
}


/*
 * this function takes 2 spinors defined on the whole lattice
 * of size VOLUME
 */
void Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in)
{
   double rho;

   error((in==NULL)||(out==NULL),1,"Dphi_flt [Dphi_flt.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi_flt [Dphi_flt.c]",
         "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice,1,"Dphi [Dphi.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

   Dphi_flt_(out, in);

   rho=+4.0f+m0;
   spinor_field_mul_add_assign_f_flt(out,rho,in);

}

void g5Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in)
{
   double rho;

   error((in==NULL)||(out==NULL),1,"g5Dphi_flt [Dphi_flt.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"g5Dphi_flt [Dphi_flt.c]",
         "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice,1,"g5Dphi [Dphi.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

   Dphi_flt_(out, in);
   
   rho=4.0f+m0;

   spinor_field_mul_add_assign_f_flt(out,rho,in);
   spinor_field_g5_assign_f_flt(out);
}

