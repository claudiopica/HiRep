/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

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
#include "spinor_field.h"
#include "geometry.h"
#include "communications.h"
#include "memory.h"

#ifdef ROTATED_SF
#include "update.h"
extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif /* ROTATED_SF */


/*
 * Init of Dphi
 */

static int init=1;
static spinor_field *gtmp=NULL;
static spinor_field *etmp=NULL;
static spinor_field *otmp=NULL;
static spinor_field *otmp2=NULL;

static void free_mem() {
    if (gtmp!=NULL) { free_spinor_field_f(gtmp); etmp=NULL; }
    if (etmp!=NULL) { free_spinor_field_f(etmp); etmp=NULL; }
    if (otmp!=NULL) { free_spinor_field_f(otmp); otmp=NULL; }
    if (otmp2!=NULL) { free_spinor_field_f(otmp2); otmp2=NULL; }
    init=1;
}

static void init_Dirac() {
    if (init) {
        gtmp=alloc_spinor_field_f(1,&glattice);
        etmp=alloc_spinor_field_f(1,&glat_even);
        otmp=alloc_spinor_field_f(1,&glat_odd);
        otmp2=alloc_spinor_field_f(1,&glat_odd);
        atexit(&free_mem);
        init=0;
    }
}


/*
 * the following variable is used to keep trace of
 * matrix-vector multiplication.
 * we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter=0;

unsigned long int getMVM() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	MVMcounter=0; /* reset counter */

	return res;
}


/* r=t*u*s */
#ifdef BC_T_THETA

#define _suNf_theta_T_multiply(r,u,s)\
    _suNf_multiply(vtmp,(u),(s));\
    _vector_mulc_f((r),eitheta[0],vtmp)
    
#define _suNf_theta_T_inverse_multiply(r,u,s)\
    _suNf_inverse_multiply(vtmp,(u),(s));\
    _vector_mulc_star_f((r),eitheta[0],vtmp)

#else

#define _suNf_theta_T_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_T_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_X_THETA

#define _suNf_theta_X_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta[1],vtmp)

#define _suNf_theta_X_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta[1],vtmp)

#else

#define _suNf_theta_X_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_X_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_Y_THETA

#define _suNf_theta_Y_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta[2],vtmp)

#define _suNf_theta_Y_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta[2],vtmp)

#else

#define _suNf_theta_Y_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_Y_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif

/* r=t*u*s */
#ifdef BC_Z_THETA

#define _suNf_theta_Z_multiply(r,u,s)\
_suNf_multiply(vtmp,(u),(s));\
_vector_mulc_f((r),eitheta[3],vtmp)

#define _suNf_theta_Z_inverse_multiply(r,u,s)\
_suNf_inverse_multiply(vtmp,(u),(s));\
_vector_mulc_star_f((r),eitheta[3],vtmp)

#else

#define _suNf_theta_Z_multiply(r,u,s) _suNf_multiply((r),(u),(s))
#define _suNf_theta_Z_inverse_multiply(r,u,s) _suNf_inverse_multiply((r),(u),(s))

#endif




/*
 * This function defines the massless Dirac operator
 * It can act on spinors defined on the whole lattice 
 * or on spinors with definite parity
 */

void Dphi_(spinor_field *out, spinor_field *in)
{

   error((in==NULL)||(out==NULL),1,"Dphi_ [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi_ [Dphi.c]",
         "Input and output fields must be different");

#ifndef CHECK_SPINOR_MATCHING
   error(out->type==&glat_even && in->type!=&glat_odd,1,"Dphi_ [Dphi.c]", "Spinors don't match! (1)");
   error(out->type==&glat_odd && in->type!=&glat_even,1,"Dphi_ [Dphi.c]", "Spinors don't match! (2)");
   error(out->type==&glattice && in->type!=&glattice,1,"Dphi_ [Dphi.c]", "Spinors don't match! (3)");
#endif

        
   ++MVMcounter; /* count matrix calls */
   if(out->type==&glattice) ++MVMcounter;
 
/************************ loop over all lattice sites *************************/
   /* start communication of input spinor field */
   start_sf_sendrecv(in);

  _PIECE_FOR(out->type,ixp) {
     if(ixp==out->type->inner_master_pieces) {
       _OMP_PRAGMA( master )
       /* wait for spinor to be transfered */
       complete_sf_sendrecv(in);
       _OMP_PRAGMA( barrier )
     }
     _SITE_FOR(out->type,ixp,ix) {

       int iy;
       suNf *up,*um;
       suNf_vector psi,chi;
       suNf_spinor *r=0,*sp,*sm;
#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
       suNf_vector vtmp;
#endif

       r=_FIELD_AT(out,ix);
 
       /******************************* direction +0 *********************************/

       iy=iup(ix,0);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f(ix,0);
      
       _vector_add_f(psi,(*sp).c[0],(*sp).c[2]);
       _suNf_theta_T_multiply(chi,(*up),psi);
      
       (*r).c[0]=chi;
       (*r).c[2]=chi;

       _vector_add_f(psi,(*sp).c[1],(*sp).c[3]);
       _suNf_theta_T_multiply(chi,(*up),psi);
         
       (*r).c[1]=chi;
       (*r).c[3]=chi;

       /******************************* direction -0 *********************************/

       iy=idn(ix,0);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f(iy,0);
      
       _vector_sub_f(psi,(*sm).c[0],(*sm).c[2]);
       _suNf_theta_T_inverse_multiply(chi,(*um),psi);
         
       _vector_add_assign_f((*r).c[0],chi);
       _vector_sub_assign_f((*r).c[2],chi);

       _vector_sub_f(psi,(*sm).c[1],(*sm).c[3]);
       _suNf_theta_T_inverse_multiply(chi,(*um),psi);
         
       _vector_add_assign_f((*r).c[1],chi);
       _vector_sub_assign_f((*r).c[3],chi);

       /******************************* direction +1 *********************************/

       iy=iup(ix,1);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f(ix,1);
      
       _vector_i_add_f(psi,(*sp).c[0],(*sp).c[3]);
       _suNf_theta_X_multiply(chi,(*up),psi);
         
       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_sub_assign_f((*r).c[3],chi);

       _vector_i_add_f(psi,(*sp).c[1],(*sp).c[2]);
       _suNf_theta_X_multiply(chi,(*up),psi);
         
       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_sub_assign_f((*r).c[2],chi);

       /******************************* direction -1 *********************************/

       iy=idn(ix,1);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f(iy,1);
      
       _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[3]);
       _suNf_theta_X_inverse_multiply(chi,(*um),psi);
         
       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_add_assign_f((*r).c[3],chi);

       _vector_i_sub_f(psi,(*sm).c[1],(*sm).c[2]);
       _suNf_theta_X_inverse_multiply(chi,(*um),psi);
         
       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_add_assign_f((*r).c[2],chi);

       /******************************* direction +2 *********************************/

       iy=iup(ix,2);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f(ix,2);
      
       _vector_add_f(psi,(*sp).c[0],(*sp).c[3]);
       _suNf_theta_Y_multiply(chi,(*up),psi);
         
       _vector_add_assign_f((*r).c[0],chi);
       _vector_add_assign_f((*r).c[3],chi);

       _vector_sub_f(psi,(*sp).c[1],(*sp).c[2]);
       _suNf_theta_Y_multiply(chi,(*up),psi);
      
       _vector_add_assign_f((*r).c[1],chi);
       _vector_sub_assign_f((*r).c[2],chi);

       /******************************* direction -2 *********************************/

       iy=idn(ix,2);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f(iy,2);
      
       _vector_sub_f(psi,(*sm).c[0],(*sm).c[3]);
       _suNf_theta_Y_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_sub_assign_f((*r).c[3],chi);

       _vector_add_f(psi,(*sm).c[1],(*sm).c[2]);
       _suNf_theta_Y_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_add_assign_f((*r).c[2],chi);

       /******************************* direction +3 *********************************/

       iy=iup(ix,3);
       sp=_FIELD_AT(in,iy);
       up=pu_gauge_f(ix,3);
      
       _vector_i_add_f(psi,(*sp).c[0],(*sp).c[2]);
       _suNf_theta_Z_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_sub_assign_f((*r).c[2],chi);

       _vector_i_sub_f(psi,(*sp).c[1],(*sp).c[3]);
       _suNf_theta_Z_multiply(chi,(*up),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_add_assign_f((*r).c[3],chi);

       /******************************* direction -3 *********************************/

       iy=idn(ix,3);
       sm=_FIELD_AT(in,iy);
       um=pu_gauge_f(iy,3);
      
       _vector_i_sub_f(psi,(*sm).c[0],(*sm).c[2]);
       _suNf_theta_Z_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[0],chi);
       _vector_i_add_assign_f((*r).c[2],chi);

       _vector_i_add_f(psi,(*sm).c[1],(*sm).c[3]);
       _suNf_theta_Z_inverse_multiply(chi,(*um),psi);

       _vector_add_assign_f((*r).c[1],chi);
       _vector_i_sub_assign_f((*r).c[3],chi);
      
       /******************************** end of loop *********************************/

       _spinor_mul_f(*r,-0.5,*r);

     } /* SITE_FOR */
   } /* PIECE FOR */
}



/*
 * this function takes 2 spinors defined on the whole lattice
 */
void Dphi(double m0, spinor_field *out, spinor_field *in)
{
   double rho;
#ifdef ROTATED_SF
   int ix,iy,iz,index;
   suNf_spinor *r, *sp;
   double SFrho;
   suNf_spinor tmp;
#endif /* ROTATED_SF */
   
   error((in==NULL)||(out==NULL),1,"Dphi [Dphi.c]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi [Dphi.c]",
         "Input and output fields must be different");

   apply_BCs_on_spinor_field(in);

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice,1,"Dphi [Dphi.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */
    
   Dphi_(out, in);

   rho=4.+m0;
   spinor_field_mul_add_assign_f(out,rho,in);
   
#ifdef ROTATED_SF
   SFrho=3.*_update_par.SF_ds+_update_par.SF_zf-4.;

	if(COORD[0] == 0) {
		for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
			index=ipt(1,ix,iy,iz);
			r=_FIELD_AT(out,index);
      sp=_FIELD_AT(in,index);
			_spinor_mul_add_assign_f(*r,SFrho,*sp);
			
			_spinor_pminus_f(tmp,*sp);
			_spinor_g5_assign_f(tmp);
			if(_update_par.SF_sign==1) {
				_spinor_i_add_assign_f(*r,tmp);
			} else {
				_spinor_i_sub_assign_f(*r,tmp);
			}
		}
	}
	if(COORD[0] == NP_T-1) {
		for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
			index=ipt(T-1,ix,iy,iz);
			r=_FIELD_AT(out,index);
      sp=_FIELD_AT(in,index);
			_spinor_mul_add_assign_f(*r,SFrho,*sp);
			
			_spinor_pplus_f(tmp,*sp);
			_spinor_g5_assign_f(tmp);
			if(_update_par.SF_sign==1) {
				_spinor_i_add_assign_f(*r,tmp);
			} else {
				_spinor_i_sub_assign_f(*r,tmp);
			}
		}
	}
#endif /* ROTATED_SF */

  apply_BCs_on_spinor_field(out);
}

void g5Dphi(double m0, spinor_field *out, spinor_field *in)
{
   double rho;
#ifdef ROTATED_SF
   int ix,iy,iz,index;
   suNf_spinor *r, *sp;
   double SFrho;
   suNf_spinor tmp;
#endif /* ROTATED_SF */

   error((in==NULL)||(out==NULL),1,"g5Dphi [Dphi.c]",
         "Attempt to access unallocated memory space");

   error(in==out,1,"g5Dphi [Dphi.c]",
         "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice,1,"g5Dphi [Dphi.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */
    
   apply_BCs_on_spinor_field(in);

   Dphi_(out, in);

   rho=4.+m0;
   spinor_field_mul_add_assign_f(out,rho,in);
   
#ifdef ROTATED_SF
   SFrho=3.*_update_par.SF_ds+_update_par.SF_zf-4.;
   
   if(COORD[0] == 0) {
     for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
       index=ipt(1,ix,iy,iz);
       r=_FIELD_AT(out,index);
       sp=_FIELD_AT(in,index);
       _spinor_mul_add_assign_f(*r,SFrho,*sp);
       
       _spinor_pminus_f(tmp,*sp);
       _spinor_g5_assign_f(tmp);
       if(_update_par.SF_sign==1) {
	 _spinor_i_add_assign_f(*r,tmp);
       } else {
	 _spinor_i_sub_assign_f(*r,tmp);
       }
     }
   }
   if(COORD[0] == NP_T-1) {
     for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
       index=ipt(T-1,ix,iy,iz);
       r=_FIELD_AT(out,index);
       sp=_FIELD_AT(in,index);
       _spinor_mul_add_assign_f(*r,SFrho,*sp);
       
       _spinor_pplus_f(tmp,*sp);
       _spinor_g5_assign_f(tmp);
       if(_update_par.SF_sign==1) {
	 _spinor_i_add_assign_f(*r,tmp);
       } else {
	 _spinor_i_sub_assign_f(*r,tmp);
       }
     }
   }
#endif /* ROTATED_SF */
   
   spinor_field_g5_assign_f(out);

   apply_BCs_on_spinor_field(out);
}



/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the even lattice
 * Dphi in = (4+m0)^2*in - D_EO D_OE in
 *
 */
void Dphi_eopre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;
  
  error((in==NULL)||(out==NULL),1,"Dphi_eopre [Dphi.c]",
	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Dphi_eopre [Dphi.c]",
	"Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even,1,"Dphi_eopre " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field(in);

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); init=0; }
  
  Dphi_(otmp, in);
  apply_BCs_on_spinor_field(otmp);
  Dphi_(out, otmp);
  
  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */
  
  spinor_field_mul_add_assign_f(out,rho,in);
  spinor_field_minus_f(out,out);
  apply_BCs_on_spinor_field(out);
}


/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the odd lattice
 * Dphi in = (4+m0)^2*in - D_OE D_EO in
 *
 */
void Dphi_oepre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;
  
  error((in==NULL)||(out==NULL),1,"Dphi_oepre [Dphi.c]",
	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Dphi_oepre [Dphi.c]",
	"Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_odd || in->type!=&glat_odd,1,"Dphi_oepre " __FILE__, "Spinors are not defined on odd lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field(in);

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); init=0; }
  
  Dphi_(etmp, in);
  apply_BCs_on_spinor_field(etmp);
  Dphi_(out, etmp);
  
  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */
  
  spinor_field_mul_add_assign_f(out,rho,in);
  spinor_field_minus_f(out,out);

  apply_BCs_on_spinor_field(out);
}



void g5Dphi_eopre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;
  
  error((in==NULL)||(out==NULL),1,"g5Dphi_eopre [Dphi.c]",
	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Dphi_eopre [Dphi.c]",
	"Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even,1,"g5Dphi_eopre " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field(in);

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); init=0; }
  
  Dphi_(otmp, in);
  apply_BCs_on_spinor_field(otmp);
  Dphi_(out, otmp);
  
  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */
  
  spinor_field_mul_add_assign_f(out,rho,in);
  spinor_field_minus_f(out,out);
  spinor_field_g5_assign_f(out);

  apply_BCs_on_spinor_field(out);
}

/* g5Dphi_eopre ^2 */
void g5Dphi_eopre_sq(double m0, spinor_field *out, spinor_field *in) {
  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); init=0; }

  g5Dphi_eopre(m0, etmp, in);
  g5Dphi_eopre(m0, out, etmp);
  
}

/* g5Dhi ^2 */
void g5Dphi_sq(double m0, spinor_field *out, spinor_field *in) {
  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); init=0; }
  
#ifdef ROTATED_SF
  /*the switch of the SF_sign is needed to take care of the antihermiticity of the boundary term of the dirac operator*/
  g5Dphi(m0, gtmp, in);
  _update_par.SF_sign = -_update_par.SF_sign;
  g5Dphi(m0, out, gtmp);
  _update_par.SF_sign = -_update_par.SF_sign;
#else
  g5Dphi(m0, gtmp, in);
  g5Dphi(m0, out, gtmp);
#endif

}

//Twisted mass operator for even odd preconditioned case
/* g5 (M^+-_ee-M_eo {M_oo}^{-1} M_oe*/
void Qhat_eopre(double m0, double mu, spinor_field* out, spinor_field *in){
  double norm = (4+m0)*(4+m0)+mu*mu;
  double rho = (4+m0)/norm;
  complex imu;
  imu.re=0;
  imu.im=-mu/norm;
  
  error((in==NULL)||(out==NULL),1,"Qhat_eopre [Dphi.c]",
  	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Qhat_eopre [Dphi.c]",
        "Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even,1,"Qhat_eopre " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field(in);

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); init=0; }
  Dphi_(otmp, in);
  apply_BCs_on_spinor_field(otmp);
  spinor_field_mul_f(otmp2,rho,otmp);
  spinor_field_g5_mulc_add_assign_f(otmp2,imu,otmp);
  Dphi_(out, otmp2);
  
  rho = -(4+m0);
  spinor_field_mul_add_assign_f(out,rho,in);
  imu.im=-mu;
  spinor_field_g5_mulc_add_assign_f(out,imu,in);

  spinor_field_minus_f(out,out);
  spinor_field_g5_assign_f(out);

  apply_BCs_on_spinor_field(out);
}


void Qhat_eopre_sq(double m0, double mu, spinor_field *out, spinor_field *in) {
  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); init=0; }
  
#ifdef ROTATED_SF
  /*the switch of the SF_sign is needed to take care of the antihermiticity of the boundary term of the dirac operator*/
  error(1,"Qhat_eopre_sq",__FILE__,"Not implemented\n");
#else
  Qhat_eopre(m0,-mu,etmp,in);
  Qhat_eopre(m0,mu,out,etmp);
#endif
}
