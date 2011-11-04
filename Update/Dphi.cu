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

/*
 * This function defines the massless Dirac operator
 * It can act on spinors defined on the whole lattice 
 * or on spinors with definite parity
 */


__global__ void Dphi_gpu(suNf_spinor* out, suNf_spinor* in, suNf* gauge, int N){
  suNf_spinor r,sn;
  SuNf_vector psi,chi;
  suNf u;
  int iy;
  int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;

  if (ix < N){
      /******************************* direction +0 *********************************/
      iy=iup(ix,0);
      u = gauge[coord_to_index(ix,0)];
      sn= in[iy];
      
      _vector_add_f(psi,sn.c[0],sn.c[2]);
      _suNf_multiply(chi,u,psi);
      
      r.c[0]=chi;
      r.c[2]=chi;

      _vector_add_f(psi,sn.c[1],sn.c[3]);
      _suNf_multiply(chi,u,psi);
      
       r.c[1]=chi;
       r.c[3]=chi;

       /******************************* direction -0 *********************************/

       iy=idn(ix,0);
       u = gauge[coord_to_index(iy,0)];
       sn = in[iy];

       _vector_sub_f(psi,sn.c[0],sn.c[2]);
       _suNf_inverse_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[0],chi);
       _vector_sub_assign_f(r.c[2],chi);

       _vector_sub_f(psi,sn.c[1],sn.c[3]);
       _suNf_inverse_multiply(chi,u,psi);
      
       _vector_add_assign_f(r.c[1],chi);
       _vector_sub_assign_f(r.c[3],chi);       

      /******************************* direction +1 *********************************/

       iy=iun(ix,1);
       u = gauge[coord_to_index(iy,1)];
       sn = in[iy];
      
       _vector_i_add_f(psi,sn.c[0],sn.c[3]);
       _suNf_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[0],chi);
       _vector_i_sub_assign_f(r.c[3],chi);

       _vector_i_add_f(psi,sn.c[1],sn.c[2]);
       _suNf_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[1],chi);
       _vector_i_sub_assign_f(r.c[2],chi);

       /******************************* direction -1 *********************************/

       iy=idn(ix,1);
       u = gauge[coord_to_index(iy,1)];
       sn = in[iy];

       _vector_i_sub_f(psi,sn.c[0],sn.c[3]);
       _suNf_inverse_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[0],chi);
       _vector_i_add_assign_f(r.c[3],chi);

       _vector_i_sub_f(psi,sn.c[1],sn.c[2]);
       _suNf_inverse_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[1],chi);
       _vector_i_add_assign_f(r.c[2],chi);

       /******************************* direction +2 *********************************/
       iy= iup(ix,2);
       u = gauge[coord_to_index(iy,2)];
       sn = in[iy];

       _vector_add_f(psi,sn.c[0],sn.c[3]);
       _suNf_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[0],chi);
       _vector_add_assign_f(r.c[3],chi);

       _vector_sub_f(psi,sn.c[1],sn.c[2]);
       _suNf_multiply(chi,u,psi);
      
       _vector_add_assign_f(r.c[1],chi);
       _vector_sub_assign_f(r.c[2],chi);

       /******************************* direction -2 *********************************/
       iy = idn(ix,2);
       u  = gauge[coord_to_index(iy,2)];
       sn = in[iy];

       _vector_sub_f(psi,sn.c[0],sn.c[3]);
       _suNf_inverse_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[0],chi);
       _vector_sub_assign_f(r.c[3],chi);

       _vector_add_f(psi,sn.c[1],sn.c[2]);
       _suNf_inverse_multiply(chi,u,psi);
      
       _vector_add_assign_f(r.c[1],chi);
       _vector_add_assign_f(r.c[2],chi);

       /******************************* direction +3 *********************************/
       iy = iup(ix,3);
       u  = gauge[coord_to_index(iy,3)];
       sn = in[iy];
      
       _vector_i_add_f(psi,sn.c[0],sn.c[2]);
       _suNf_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[0],chi);
       _vector_i_sub_assign_f(r.c[2],chi);

       _vector_i_sub_f(psi,sn.c[1],sn.c[3]);
       _suNf_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[1],chi);
       _vector_i_add_assign_f(r.c[3],chi);

       /******************************* direction -3 *********************************/

       iy = idn(ix,3);
       u  = gauge[coord_to_index(iy,3)];
       sn = in[iy];
      
       _vector_i_sub_f(psi,sn.c[0],sn.c[2]);
       _suNf_inverse_multiply(chi,u,psi);
      
       _vector_add_assign_f(r.c[0],chi);
       _vector_i_add_assign_f(r.c[2],chi);

       _vector_i_add_f(psi,sn.c[1],sn.c[3]);
       _suNf_inverse_multiply(chi,u,psi);

       _vector_add_assign_f(r.c[1],chi);
       _vector_i_sub_assign_f(r.c[3],chi);
       /******************************** end of directions *********************************/      

       _spinor_mul_f(r,-0.5,r);
    }
 }
 


/*
 * this function takes 2 spinors defined on the whole lattice
 */
 void Dphi(double m0, spinor_field *out, spinor_field *in)
{
   double rho;
   int N,grid;

   error((in==NULL)||(out==NULL),1,"Dphi [Dphi.cu]",
         "Attempt to access unallocated memory space");
   
   error(in==out,1,"Dphi [Dphi.cu]",
         "Input and output fields must be different");


#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice,1,"Dphi [Dphi.cu]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */
   
   N = out->type->master_end[0] -out->type->master_start[0] + 1 ;
   grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
   
   Dphi_gpu<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out),START_SP_ADDRESS_GPU(in), u_gauge->gpu_ptr,N);
   rho=4.+m0;
   spinor_field_mul_add_assign_f(out,rho,in);
   
}

void g5Dphi(double m0, spinor_field *out, spinor_field *in)
{
  double rho;
  int N,grid;

  error((in==NULL)||(out==NULL),1,"g5Dphi [Dphi.cu]",
	"Attempt to access unallocated memory space");

  error(in==out,1,"g5Dphi [Dphi.cu]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice,1,"g5Dphi [Dphi.cu]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

   N = out->master_end[0] -out->type->master_start[0] + 1 ;
   grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
   
   Dphi_gpu<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out), START_SP_ADDRESS(in),u_gauge->gpu_ptr,N);

   rho=4.+m0;
   spinor_field_mul_add_assign_f(out,rho,in);
   spinor_field_g5_assign_f(out);
}

static int init=1;
static spinor_field *gtmp=NULL;
static spinor_field *etmp=NULL;
static spinor_field *otmp=NULL;

static void free_mem() {
  if (gtmp!=NULL) { 
    free_spinor_field_gpu(gtmp);
    free_spinor_field(gtmp); 
    etmp=NULL; 
  }
  if (etmp!=NULL) { 
    free_spinor_field_gpu(etmp); 
    free_spinor_field(etmp); 
    etmp=NULL; 
  }
  if (otmp!=NULL) { 
    free_spinor_field_gpu(otmp);
    free_spinor_field(otmp); 
    otmp=NULL; 
  }
  init=1;
}

static void init_Dirac() {
  if (init) {
    gtmp=alloc_spinor_field_f(1,&glattice);
    alloc_spinor_field_f_gpu(1,gtmp);
    etmp=alloc_spinor_field_f(1,&glat_even);
    alloc_spinor_field_f_gpu(1,etmp);
    otmp=alloc_spinor_field_f(1,&glat_odd);
    alloc_spinor_field_f_gpu(1,otmp);
    atexit(&free_mem);
    init=0;
  }
}

/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the even lattice
 * Dphi in = (4+m0)^2*in - D_EO D_OE in
 *
 */
void Dphi_eopre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;
  int N,grid;
  
  error((in==NULL)||(out==NULL),1,"Dphi_eopre [Dphi.cu]",
	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Dphi_eopre [Dphi.cu]",
	"Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even,1,"Dphi_eopre " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); }
  
  start1 = out->type->master_start[0];
  start2 = otmp->type->master_start[0];
  N = out->type->master_end[0] -out->type->master_start[0] + 1 ;
  grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
  
  Dphi_gpu<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(otmp), START_SP_ADDRESS_GPU(in), u_gauge->gpu_ptr,N);
  Dphi_gpu<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out), START(SP_ADDRESS_GPU(otmp), u_gauge->gpu_ptr,N);
  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */
  
  spinor_field_mul_add_assign_f(out,rho,in);
  spinor_field_minus_f(out,out);
}


/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the odd lattice
 * Dphi in = (4+m0)^2*in - D_OE D_EO in
 *
 */
void Dphi_oepre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;
  int N,grid;
  
  error((in==NULL)||(out==NULL),1,"Dphi_oepre [Dphi.cu]",
	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Dphi_oepre [Dphi.cu]",
	"Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_odd || in->type!=&glat_odd,1,"Dphi_oepre " __FILE__, "Spinors are not defined on odd lattice!");
#endif /* CHECK_SPINOR_MATCHING */


  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac();}
  
  N = out->type->master_end[0] -out->type->master_start[0] + 1 ;
  grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
  
  Dphi_gpu<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(etmp), START_SP_ADDRESS_GPU(in), u_gauge->gpu_ptr,N);
  Dphi_gpu<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out), START_SP_ADDRESS_GPU(etmp), u_gauge->gpu_ptr,N);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */
  
  spinor_field_mul_add_assign_f(out,rho,in);
  spinor_field_minus_f(out,out);


}



void g5Dphi_eopre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;
  int N,grid;

  error((in==NULL)||(out==NULL),1,"g5Dphi_eopre [Dphi.cu]",
	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Dphi_eopre [Dphi.cu]",
	"Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even,1,"g5Dphi_eopre " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

#if defined(BASIC_SF) || defined(ROTATED_SF)
  SF_spinor_bcs(in);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac();}
  
  N = out->type->master_end[0] -out->type->master_start[0] + 1 ;
  grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
  
  Dphi_gpu<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(otmp), START_SP_ADDRESS_GPU(in), u_gauge->gpu_ptr,N);
  Dphi_gpu<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out), START_SP_ADDRESS_GPU(otmp), u_gauge->gpu_ptr,N);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */
  
  spinor_field_mul_add_assign_f(out,rho,in);
  spinor_field_minus_f(out,out);
  spinor_field_g5_assign_f(out);

#if defined(BASIC_SF) || defined(ROTATED_SF)
  SF_spinor_bcs(out);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */
  
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
  if (init) { init_Dirac();  }
  

  g5Dphi(m0, gtmp, in);
  g5Dphi(m0, out, gtmp);

}
