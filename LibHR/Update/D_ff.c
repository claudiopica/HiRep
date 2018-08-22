/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Jarno Rantaharju                        *   
* All rights reserved.                                                      * 
\***************************************************************************/


#include "stddef.h"
#include "logger.h"
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



void spinor_sigma_pi_rho_div_assign(spinor_field *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field *in);
void spinor_scalarfield_mult_add_assign(spinor_field *out,scalar_field *sigma,double rho, spinor_field *in);
void spinor_scalarfield_ig5_mult_add_assign(spinor_field *out,scalar_field *pi, spinor_field *in);
void spinor_sigma_pi_dagger_rho_div_assign(spinor_field *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field *in);
void spinor_scalarfield_mig5_mult_add_assign(spinor_field *out,scalar_field *pi, spinor_field *in);

static double static_mass=0.;
static double static_shift=0.;

void set_ff_dirac_mass(double mass) {
  static_mass=mass;
}

void set_ff_dirac_shift(double shift){
  static_shift=shift;
}


//Fields for temporary storage, from Dphi.c
extern int init_dirac;
extern spinor_field *gtmp;
extern spinor_field *etmp;
extern spinor_field *otmp;
extern spinor_field *otmp2;

void init_Dirac();




//NOTE: With four fermion auxiliary fields, there is a diagonal
//      term: A = sigma + i*g5*pi.
//      The even odd preconditioned matrix is
//      Dphi in = (4+m0+A_E)*in - D_EO 1/(4+m_0+A_O) D_OE in
//
//      In addition there is the A_0 term in the action,
//      which is independent of fermion fields, S_O = log Tr (A_O^2)
void Dphi_eopre_4f(double m0, spinor_field *out, spinor_field *in, double shift)
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
  if (init_dirac) { init_Dirac(); init_dirac=0; }
  
  rho=4.0+m0;

  Dphi_(otmp, in);
  spinor_sigma_pi_rho_div_assign(otmp,ff_sigma,ff_pi,rho, otmp);
  apply_BCs_on_spinor_field(otmp);
  Dphi_(out, otmp);
  spinor_field_minus_f(out,out);
  
  spinor_scalarfield_mult_add_assign(out,ff_sigma,rho,in);
  spinor_scalarfield_ig5_mult_add_assign(out,ff_pi,in);

  spinor_field_mul_add_assign_f(out,shift,in);

  apply_BCs_on_spinor_field(out);
}


//With pi!=0, this is not hermitiean even with after multiplication with g5
//We need to invert D^\dagger * D
//Define d^\dagger here, using g5 on the original D
void Dphi_eopre_4f_dagger(double m0, spinor_field *out, spinor_field *in, double shift)
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
  if (init_dirac) { init_Dirac(); init_dirac=0; }
  
  rho=4.0+m0;

  //Get dagger with g5, for now
  spinor_field_g5_assign_f(in);
  Dphi_(otmp, in);
  spinor_field_g5_assign_f(in);
  spinor_field_g5_assign_f(otmp);
  apply_BCs_on_spinor_field(otmp);
  spinor_sigma_pi_dagger_rho_div_assign(otmp,ff_sigma,ff_pi,rho, otmp);
  spinor_field_g5_assign_f(otmp);
  Dphi_(out, otmp);
  spinor_field_g5_assign_f(out);
  spinor_field_minus_f(out,out);
  

  spinor_scalarfield_mult_add_assign(out,ff_sigma,rho,in);
  spinor_scalarfield_mig5_mult_add_assign(out,ff_pi,in);

  spinor_field_mul_add_assign_f(out,shift,in);

  apply_BCs_on_spinor_field(out);
}


/* Dphi_4f^dagger * Dphi_4f */
void Dphieopre_4f_sq(double m0, spinor_field *out, spinor_field *in, double shift) {
  /* alloc memory for temporary spinor field */
  if (init_dirac) { init_Dirac(); init_dirac=0; }

  Dphi_eopre_4f(m0, etmp, in, shift);
  Dphi_eopre_4f_dagger(m0, out, etmp, shift);

  
}
void Dphieopre_4f_DDdagger(double m0, spinor_field *out, spinor_field *in, double shift) {
  /* alloc memory for temporary spinor field */
  if (init_dirac) { init_Dirac(); init_dirac=0; }
  Dphi_eopre_4f_dagger(m0, etmp, in, shift);
  Dphi_eopre_4f(m0, out, etmp, shift);
}


void Dphi_4f(double m0, spinor_field *out, spinor_field *in)
{
   double rho;
  
  error((in==NULL)||(out==NULL),1,"Dphi_eopre [Dphi.c]",
	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Dphi_eopre [Dphi.c]",
	"Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glattice || in->type!=&glattice,1,"Dphi_eopre " __FILE__, "Spinors are not defined on full lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field(in);

  /* alloc memory for temporary spinor field */
  if (init_dirac) { init_Dirac(); init_dirac=0; }
  
  rho=4.0+m0;

  Dphi_(out, in);

  spinor_scalarfield_mult_add_assign(out,ff_sigma,rho,in);
  spinor_scalarfield_ig5_mult_add_assign(out,ff_pi,in);
  apply_BCs_on_spinor_field(out);
}

void Dphi_4f_dagger(double m0, spinor_field *out, spinor_field *in)
{
   double rho;
  
  error((in==NULL)||(out==NULL),1,"Dphi_eopre [Dphi.c]",
	"Attempt to access unallocated memory space");
  
  error(in==out,1,"Dphi_eopre [Dphi.c]",
	"Input and output fields must be different");
  
#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glattice || in->type!=&glattice,1,"Dphi_eopre " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  apply_BCs_on_spinor_field(in);

  /* alloc memory for temporary spinor field */
  if (init_dirac) { init_Dirac(); init_dirac=0; }
  
  rho=4.0+m0;

  //Get dagger with g5, for now
  spinor_field_g5_assign_f(in);
  Dphi_(out, in);
  spinor_field_g5_assign_f(in);
  spinor_field_g5_assign_f(out);
  apply_BCs_on_spinor_field(otmp);  

  spinor_scalarfield_mult_add_assign(out,ff_sigma,rho,in);
  spinor_scalarfield_mig5_mult_add_assign(out,ff_pi,in);
  apply_BCs_on_spinor_field(out);
}

void Dphi_4f_sq(double m0, spinor_field *out, spinor_field *in) {
  /* alloc memory for temporary spinor field */
  if (init_dirac) { init_Dirac(); init_dirac=0; }

  Dphi_4f(m0, gtmp, in);
  Dphi_4f_dagger(m0, out, gtmp);
}





//Dirac operators with four fermion term
void Dff(spinor_field *out, spinor_field *in){
  //printf(" %g %g\n", static_mass, static_shift);
#ifdef UPDATE_EO
  Dphi_eopre_4f(static_mass, out, in, static_shift);
#else
  // Not implemented
#endif
}

void Dff_dagger(spinor_field *out, spinor_field *in){
  //printf("d %g %g\n", static_mass, static_shift);
#ifdef UPDATE_EO
  Dphi_eopre_4f_dagger(static_mass, out, in, static_shift);
#else
  // Not implemented
#endif
}

void Dff_sq(spinor_field *out, spinor_field *in){
  //printf("sq %g %g\n", static_mass, static_shift);
#ifdef UPDATE_EO
  Dphieopre_4f_sq(static_mass, out, in, static_shift);
#else
  // Not implemented
#endif
}




