/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef DIRAC_H
#define DIRAC_H

#include "suN_types.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

void Dphi_cpu_(spinor_field *out, spinor_field *in);
void Dphi_cpu(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_cpu(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_sq_cpu(double m0, spinor_field *out, spinor_field *in);
#ifdef WITH_GPU
  void Dphi_gpu_(spinor_field *out, spinor_field *in);
  void Dphi_gpu(double m0, spinor_field *out, spinor_field *in);
  void g5Dphi_gpu(double m0, spinor_field *out, spinor_field *in);
  void g5Dphi_sq_gpu(double m0, spinor_field *out, spinor_field *in);
#endif
void (*Dphi_) (spinor_field *out, spinor_field *in);
void (*Dphi) (double m0, spinor_field *out, spinor_field *in);
void (*g5Dphi) (double m0, spinor_field *out, spinor_field *in);
void (*g5Dphi_sq) (double m0, spinor_field *out, spinor_field *in);

void Dphi_flt_cpu_(spinor_field_flt *out, spinor_field_flt *in);
void Dphi_flt_cpu(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_flt_cpu(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_sq_flt_cpu(double m0, spinor_field_flt *out, spinor_field_flt *in);
#ifdef WITH_GPU
  void Dphi_flt_gpu_(spinor_field_flt *out, spinor_field_flt *in);
  void Dphi_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in);
  void g5Dphi_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in);
  void g5Dphi_sq_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in);
#endif
void Dphi_flt_(spinor_field_flt *out, spinor_field_flt *in);
void Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

unsigned long int getMVM_cpu();
#ifdef WITH_GPU
  unsigned long int getMVM_gpu();
#endif
unsigned long int (*getMVM) ();
unsigned long int getMVM_flt();

// Dirac operators with clover term
void Cphi(double, spinor_field *, spinor_field *);
void g5Cphi(double, spinor_field *, spinor_field *);
void g5Cphi_sq(double, spinor_field *, spinor_field *);
void Cphi_eopre(double, spinor_field *, spinor_field *);
void g5Cphi_eopre(double, spinor_field *, spinor_field *);
void g5Cphi_eopre_sq(double, spinor_field *, spinor_field *);
void Cphi_diag(double, spinor_field *, spinor_field *);
void Cphi_diag_inv(double, spinor_field *, spinor_field *);

/* Even/Odd preconditioned matrix */
void Dphi_eopre_cpu(double m0, spinor_field *out, spinor_field *in);
void Dphi_oepre_cpu(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_eopre_cpu(double m0, spinor_field *out, spinor_field *in);
void g5Dphi_eopre_sq_cpu(double m0, spinor_field *out, spinor_field *in);
#ifdef WITH_GPU
  void Dphi_eopre_gpu(double m0, spinor_field *out, spinor_field *in);
  void Dphi_oepre_gpu(double m0, spinor_field *out, spinor_field *in);
  void g5Dphi_eopre_gpu(double m0, spinor_field *out, spinor_field *in);
  void g5Dphi_eopre_sq_gpu(double m0, spinor_field *out, spinor_field *in);
#endif //WITH_GPU
void (*Dphi_eopre) (double m0, spinor_field *out, spinor_field *in);
void (*Dphi_oepre) (double m0, spinor_field *out, spinor_field *in);
void (*g5Dphi_eopre) (double m0, spinor_field *out, spinor_field *in);
void (*g5Dphi_eopre_sq) (double m0, spinor_field *out, spinor_field *in);
void Dphi_eopre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphi_oepre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_eopre_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void g5Dphi_eopre_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

// Twisted mass operators
void Qhat_eopre(double m0, double mu, spinor_field *out, spinor_field *in);
void Qhat_eopre_sq(double m0, double mu, spinor_field *out, spinor_field *in);

/* Dirac operators used in the Update */
void set_dirac_mass(double mass); // this is the mass used in the following operators
double get_dirac_mass();
void set_twisted_mass(double mass); // this is the twisted mass used in the twisted mass operators (Q)
void H(spinor_field *out, spinor_field *in);
void H_flt(spinor_field_flt *out, spinor_field_flt *in);
void H2(spinor_field *out, spinor_field *in);
void D(spinor_field *out, spinor_field *in);
void D_flt(spinor_field_flt *out, spinor_field_flt *in);

void Qtm_p(spinor_field *out, spinor_field *in);
void Qtm_m(spinor_field *out, spinor_field *in);
void QpQm_tm(spinor_field *out, spinor_field *in);

void Qtm_p_alt(spinor_field *out, spinor_field *in);
void Qtm_m_alt(spinor_field *out, spinor_field *in);
void QpQm_tm_alt(spinor_field *out, spinor_field *in);

// Wrapper to invert the twisted mass operator Qtm_p
void tm_invert(spinor_field *out, spinor_field *in, mshift_par *mpar);
void tm_invert_alt(spinor_field *out, spinor_field *in, mshift_par *mpar);

/* Dirac operators with a four fermion interaction */
void Dphi_eopre_4f(double m0, spinor_field *out, spinor_field *in, double shift);
void Dphi_eopre_4f_dagger(double m0, spinor_field *out, spinor_field *in, double shift);
void Dphieopre_4f_sq(double m0, spinor_field *out, spinor_field *in, double shift);

void Dphi_eopre_4f_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphi_eopre_4f_dagger_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphieopre_4f_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

/* Dirac operators used in update */
void set_ff_dirac_mass(double mass);  // this is the mass used in the following operators
void set_ff_dirac_shift(double mass); // The shift added to four fermion Hasenbush-Dirac operators (Dff, Dff_dagger and Df_sq)
void Dff(spinor_field *out, spinor_field *in);
void Dff_dagger(spinor_field *out, spinor_field *in);
void Dff_sq(spinor_field *out, spinor_field *in);

typedef enum
{
  DIRECT,
  DAGGER
} tw_D_type;

void Dxx_tw_inv(double mass, double twmass, spinor_field *out, spinor_field *in, tw_D_type tw_type);
void g5Dphi_eopre_tw(double m0, double mu, spinor_field *out, spinor_field *in, tw_D_type tw_type);
void g5Dphi_eopre_tw_sq(double m0, double mu, spinor_field *out, spinor_field *in);
#ifdef __cplusplus
}
#endif

#endif
