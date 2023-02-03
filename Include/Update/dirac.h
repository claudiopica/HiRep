/***************************************************************************\
* Copyright (c) 2008-2022, Claudio Pica, Sofie Martins                      *
* All rights reserved.                                                      *
\***************************************************************************/

/// Headerfile for:
/// - Dphi.c
/// - Dphi_flt.c
/// - Dphi_gpu.cu
/// - D_update.c
/// - D_ff.c

#ifndef DIRAC_H
#define DIRAC_H

#include "spinor_field.h"
#include "Inverters/linear_solvers.h"
#include "Utils/generics.h"

#ifdef __cplusplus
    extern "C" {
#endif

 enum DIRECTION {
      UP = 0,
      DOWN = 1
};

#define _SPINOR_FIELD_TYPE spinor_field
#define _SUFFIX 
#include "TMPL/dirac.h.tmpl"

#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SUFFIX _flt
#include "TMPL/dirac.h.tmpl"

// Dirac operator with twisted mass
void Qhat_eopre(double m0, double mu, spinor_field *out, spinor_field *in);
void Qhat_eopre_sq(double m0, double mu, spinor_field *out, spinor_field *in);

typedef enum {
  DIRECT,
  DAGGER
} tw_D_type;

void Dxx_tw_inv(double mass, double twmass, spinor_field *out, spinor_field *in,tw_D_type tw_type);
void g5Dphi_eopre_tw(double m0, double mu, spinor_field *out, spinor_field *in,tw_D_type tw_type);
void g5Dphi_eopre_tw_sq(double m0, double mu, spinor_field *out, spinor_field *in);

#if (NG==3) && defined(REPR_FUNDAMENTAL)
void Dphi_fused_(spinor_field *out, spinor_field *in); //TODO: should we remove this?
#endif

// Clover operators
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
void Cphi(double mass, spinor_field *dptr, spinor_field *sptr);
void g5Cphi(double mass, spinor_field *dptr, spinor_field *sptr);
void g5Cphi_sq(double mass, spinor_field *dptr, spinor_field *sptr);
void Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr);
void g5Cphi_eopre(double mass, spinor_field *dptr, spinor_field *sptr);
void g5Cphi_eopre_sq(double mass, spinor_field *dptr, spinor_field *sptr);
void Cphi_diag(double mass, spinor_field *dptr, spinor_field *sptr);
void Cphi_diag_inv(double mass, spinor_field *dptr, spinor_field *sptr);
#endif

// D_update.c
void set_dirac_mass(double mass); // this is the mass used in the following operators
double get_dirac_mass(void);
void set_twisted_mass(double mu); // this is the twisted mass used in the twisted mass operators (Q)
void Hoo(spinor_field *out, spinor_field *in);
void Hoo2(spinor_field *out, spinor_field *in);
void H2(spinor_field *out, spinor_field *in);
void H(spinor_field *out, spinor_field *in);
void H_flt(spinor_field_flt *out, spinor_field_flt *in);
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

// D_ff.c
/* Dirac operators with a four fermion interaction */
void Dphi_eopre_4f(double m0, spinor_field *out, spinor_field *in, double shift);
void Dphi_eopre_4f_dagger(double m0, spinor_field *out, spinor_field *in, double shift);
void Dphieopre_4f_sq(double m0, spinor_field *out, spinor_field *in, double shift);

void Dphi_eopre_4f_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
// void Dphi_eopre_4f_dagger_flt(double m0, spinor_field_flt *out, spinor_field_flt *in); //TODO: this is not defined
// void Dphieopre_4f_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in); //TODO: this is not defined 
// void Dphieopre_4f_DDdagger(double m0, spinor_field *out, spinor_field *in, double shift); //TODO: this is defined but it was not in the header

/* Dirac operators used in update */
void set_ff_dirac_mass(double mass);  // this is the mass used in the following operators
void set_ff_dirac_shift(double mass); // The shift added to four fermion Hasenbush-Dirac operators (Dff, Dff_dagger and Df_sq)

void Dff(spinor_field *out, spinor_field *in);
void Dff_dagger(spinor_field *out, spinor_field *in);
void Dff_sq(spinor_field *out, spinor_field *in);

// void Dphi_4f(double m0, spinor_field *out, spinor_field *in); //TODO: this is defined but it was not in the header
// void Dphi_4f_dagger(double m0, spinor_field *out, spinor_field *in); //TODO: this is defined but it was not in the header
// void Dphi_4f_sq(double m0, spinor_field *out, spinor_field *in); //TODO: this is defined but it was not in the header

#ifdef __cplusplus
    }
#endif
#endif
