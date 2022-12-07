// Twisted mass operators
#ifdef __cplusplus
  extern "C" {
#endif

void Qhat_eopre(double m0, double mu, spinor_field *out, spinor_field *in);
void Qhat_eopre_sq(double m0, double mu, spinor_field *out, spinor_field *in);

void set_twisted_mass(double mass); // this is the twisted mass used in the twisted mass operators (Q)

void Qtm_p(spinor_field *out, spinor_field *in);
void Qtm_m(spinor_field *out, spinor_field *in);
void QpQm_tm(spinor_field *out, spinor_field *in);

void Qtm_p_alt(spinor_field *out, spinor_field *in);
void Qtm_m_alt(spinor_field *out, spinor_field *in);
void QpQm_tm_alt(spinor_field *out, spinor_field *in);

// Wrapper to invert the twisted mass operator Qtm_p
void tm_invert(spinor_field *out, spinor_field *in, mshift_par *mpar);
void tm_invert_alt(spinor_field *out, spinor_field *in, mshift_par *mpar);

void H(spinor_field *out, spinor_field *in);
void H_flt(spinor_field_flt *out, spinor_field_flt *in);
void H2(spinor_field *out, spinor_field *in);
void D(spinor_field *out, spinor_field *in);
void D_flt(spinor_field_flt *out, spinor_field_flt *in);

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

