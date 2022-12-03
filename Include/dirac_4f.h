/* Dirac operators with a four fermion interaction */
void Dphi_eopre_4f(double m0, spinor_field *out, spinor_field *in, double shift);
void Dphi_eopre_4f_dagger(double m0, spinor_field *out, spinor_field *in, double shift);
void Dphieopre_4f_sq(double m0, spinor_field *out, spinor_field *in, double shift);

void Dphi_eopre_4f_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphi_eopre_4f_dagger_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);
void Dphieopre_4f_sq_flt(double m0, spinor_field_flt *out, spinor_field_flt *in);

