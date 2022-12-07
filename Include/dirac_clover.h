// Dirac operators with clover term
#ifdef __cplusplus
    extern "C" {
#endif

void Cphi(double, spinor_field *, spinor_field *);
void g5Cphi(double, spinor_field *, spinor_field *);
void g5Cphi_sq(double, spinor_field *, spinor_field *);
void Cphi_eopre(double, spinor_field *, spinor_field *);
void g5Cphi_eopre(double, spinor_field *, spinor_field *);
void g5Cphi_eopre_sq(double, spinor_field *, spinor_field *);
void Cphi_diag(double, spinor_field *, spinor_field *);
void Cphi_diag_inv(double, spinor_field *, spinor_field *);

#ifdef __cplusplus
    }
#endif
