#include "dirac.h"


static double static_mass=0.;

void set_dirac_mass(double mass) {
  static_mass=mass;
}

/* this is the basic operator used in the update */
void H2(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
  g5Dphi_eopre_sq(static_mass, out, in);
#else
  g5Dphi_sq(static_mass, out, in);
#endif
}

void H(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
  g5Dphi_eopre(static_mass, out, in);
#else
  g5Dphi(static_mass, out, in);
#endif
}

void H_flt(spinor_field_flt *out, spinor_field_flt *in){
#ifdef UPDATE_EO
  g5Dphi_eopre_flt(static_mass, out, in);
#else
  g5Dphi_flt(static_mass, out, in);
#endif
}

void D(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
  Dphi_eopre(static_mass, out, in);
#else
  Dphi(static_mass, out, in);
#endif
}


void D_flt(spinor_field_flt *out, spinor_field_flt *in){
#ifdef UPDATE_EO
  Dphi_eopre_flt((float)(static_mass), out, in);
#else
  Dphi_flt((float)(static_mass), out, in);
#endif
}
