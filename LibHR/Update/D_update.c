#include "dirac.h"
#include "linear_algebra.h"
#include "stddef.h"
#include "memory.h"
#include "global.h"
#include "logger.h"

static double static_mass=0.;
static double static_mu=0.;


void set_dirac_mass(double mass) {
  static_mass=mass;
}

double get_dirac_mass(){
	return static_mass;
}

void set_twisted_mass(double mu){
  static_mu=mu;
}

/* this is the basic operator used in the update */
void Hoo(spinor_field *out, spinor_field *in){
#ifdef WITH_CLOVER
	Cphi_diag(static_mass, out, in);
	spinor_field_g5_assign_f(out);
#endif
}

void Hoo2(spinor_field *out, spinor_field *in){
#ifdef WITH_CLOVER
	Cphi_diag(static_mass, out, in);
	spinor_field_g5_assign_f(out);
	Cphi_diag(static_mass, out, out);
	spinor_field_g5_assign_f(out);
#endif
}

void H2(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
#ifdef WITH_CLOVER
	g5Cphi_eopre_sq(static_mass, out, in);
#else
	g5Dphi_eopre_sq(static_mass, out, in);
#endif
#else
#ifdef WITH_CLOVER
	g5Cphi_sq(static_mass, out, in);
#else
	g5Dphi_sq(static_mass, out, in);
#endif
#endif
}

void H(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
#ifdef WITH_CLOVER
	g5Cphi_eopre(static_mass, out, in);
#else
	g5Dphi_eopre(static_mass, out, in);
#endif
#else
#ifdef WITH_CLOVER
	g5Cphi(static_mass, out, in);
#else
	g5Dphi(static_mass, out, in);
#endif
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
#ifdef WITH_CLOVER
	Cphi_eopre(static_mass, out, in);
#else
	Dphi_eopre(static_mass, out, in);
#endif
#else
#ifdef WITH_CLOVER
	Cphi(static_mass, out, in);
#else
	Dphi(static_mass, out, in);
#endif
#endif
}


void D_flt(spinor_field_flt *out, spinor_field_flt *in){
#ifdef UPDATE_EO
  Dphi_eopre_flt((float)(static_mass), out, in);
#else
  Dphi_flt((float)(static_mass), out, in);
#endif
}


void Qtm_p(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
    Qhat_eopre(static_mass, static_mu,out, in);
#else
  complex imu;
  imu.re=0;imu.im=static_mu;
  g5Dphi(static_mass, out, in);
  spinor_field_mulc_add_assign_f(out,imu,in);
#endif
}

void Qtm_m(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
    Qhat_eopre(static_mass, -static_mu,out, in);
#else
  complex imu;
  imu.re=0;imu.im=-static_mu;
  g5Dphi(static_mass, out, in);
  spinor_field_mulc_add_assign_f(out,imu,in);
#endif
}

void QpQm_tm(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
    Qhat_eopre_sq(static_mass,static_mu,out,in);
#else
  g5Dphi_sq(static_mass, out, in);
  spinor_field_mul_add_assign_f(out,static_mu*static_mu,in);
#endif
}


void Qtm_p_alt(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
	complex imu;
	imu.re=0;imu.im=static_mu;
#ifdef WITH_CLOVER
	g5Cphi_eopre(static_mass,out,in);
#else
	g5Dphi_eopre(static_mass,out,in);
#endif
	spinor_field_mulc_add_assign_f(out,imu,in);
#else
	complex imu;
	imu.re=0;imu.im=static_mu;
#ifdef WITH_CLOVER
	g5Cphi(static_mass,out,in);
#else
	g5Dphi(static_mass, out, in);
#endif
	spinor_field_mulc_add_assign_f(out,imu,in);
#endif
}

void Qtm_m_alt(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
	complex imu;
	imu.re=0;imu.im=-static_mu;
#ifdef WITH_CLOVER
	g5Cphi_eopre(static_mass,out,in);
#else
	g5Dphi_eopre(static_mass,out, in);
#endif
	spinor_field_mulc_add_assign_f(out,imu,in);
#else
	complex imu;
	imu.re=0;imu.im=-static_mu;
#ifdef WITH_CLOVER
	g5Cphi(static_mass,out,in);
#else
	g5Dphi(static_mass, out, in);
#endif
	spinor_field_mulc_add_assign_f(out,imu,in);
#endif
}

void QpQm_tm_alt(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
#ifdef WITH_CLOVER
	g5Cphi_eopre_sq(static_mass,out,in);
#else
	g5Dphi_eopre_sq(static_mass, out, in);
#endif
	spinor_field_mul_add_assign_f(out,static_mu*static_mu,in);
#else
#ifdef WITH_CLOVER
	g5Cphi_sq(static_mass,out,in);
#else
	g5Dphi_sq(static_mass, out, in);
#endif
	spinor_field_mul_add_assign_f(out,static_mu*static_mu,in);
#endif
}

//This inverts Qtm_p
void tm_invert(spinor_field* out, spinor_field *in, mshift_par* mpar){
  static spinor_field* tmp=NULL;
  if (tmp==NULL){
#ifndef UPDATE_EO
    tmp =  alloc_spinor_field_f(1,&glattice);
#else
    tmp =  alloc_spinor_field_f(1,&glat_even);
#endif
  }
  lprintf("tm_invert",50,"mu: %g, m: %g\n",static_mu,static_mass);
  spinor_field_zero_f(tmp);
  cg_mshift(mpar,QpQm_tm,in,tmp);
  Qtm_m(out,tmp);
}

//This inverts Qtm_p_alt
void tm_invert_alt(spinor_field* out, spinor_field *in, mshift_par* mpar){
  static spinor_field* tmp=NULL;
  if (tmp==NULL){
#ifndef UPDATE_EO
    tmp =  alloc_spinor_field_f(1,&glattice);
#else
    tmp =  alloc_spinor_field_f(1,&glat_even);
#endif
  }
  lprintf("tm_invert_alt",50,"mu: %g, m: %g\n",static_mu,static_mass);
  spinor_field_zero_f(tmp);
  cg_mshift(mpar,QpQm_tm_alt,in,tmp);
  Qtm_m_alt(out,tmp);
}



