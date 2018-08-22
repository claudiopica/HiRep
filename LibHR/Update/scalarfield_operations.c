/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include <string.h>
#include "spinor_field.h"
#include "update.h"
#include "random.h"

void scalar_field_copy(scalar_field *s1, scalar_field *s2)
{ 
  error(s1==NULL,1,"scalar_field_copy [scalarfield_operations.c]",
	"Attempt to access unallocated memory space");
  error(s2==NULL,1,"scalar_field_copy [scalarfield_operations.c]",
	"Attempt to access unallocated memory space");
  
  _MASTER_FOR(&glattice,ix) {
    double c = *_FIELD_AT(s2,ix);
    *_FIELD_AT(s1,ix) = c;
  }
}

//s=-s
void flip_scalar_field(scalar_field *s)
{
  error(s==NULL,1,"flip_scalar_field [scalarfield_operations.c]",
	"Attempt to access unallocated memory space");
  
  _MASTER_FOR(&glattice,ix) {
    double c = *_FIELD_AT(s,ix);
    *_FIELD_AT(s,ix) = -c;
  }
}

//s=c
void set_scalar_field(scalar_field *s, double c)
{
  error(s==NULL,1,"set_scalar_field [scalarfield_operations.c]",
	"Attempt to access unallocated memory space");
  
  _MASTER_FOR(&glattice,ix) {
    *_FIELD_AT(s,ix) = c;
  }
}

void gaussian_scalar_field(scalar_field *s)
{
  error(s==NULL,1,"gaussian_scalar_field [scalarfield_operations.c]",
	"Attempt to access unallocated memory space");   
  
  _MASTER_FOR(&glattice,ix) {
    double *ptr = _FIELD_AT(s,ix);
    gauss(ptr,1);
  }
}





// out += (sigma(x)+rho)*in 
void spinor_scalarfield_mult_add_assign(spinor_field *out,scalar_field *sigma,double rho, spinor_field *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor s,r;
     double k = *_FIELD_AT(sigma,i)+rho;
     s = *_FIELD_AT(in,i);
     r = *_FIELD_AT(out,i);
     _spinor_mul_add_assign_f(r,k,s);
     //if(i==3002) printf("mult sigma: %g\n",k);
     *_FIELD_AT(out,i) = r;
  }
}

//out += i*g5*pi(x)*in
void spinor_scalarfield_ig5_mult_add_assign(spinor_field *out,scalar_field *pi, spinor_field *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor s,r;
     double k = *_FIELD_AT(pi,i);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = k; 
     s = *_FIELD_AT(in,i);
     r = *_FIELD_AT(out,i);
     _spinor_g5_assign_f(s);
     _spinor_mulc_add_assign_f(r,z,s);
     //if(i==3002) printf("mult pi: %g\n",k);
     *_FIELD_AT(out,i) = r;
  }
}
//out += -i*g5*pi(x)*in
void spinor_scalarfield_mig5_mult_add_assign(spinor_field *out,scalar_field *pi, spinor_field *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor s,r;
     double k = *_FIELD_AT(pi,i);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = -k; 
     s = *_FIELD_AT(in,i);
     r = *_FIELD_AT(out,i);
     _spinor_g5_assign_f(s);
     _spinor_mulc_add_assign_f(r,z,s);
     //if(i==3002) printf("mult pi: %g\n",k);
     *_FIELD_AT(out,i) = r;
  }
}
// out(x) = -1/(sigma(x)+rho + i*g5*pi(x)) in(x) = (-sigma(x)+rho + i*g5*pi(x))/((sigma(x)+rho)^2 + pi(x)^2)
void spinor_sigma_pi_rho_minus_div_assign(spinor_field *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor s,g5s,r;
     double tm = *_FIELD_AT(sigma,i)+rho;
     double tp = *_FIELD_AT(pi,i);
     double tsq = 1./(tm*tm+tp*tp);
     s = *_FIELD_AT(in,i);
     _spinor_mul_f(r,-tm*tsq,s);
     _spinor_g5_f(g5s,s);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = tp*tsq;
     _spinor_mulc_add_assign_f(r,z,g5s);
     //if(i==8000) printf("div sigma: %g  %g\n",tm,tp);
     *_FIELD_AT(out,i) = r;
  }
}
// out(x) = 1/(sigma(x)+rho + i*g5*pi(x)) in(x) = (sigma(x)+rho - i*g5*pi(x))/((sigma(x)+rho)^2 + pi(x)^2)
void spinor_sigma_pi_rho_div_assign(spinor_field *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor s,g5s,r;
     double tm = *_FIELD_AT(sigma,i)+rho;
     double tp = *_FIELD_AT(pi,i);
     double tsq = 1./(tm*tm+tp*tp);
     s = *_FIELD_AT(in,i);
     _spinor_mul_f(r,tm*tsq,s);
     _spinor_g5_f(g5s,s);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = -tp*tsq;
     _spinor_mulc_add_assign_f(r,z,g5s);
     //if(i==8000) printf("div sigma: %g  %g\n",tm,tp);
     *_FIELD_AT(out,i) = r;
  }
}
// out(x) = 1/(sigma(x)+rho - i*g5*pi(x)) in(x) = (sigma(x)+rho + i*g5*pi(x))/((sigma(x)+rho)^2 + pi(x)^2)
void spinor_sigma_pi_dagger_rho_div_assign(spinor_field *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor s,g5s,r;
     double tm = *_FIELD_AT(sigma,i)+rho;
     double tp = *_FIELD_AT(pi,i);
     double tsq = 1./(tm*tm+tp*tp);
     s = *_FIELD_AT(in,i);
     _spinor_mul_f(r,tm*tsq,s);
     _spinor_g5_f(g5s,s);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = tp*tsq;
     _spinor_mulc_add_assign_f(r,z,g5s);
     //if(i==8000) printf("div sigma: %g  %g\n",tm,tp);
     *_FIELD_AT(out,i) = r;
  }
}
// out(x) = (sigma(x)+rho + i*g5*pi(x))
void spinor_sigma_pi_dagger_rho_assign(spinor_field *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor s,g5s,r;
     double tm = *_FIELD_AT(sigma,i)+rho;
     double tp = *_FIELD_AT(pi,i);
     s = *_FIELD_AT(in,i);
     _spinor_mul_f(r,tm,s);
     _spinor_g5_f(g5s,s);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = tp;
     _spinor_mulc_add_assign_f(r,z,g5s);
     //if(i==8000) printf("div sigma: %g  %g\n",tm,tp);
     *_FIELD_AT(out,i) = r;
  }
}
// out(x) = sigma(x) in(x)
void spinor_sigma_assign(spinor_field *out,scalar_field *sigma, spinor_field *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor s,r;
     double tm = *_FIELD_AT(sigma,i);
     s = *_FIELD_AT(in,i);
     _spinor_mul_f(r,tm,s);
     *_FIELD_AT(out,i) = r;
  }
}




// in = (sigma(x)+rho)*in 
void spinor_scalarfield_mult_add_assign_flt(spinor_field_flt *out,scalar_field *sigma,double rho, spinor_field_flt *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor_flt s,r;
     double k = *_FIELD_AT(sigma,i)+rho;
     s = *_FIELD_AT(in,i);
     r = *_FIELD_AT(out,i);
     _spinor_mul_add_assign_f(r,k,s);
     //if(i==3002) printf("mult sigma flt: %g\n",k);
     *_FIELD_AT(out,i) = r;
  }
}

//in += i*g5*pi(x)*in
void spinor_scalarfield_ig5_mult_add_assign_flt(spinor_field_flt *out,scalar_field *pi, spinor_field_flt *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor_flt s,r;
     double k = *_FIELD_AT(pi,i);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = k; 
     s = *_FIELD_AT(in,i);
     r = *_FIELD_AT(out,i);
     _spinor_g5_assign_f(s);
     _spinor_mulc_add_assign_f(r,z,s);
     //if(i==3002) printf("mult pi flt: %g\n",k);
     *_FIELD_AT(out,i) = r;
  }
}
//in -= i*g5*pi(x)*in
void spinor_scalarfield_mig5_mult_add_assign_flt(spinor_field_flt *out,scalar_field *pi, spinor_field_flt *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor_flt s,r;
     double k = *_FIELD_AT(pi,i);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = -k; 
     s = *_FIELD_AT(in,i);
     r = *_FIELD_AT(out,i);
     _spinor_g5_assign_f(s);
     _spinor_mulc_add_assign_f(r,z,s);
     //if(i==3002) printf("mult pi flt: %g\n",k);
     *_FIELD_AT(out,i) = r;
  }
}


// out(x) = 1/(sigma(x)+rho + i*g5*pi(x)) in(x) = (sigma(x)+rho - i*g5*pi(x))/((sigma(x)+rho)^2 + pi(x)^2)
void spinor_sigma_pi_rho_div_assign_flt(spinor_field_flt *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field_flt *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor_flt s,g5s,r;
     double tm = *_FIELD_AT(sigma,i)+rho;
     double tp = *_FIELD_AT(pi,i);
     double tsq = 1./(tm*tm+tp*tp);
     s = *_FIELD_AT(in,i);
     _spinor_mul_f(r,tm*tsq,s);
     _spinor_g5_f(g5s,s);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = -tp*tsq;
     _spinor_mulc_add_assign_f(r,z,g5s);
     //if(i==8000) printf("div sigma: %g  %g\n",tm,tp);
     *_FIELD_AT(out,i) = r;
  }
}
// out(x) = 1/(sigma(x)+rho - i*g5*pi(x)) in(x) = (sigma(x)+rho + i*g5*pi(x))/((sigma(x)+rho)^2 + pi(x)^2)
void spinor_sigma_pi_dagger_rho_div_assign_flt(spinor_field_flt *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field_flt *in){
  _MASTER_FOR(in->type,i) {
     suNf_spinor_flt s,g5s,r;
     double tm = *_FIELD_AT(sigma,i)+rho;
     double tp = *_FIELD_AT(pi,i);
     double tsq = 1./(tm*tm+tp*tp);
     s = *_FIELD_AT(in,i);
     _spinor_mul_f(r,tm*tsq,s);
     _spinor_g5_f(g5s,s);
     complex z;
     _complex_re(z) = 0;
     _complex_im(z) = tp*tsq;
     _spinor_mulc_add_assign_f(r,z,g5s);
     //if(i==8000) printf("div sigma: %g  %g\n",tm,tp);
     *_FIELD_AT(out,i) = r;
  }
}




