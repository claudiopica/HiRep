#ifndef RATIONAL_FUNCTIONS_H
#define RATIONAL_FUNCTIONS_H

#include "inverters.h"
#include "suN.h"

typedef struct _rational_app {
  unsigned int order;
  double error;
  double *a;
  double *b;
} rational_app;

/* converts the coef from the root/poles form to the partial fraction exp
 * before : r(x)=a[0]*(x-a[1])/(x-b[0])*...*(x-a[n+1])/(x-b[n])
 * after  : r(x)=a[0]+a[1]/(x-b[0])+a[2]/(x-b[1])+...+a[n+1)/(x-b[n])
 */
void r_app_rp2pfe(rational_app *app);

/* converts the coef from the root/poles form to the partial fraction exp
 * as r_app_rp3pfe but for the inverse function
 */
void r_app_rp2pfe_inv(rational_app *app);

/* approximation for f(x)=x^(-1/2) */
/* the order to use is written into app */
void inv_sqrt_rp(rational_app *app, double min, double max);

/* approximation for f(x)=x^(-1/4) */
/* the order to use is written into app */
void inv_fourrt_rp(rational_app *app, double min, double max);

void inv_sqrt_coef(rational_app *app, double min, double max);
void sqrt_coef(rational_app *app, double min, double max);
void inv_fourrt_coef(rational_app *app, double min, double max);
void fourrt_coef(rational_app *app, double min, double max);

void rational_func(rational_app *coef, spinor_operator Q, suNf_spinor *out, suNf_spinor *in);


#endif
