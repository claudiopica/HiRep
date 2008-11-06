/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define APP_ORD 2*25+1

typedef double real;

#define MAXREC 50
static real arithgeom3(real z, real a, real b) {
  real diff[MAXREC], sum[MAXREC], sn;
  real pb = -1.; int n=MAXREC;
  while (n && b>pb) {
    --n; 
    sum[n]=a+b; diff[n]=a-b;
    pb=b; b=sqrt(a*b); a=0.5*(a+pb); 
  }
  if(n==0 && b>pb) {
    printf("ERRORE: raggiunto numero massimo di iterazioni in arithgeom3!!\n");
    return 0.;
  }
  sn=sin(z*a);
  for(;n<MAXREC;++n)
    sn=(sum[n]+diff[n])*sn/(sum[n]+diff[n]*sn*sn);
  return sn;
}

static real sn(real z, real kp) {
  real diff[MAXREC], sum[MAXREC], sn;
  real pb = -1.; real a=1.;
  int n=MAXREC; 
  while (n && kp>pb) {
    --n; 
    sum[n]=a+kp; diff[n]=a-kp;
    pb=kp; kp=sqrt(a*kp); a=0.5*(a+pb); 
  }
  if(n==0 && kp>pb) {
    printf("ERRORE: raggiunto numero massimo di iterazioni in sn!!\n");
    return 0.;
  }
  sn=sin(z*a);
  for(;n<MAXREC;++n)
    sn=(sum[n]+diff[n])*sn/(sum[n]+diff[n]*sn*sn);
  return sn;
}


static real arithgeom_old(real z, real a, real b) {
  static real pb = -1.;
  real xi;

  if(b<=pb) {pb=-1.; return sin(z*a);}
  pb = b;
  xi = arithgeom_old( z, 0.5*(a+b), sqrt(a*b));
  return 2.*a*xi/((a+b)+(a-b)*xi*xi);
}

static real sn_old(real z, real kp) {
  return arithgeom_old(z, 1., kp);
}


/*
static real K_rec(real kp) {
  real agm;
  arithgeom(.5, 1., kp, &agm);
  return  M_PI_2 / agm;
}
*/

static real K(real kp) {
  real pb = -1.; real a=1.;
  while (kp>pb) {
    pb=kp; kp=sqrt(a*kp); a=0.5*(a+pb); 
  }
  return M_PI_2/a;
}

static real dn(real z, real kp) {
  real agm;
  agm = sn(z, kp);
  return sqrt(1.-(1.-kp*kp)*agm*agm);  
}

/*
static real arithgeom2( real z, real a, real b, real *sn) {
  static real pb = -1.;
  real agm;
  
  if(b<=pb) { pb=-1; *sn = sin (z*a); return a;}
  pb = b;
  agm = arithgeom( z, 0.5*(a+b), sqrt(a*b), sn);
  *sn = 2.*a*(*sn)/((a+b)+(a-b)*(*sn)*(*sn));
  return agm;
}

#define MAXREC 100
static real arithgeom3( real z, real a, real b, real *sn) {
  real diff[MAXREC], sum[MAXREC];
  real pb = -1.; int n=MAXREC;
  while (n && b>pb) {
    --n; 
    sum[n]=a+b; diff[n]=a;
    pb=b; b=sqrt(a*b); a=0.5*(a+b); 
  }
  if(n==0) {
    printf("ERRORE: raggiunto numero massimo di iterazioni in arithgeom3!!\n");
    return 0.;
  }
  *sn=sin(z*a);
  for(;n<MAXREC;++n)
    *sn=(sum[n]+diff[n])*(*sn)/(sum[n]+diff[n]*(*sn)*(*sn));
  return a;
}
*/

void zolotarev_coef(int n, real k, real *a, real *b, real *delta){
  real k2, kp, Kp_n, Kp_n2, inv_lambda, inv_xi_bar;
  real z, sn_m, cn_m, inv_c_m, inv_cp_m;

  a[0] = 1.;
  k2=k*k;
  Kp_n = K(k)/((real) n);
  Kp_n2 = 2.*Kp_n; /* Kp over n mult by 2 */

  inv_xi_bar = dn(Kp_n, k);
  inv_lambda = 1./inv_xi_bar;
  inv_xi_bar *= inv_xi_bar; /* inv_xi_bar is squared */

  n >>= 1;
  for (;n>0; --n) {
    z = Kp_n2 * ((real) n);
    sn_m=sn(z, k);
    sn_m *= sn_m; cn_m = 1. - sn_m;
    inv_c_m = -sn_m / cn_m;
    a[n] = k2*inv_c_m;

    z -= Kp_n;
    sn_m=sn(z, k);
    sn_m *= sn_m; cn_m = 1. - sn_m;
    inv_cp_m = -sn_m / cn_m;
    b[n-1] = k2*inv_cp_m;

    inv_lambda *= (inv_xi_bar*inv_c_m-1.)/(inv_xi_bar*inv_cp_m-1.);
    
    a[0] *= (inv_cp_m-1.)/(inv_c_m-1.);
  }

  inv_lambda *= a[0];
  a[0] *= 2./ ((1.+inv_lambda)*k);
  *delta = (inv_lambda-1.) / (1.+inv_lambda);
  
}

real zolotarev_invsqrt(real x, int n, real *a, real *b){
  real res=a[0];
  n >>= 1;
  for (;n>0; --n) {
    /*    printf("%d\n",n); */
    res *= (x-a[n])/(x-b[n-1]);
  }
  return res;
}

real zolotarev_sign(real x, int n, real *a, real *b){
  real res=x*a[0];
  x*=x;
  n >>= 1;
  for (;n>0; --n) {
    /*    printf("%d\n",n); */
    res *= (x-a[n])/(x-b[n-1]);
  }
  return res;
}

void norm_coef(real min, real max, int n, real *a, real *b, real *delta) {
  real k=sqrt(min/max);
  zolotarev_coef(n, k, a, b, delta);
  a[0] /= sqrt(max);
  n >>= 1;
  for (;n>0; --n) {
    a[n] *= max;
    b[n-1] *= max;
  }
}

void partial_frac_coef(int n, real *a, real*b) {
  int order = n>>1;
  int i, k;
  real ap[order];
  
  for (i=0; i<order; ++i)
    ap[i] = a[i+1];

  for (i=0; i<order; ++i) {
    a[i+1] = a[0];
    for (k=0; k<order; ++k) {
      a[i+1] *= (b[i]-ap[k]);
      if (k!=i)
	a[i+1] /= (b[i]-b[k]);
    }
  }

}

real pfrac_eval(real x, int n, real *a, real *b) {
  real res = a[0];
  n >>= 1;
  for (;n>0; --n) {
    res += a[n]/(x-b[n-1]);
  }
  return res;
  
}

int main()
{
  real a[(APP_ORD>>1)+1], b[(APP_ORD>>1)];
  real delta, fn, x, app, app2;
  /* real k=.003; */
  real min=0.4;
  real min2, delta2, err, err2;
  real max=64.0;
  int np, i, j;
  real threshold = 1.e-15;
	real diff; /* for checks */

  /*
  k=sqrt(k);
  zolotarev_coef(APP_ORD, k, a, b, &delta);
  */

  /*  
  for (min/=1.5;min>1e-10;min/=1.2) {
    printf("[%g,%g]\t", min, max);
    for (np=3; np<APP_ORD+1; np+=2) {
    norm_coef(min, max, np, a, b, &delta);
    if (delta<threshold) {
      printf("%d\t%g\n", np>>1, delta);
      break;
    }
  }
  }
  */

  min = 60.;
  min2 = 1.e-10;
  delta2 = min - min2;

  np = APP_ORD;
  norm_coef(min, max, np, a, b, &err);
  norm_coef(min2, max, np, a, b, &err2);
  if (err>threshold || err2<threshold) {
    printf("ERRORE: limiti sbagliati!!!\n");
    exit(1);
  }

  do {
    /* printf("[%g,%g]\t", min, max); */
    real err3;
    real mean = (min + min2)/2.;
    norm_coef(mean, max, np, a, b, &err3);

    if (err3>threshold) {
      min2 = mean;
      err2 = err3;
    } else {
      min = mean;
      err = err3;
    }
    delta2 = min - min2; 
    printf("Min = %g ; Min2 = %g ; Delta = %g\n", min, min2, delta2);

  } while (delta2>2.e-15);
  
  printf("[%1.16e,%g]\t", min, max);
  printf("%d\t%g\n", np>>1, err);
  printf("[%1.16e,%g]\t", min2, max);
  printf("%d\t%1.16e\n", np>>1, err2);

  norm_coef(min, max, np, a, b, &delta);

  partial_frac_coef(np, a, b);
 
  printf("Ra_a0 = %1.16e\n",a[0]);
  for(i=0;i<(np>>1);++i) {
    /*    printf("[%d] ", (np>>1)-i-1);*/
    printf("Ra_a[%d] = %1.16e ; ", (np>>1)-i-1,a[i+1]);
    printf("Ra_b[%d] = %1.16e ;\n", (np>>1)-i-1,b[i]);
  }

   
  printf("[%1.16e,%g]\t", min, max);
  delta = max-min;
  for (np=0; np<10000; ++np) {
    x = min + ((real) np) / ((real)9999) * delta;
    fn = 1./ sqrt(x);
    /*  
    app = zolotarev_invsqrt(x, APP_ORD, a, b);
    app2 = zolotarev_sign(x, APP_ORD, a, b);
    */
    diff = pfrac_eval(x, APP_ORD, a, b);
    printf("%g\t%g\t%g\n", x, delta, (fn-diff)/fabs(fn));
  }
  
  return 0;
}
