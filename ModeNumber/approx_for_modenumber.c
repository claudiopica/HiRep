#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include <stdio.h>


#define NMAX 1000



/* Chebyshev polynomials */
/* T[i,z] = 2zT[i-1,z] - T[i-2,z] */
static void Chebyshev(double *ret, int n, double z) {
  ret[0] = 1.;
  ret[1] = z;
  for(int i=2; i<=n; i++)
    ret[i] = 2.*z*ret[i-1]-ret[i-2];
}

/* Derivative of Chebyshev polynomials */
/* T'[i,z] = 2T[i-1,z] + 2zT'[i-1,z] - T'[i-2,z] */
static void DerChebyshev(double *der, double *T, int n, double z) {
  Chebyshev(T,n,z);
  der[0] = 0.;
  der[1] = 1.;
  for(int i=2; i<=n; i++)
    der[i] = 2.*T[i-1] + 2.*z*der[i-1] - der[i-2];
}

/* Approximating polinomial */
static double P(int n, double *c, double y, double epsilon) {
  double z=(2.*y-1.-epsilon)/(1.-epsilon);
  double T[n+1];
  double ret = 0.;
  
  Chebyshev(T,n,z);
  for(int i=0; i<=n; i++)
    ret += c[i]*T[i];
  
  return ret;
}


/* Derivative of the approximating polinomial */
static double derP(int n, double *c, double y, double epsilon) {
  double z=(2.*y-1.-epsilon)/(1.-epsilon);
  double T[n+1], derT[n+1];
  double ret = 0.;
  
  DerChebyshev(derT,T,n,z);
  for(int i=0; i<=n; i++)
    ret += c[i]*derT[i];
  
  return ret*2./(1.-epsilon);
}


/* Error function */
static double errh(int n, double *c, double y, double epsilon) {
  return 1.-sqrt(y)*P(n,c,y,epsilon);
}


/* Derivative of the error function */
static double dererrh(int n, double *c, double y, double epsilon) {
  return -.5/sqrt(y)*P(n,c,y,epsilon)-sqrt(y)*derP(n,c,y,epsilon);
}


/* 
Solve [ sqrt(y[m])*P(c,y[m]) + (-1)^m u == 1. , {c,u} ]
*/
static void solvesystem(double *c, double *u, int n, double *y, double epsilon) {
  gsl_matrix *A = gsl_matrix_alloc(n+2,n+2);
  gsl_vector *b = gsl_vector_alloc(n+2);
  gsl_vector *x = gsl_vector_alloc(n+2);
  double T[n+1];
  double z, sy;
  double sign=1.;
  
  for(int m=0; m<=n+1; m++) {
    z = (2.*y[m]-1.-epsilon)/(1.-epsilon);
    sy = sqrt(y[m]);
    Chebyshev(T,n,z);
    for(int i=0; i<=n; i++)
      gsl_matrix_set(A,m,i,sy*T[i]);
    gsl_matrix_set(A,m,n+1,sign);
    sign = -sign;
    gsl_vector_set(b,m,1.);
  }

  gsl_linalg_HH_solve(A,b,x);
  for(int i=0; i<=n; i++) {
    c[i] = gsl_vector_get(x,i);
  }
  *u = gsl_vector_get(x,n+1);

  for(int m=0; m<=n+1; m++) {
    z = (2.*y[m]-1.-epsilon)/(1.-epsilon);
    sy = sqrt(y[m]);
    Chebyshev(T,n,z);
    double X=1.;
    for(int i=0; i<=n; i++)
      X = X-c[i]*sy*T[i];
  }

}


/* Zero of the (derivative of the) error function, assuming that exactly one zero exists */
static double zero(double (*fptr)(int,double*,double,double), int n, double *c, double epsilon, double a, double b, double prec) {
  double h1, h2, hmid, mid;
  double prec2 = prec*prec;
  h1 = fptr(n,c,a,epsilon);
  h2 = fptr(n,c,b,epsilon);
  if(h1*h1<prec2) return a;
  if(h2*h2<prec2) return b;
  int counter=0;
  while(1) {
    mid = (a+b)/2.;
    counter++;
    hmid = fptr(n,c,mid,epsilon);
    if(hmid*hmid<prec2) return mid;
    if(h1*hmid>0.) {
      a=mid;
      h1=hmid;
    } else {
      b=mid;
      h2=hmid;
    }
  }
  
  return 0.;
}



static void initialguess(double *y, int n, double epsilon) {
  for(int m=0; m<=n+1; m++) y[m] = (-(1.-epsilon)*cos(m*M_PI/(n+1.))+1.+epsilon)/2.;
}




static int squareroot(double delta, double epsilon, double *c, int *order, double *err) {
  int n = 4;

  double y[NMAX+2];
  double zeroes[NMAX+1];
  double u = 0.;
  double error = 1.;
  double tmp;

  int counter;
  
  while(error > delta) {
    
    n++;
    error = 1.;
    u = 0.;
    counter=0;

    initialguess(y,n,epsilon);

    while(error-fabs(u) > fabs(u)*.01) {
      solvesystem(c,&u,n,y,epsilon);

      for(int i=0; i<=n; i++) {
        zeroes[i] = zero(&errh,n,c,epsilon,y[i],y[i+1],u/1000.);
      }

      error = fabs(errh(n,c,y[0],epsilon));
      for(int i=0; i<n; i++) {
        y[i+1] = zero(&dererrh,n,c,epsilon,zeroes[i],zeroes[i+1],u/(zeroes[i+1]-zeroes[i])/1000.);
        tmp = fabs(errh(n,c,y[i+1],epsilon));
        if(tmp>error) error=tmp;
      }
      tmp = fabs(errh(n,c,y[n+1],epsilon));
      if(tmp>error) error=tmp;
      counter++;
    }

  }

  *order = n;
  *err = error;
  
  return counter;  
}



static double h(double x, double epsilon, int order, double *c) {
  double b0,b1,b2;
  double z=(2.*x*x-1.-epsilon)/(1.-epsilon);
  
  b0=b1=0.;
  for(int n=order; n>=0; n--) {
    b2=b1;
    b1=b0;
    b0 = c[n] + 2.*z*b1 - b2;
  }
  return .5 - .5*x*(b0 - b1*z);
}

typedef struct {
  double epsilon;
  int order;
  double *c;
} error2_pars_type;
static double error2(double x, void *pars) {
  double hv = h(x,((error2_pars_type*)pars)->epsilon,((error2_pars_type*)pars)->order,((error2_pars_type*)pars)->c);
  return hv*hv*hv*hv/((1.-x)*sqrt(1.-x*x));
}



double star(double delta, double epsilon, int order, double *c) {
  double intres, interr;
  double seps = sqrt(epsilon);
  double ret;
  
  error2_pars_type error2_pars;
  error2_pars.epsilon = epsilon; error2_pars.order = order; error2_pars.c = c;
  gsl_function gsl_error2={&error2,&error2_pars};
  
  gsl_integration_workspace *gsl_ws_int = gsl_integration_workspace_alloc(1000000);
  gsl_integration_qag(&gsl_error2, -seps, seps, 1.e-2, 0, 1000000, GSL_INTEG_GAUSS15, gsl_ws_int, &intres, &interr);
  ret = sqrt((1.-seps)/(1.+seps)) + intres;
  ret = ret*ret;
  gsl_integration_workspace_free(gsl_ws_int);

  return ret;
}



int main(int argc, char* argv[]) {
  double epsilon;
  double delta;
  int order;
  double c[NMAX+1];
  double err;

  if(argc != 3) {
    printf("Usage: %s epsilon delta\n",argv[0]);
    return -1;
  }
  
  epsilon = atof(argv[1]);
  delta = atof(argv[2]);

  squareroot(delta, epsilon, c, &order, &err);

/*
  for(int i=0;i<5001;i++) {
    double x = 2.*i/5000.-1.;
    printf("%e %e SIGN\n",x,h(x));
  }
*/
  
  printf("%d\n", order);
  printf("%.16e\n", epsilon);
  printf("%.16e\n\n", err);
  printf("%.16e\n\n", star(delta, epsilon, order, c));
  for(int i=0; i<=order; i++)
    printf("%.16e\n", c[i]);
  
  return 0;
}




