#include <iostream>
#include <cmath>
#include "effective_mass.h"

using namespace std;

double fit(int left, int right, Corr_t *data) {
  double res=0., w=0.;
  for (int i=left;i<=right;++i) {
    double val=data->d[i].avr().val;
    double var=data->d[i].stderr().val;
    if (var==0. || val==0.) {cerr<<"Zero error: "<<val<<" "<<var<<endl;}
    var*=var; var=1./var;
    res+=val*var;
    w+=var;
  }
  res/=w;
  return res;
}


/* Effective mass defined as in JHEP 0702:082,2007 [hep-lat/0701009] */

double effm_err() {
  //cerr<<"Error in the evaluation of the effective mass!!\n";
  return -1.;
}

/* this function solve the equation
 * cosh[ (a+1) M ] = k cosh[ a M ]
 */
double plain_eff_mass(int a, double k) {
  double M,f,d,oldM;
  //int i=0;

  /* check for impossible equations */
  if (k<1. && a>=0) return effm_err();
  if (k>1. && a<0) return effm_err();
  if (a==0) return acosh(k);
  if (a==-1) return acosh(1./k);

  if (a>0) M=log(k);
  else M=log(1./k);

  if(fabs(((double)a)*M)>26) return M;

  oldM=M+1.;

  f=cosh(((double)(a+1))*M)-k*cosh(((double)a)*M);
  while (fabs(f)>1.e-14 && fabs(oldM-M)>1.e-14) {
    d=((double)(a+1))*sinh(((double)(a+1))*M)-k*((double)a)*sinh(((double)a)*M);
    oldM=M;
    M-=f/d;
    f=cosh(((double)(a+1))*M)-k*cosh(((double)a)*M);
    //cout<<"M="<<M<<" f="<<fabs(f)<<" diff="<<fabs(oldM-M)<<endl;
    //++i;
  }

  //cout<<"iters="<<i<<endl;

  return M;

}

double hc(int t, double m, int Lt) {
  return (exp(-m*double(Lt-t))+exp(-m*double(t)));
}



/* Effective masses using Prony's method
 * one and two effectie masses only
 */


double bin(int n, int k) {
  int i, b, c;
  c=n-k;
  if (k>c) {b=k; k=c; c=b;} /* make k<= n/2 */
  if (k==0) return 1.;
  b=(++c);
  for (i=2;i<=k;++i) {
    b*=(++c);
    b/=i;
  }
  return (double)(b);
}

void build_y(double *c, int n, double *dst) {
  double y[n];
  for (int i=0;i<=n/2;++i){
    int j=n/2-i;
    double z=2.;
    y[j]=0.;
    for (int k=0;k<=i;++k){
      y[j]+=bin(i,k)*c[n/2+i-2*k];
      z*=.5;
    }
    y[j]*=z;
  }
  dst[0]=y[0];
  for (int i=1;i<=n/2;++i) {
    dst[i]=dst[n-i]=y[i];
  }
}

int centered_prony_eff_mass_1(double *y, double *m1) {
  double y1,y2;
  double x1;
  y1=y[1]; y2=y[0];
  
  x1=y2/y1;

  if(x1<1.) {
    //cerr<<"x1<1 !\n";
    *m1=-1.;
    return 0;
  } else {
    *m1=acosh(x1);
  }

return 1;
}

int centered_prony_eff_mass_2(double *y, double *m1, double *m2) {
  double y1,y2,y3,y4;
  double a,b,c,d,x1,x2;
  y1=y[3]; y2=y[2]; y3=y[1]; y4=y[0];
  
  a=y1*y3-y2*y2;
  b=y2*y3-y1*y4;
  c=y2*y4-y3*y3;

  d=b*b-4.*a*c;

  if(d<0.) {
    //cerr<<"Delta <0 !\n";
    *m1=*m2=-1.;
    return 0;
  }

  x1=sqrt(d);
  x2=(-b-x1)/(2.*a);
  x1-=b; x1/=(2.*a);

  if(x2<x1) {//swap x1 <=> x2
    d=x2;
    x2=x1;
    x1=d;
  }
  
  if(x1<1.) {
    //cerr<<"x1<1 !\n";
    *m1=-1.;
  } else {
    *m1=acosh(x1);
  }

  if(x2<1.) {
    //cerr<<"x2<1 !\n";
    *m2=-1.;
  } else {
    *m2=acosh(x2);
  }

  if (*m1<0. && *m2<0.) { return 0; }
  if (*m1<0. && *m2>0.) {
    *m1=*m2; *m2=-1;
  } 

  return 1;

}





int shifted_prony_eff_mass_1(double* C, int t, int tmax, double* m1, int lt) {
	double y[2], x1;
	int q, tau;
	
	q = tmax-t;
	tau = tmax;
	
  for (int i=0;i<2;++i){
    double z=2.;
    y[i]=0.;
    for (int k=0;k<=q+i;++k){
      int s=(t-i+2*k+lt)%lt;
      y[i]+=bin(q+i,k)*C[(s<=lt/2)?s:(lt-s)];
      z*=.5;
    }
    y[i]*=z;
  }

    
  x1=y[1]/y[0];

  if(x1<1.) {
    //cerr<<"x1<1 !\n";
    *m1=-1.;
    return 0;
  } else {
    *m1=acosh(x1);
  }

  return 1;

}


int shifted_prony_eff_mass_2(double* C, int t, int tmax, double* m1, double *m2, int lt) {
	double y[4];
	double a,b,c,d;
	double x1,x2;
	int q, tau;
	
	q = tmax-t;
	tau = tmax;
	
  for (int i=0;i<4;++i){
    double z=2.;
    y[i]=0.;
    for (int k=0;k<=q+i;++k){
      int s=(t-i+2*k+lt)%lt;
      y[i]+=bin(q+i,k)*C[(s<=lt/2)?s:(lt-s)];
      z*=.5;
    }
    y[i]*=z;
  }
  
  a=y[0]*y[2]-y[1]*y[1];
  b=y[1]*y[2]-y[0]*y[3];
  c=y[1]*y[3]-y[2]*y[2];
  
  d=b*b-4.*a*c;

  if(d<0.) {
    //cerr<<"Delta <0 !\n";
    *m1=*m2=-1.;
    return 0;
  }

  x1=sqrt(d);
  x2=(-b-x1)/(2.*a);
  x1-=b; x1/=(2.*a);

  if(x2<x1) {//swap x1 <=> x2
    d=x2;
    x2=x1;
    x1=d;
  }
  
  if(x1<1.) {
    //cerr<<"x1<1 !\n";
    *m1=-1.;
  } else {
    *m1=acosh(x1);
  }

  if(x2<1.) {
    //cerr<<"x2<1 !\n";
    *m2=-1.;
  } else {
    *m2=acosh(x2);
  }

  if (*m1<0. && *m2<0.) { return 0; }
  if (*m1<0. && *m2>0.) {
    *m1=*m2; *m2=-1;
  } 

  return 1;
}


