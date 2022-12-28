/*******************************************************************************
 *
 * File  diag_hmat.c
 *
 * Function to diagonalize NxN hermitean matrix 
 *
 * Ari Hietanen
 *
 *******************************************************************************/

#ifdef GAUGE_SON

#include "utils.h"
#include "libhr_core.h"
#include "io.h"

void diag_hmat(suNg *hmat, double *diag){
  double offdiag[NG];
  tridiagonalize(hmat,diag,offdiag);
  diag_tridiag(hmat,diag,offdiag);
}

void tridiagonalize(suNg *hmat, double *diag, double* roffdiag){
  suNg_vector offdiag,realtrans;
  hr_complex f,g,ctmp;
  double h,scale,fa,hh;
  int i,j,k,l;
  for (i=NG-1;i>0;i--){//loop over the components of column vector 
    l = i-1;
    h = scale = 0.0;
    if (l>0){ //Not the last iteration
      for (k=0;k<l+1;k++) //i'th row, called x in comments
	scale += fabs(hmat->c[i*NG+k]);
      if (scale == 0.0) ///Skip transformation
        {
          offdiag.c[i] = hmat->c[i*NG+l];
        }
      else{
	for (k=0;k<l+1;k++){
	  hmat->c[i*NG+k]*=1.0/scale; //Use scaled matrix for transformation
	  h += hmat->c[i*NG+k]*hmat->c[i*NG+k]; //_complex_prod_re(hmat->c[i*NG+k],hmat->c[i*NG+k]); //calculate normalization h=|x|^2
	}
        f = hmat->c[i*NG+l]; //f=x[0]
        fa = sqrt(creal(f)*creal(f));
	if (fa!=0){ // g = -x[0]*|x|/|x[0]|
	  _complex_mulr(g,-sqrt(h)/fa,f);
	}
	else{
	  g = -sqrt(h);
	}
	_complex_mulr(offdiag.c[i],scale,g); //stores off-diagonal values
	h -= _complex_prod_re(f,g); // h = 1/2 |u|^2; u = x - e[l]*|x|
	hmat->c[i*NG+l]=creal(f)-cimag(g);
	_complex_0(f);
	for (j=0;j<l+1;j++){
	  hmat->c[j*NG+i] = hmat->c[i*NG+j]/h; //Keep matrix hermitean
	  _complex_0(g);
	  for (k=0;k<j+1;k++){
	    g+=hmat->c[i*NG+k]*hmat->c[j*NG+k];
	  }
	  for (k=j+1;k<l+1;k++){
	    g+=hmat->c[i*NG+k]*hmat->c[k*NG+j];
	  }
	  _complex_mulr(offdiag.c[j],1.0/h,g);
	  _complex_mulr_assign(f,hmat->c[i*NG+j],offdiag.c[j]);
	}

	hh = creal(f)/(h+h);
	for (j=0;j<l+1;j++){
	  f = hmat->c[i*NG+j]+I*cimag(f);
	  _complex_mulr(ctmp,hh,f);
	  _complex_sub(g,offdiag.c[j],ctmp);
	  offdiag.c[j] = g;
	  for (k=0;k<j+1;k++){
	    _complex_mulr(ctmp,hmat->c[i*NG+k],g);
            _complex_star(ctmp,ctmp);
	    _complex_mul_star_assign(ctmp,offdiag.c[k],f);
	    hmat->c[j*NG+k]-=creal(ctmp);
	  }
	}
      }
    }
    else{
      offdiag.c[i]=hmat->c[i*NG+l]+I*cimag(offdiag.c[i]);
    }
    diag[i]=h;
  }
  diag[0] = 0.0;
  _complex_0(offdiag.c[0]);
  for (i=0;i<NG;i++){
    l=i;
    if (diag[i] != 0.0){
      for (j=0;j<l;j++){
	_complex_0(g);
	for (k=0;k<l;k++){
	  g+=hmat->c[NG*i+k]*hmat->c[NG*k+j];
	}
	for (k=0;k<l;k++){
	  _complex_mulr(ctmp,hmat->c[NG*k+i],g);
	  hmat->c[NG*k+j]-=creal(ctmp);	  
	}
      }
    }
    diag[i] = hmat->c[i*NG+i];
    hmat->c[i*NG+i]=1.;
    for (j=0;j<i;j++){
      hmat->c[j*NG+i]=0.;
      hmat->c[i*NG+j]=0.;
    }
  }
  _complex_1(realtrans.c[0]);
  roffdiag[0]=0.0;
  for (i=1;i<NG;i++){
    roffdiag[i] = sqrt(_complex_prod_re(offdiag.c[i],offdiag.c[i]));
    if (roffdiag[i] == 0.0){
      _complex_1(realtrans.c[i]);
    }
    else{
      _complex_div(ctmp,realtrans.c[i-1],offdiag.c[i]);
      _complex_mulr(realtrans.c[i],roffdiag[i],ctmp);
    }
  }
  for (i=0;i<NG;i++){
    for (j=0;j<NG;j++){
      _complex_mulr(ctmp,hmat->c[NG*j+i],realtrans.c[i]);
      _complex_star(ctmp,ctmp);
      hmat->c[NG*j+i]=creal(ctmp);
    }
  }
  for (i=1;i<NG;i++) roffdiag[i-1]=roffdiag[i];
  roffdiag[NG-1] = 0.0;
}

void diag_tridiag(suNg* hmat, double *diag, double* offdiag){
  int i,k,l,m,iter;
  double g,r,s,c,p,f,b,dd;
  double  ctmp;
  for (l=0;l<NG;l++){
    iter = 0;
    do {
      for (m=l;m<NG-1;m++){
	dd = fabs(diag[m]) + fabs(diag[m+1]);
	if (fabs(offdiag[m]) + dd == dd) break;
      }
      if (m!=l){
	if (++iter == 31) {return;} //Too many iterations
	g = (diag[l+1]-diag[l])/(2.0*offdiag[l]);
	r = pythag(g,1.0);
	g = diag[m] - diag[l] + offdiag[l]/(g + sign(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--){
	  f = s*offdiag[i];
	  b = c*offdiag[i];
	  r = pythag(f,g);
	  offdiag[i+1]=r;
	  if (r==0.0){
	    diag[i+1] -= p;
	    offdiag[m] = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = diag[i+1]-p;
	  r= (diag[i]-g)*s + 2.0*c*b;
	  p=s*r;
	  diag[i+1] = g + p;
	  g = c*r - b;
	  for (k=0;k<NG;k++){
            ctmp = hmat->c[k*NG+i+1];
	    hmat->c[k*NG+i+1]=s*hmat->c[k*NG+i];
	    hmat->c[k*NG+i+1]+=c*ctmp;
	    hmat->c[k*NG+i]=c*hmat->c[k*NG+i];
	    hmat->c[k*NG+i]-=s*ctmp;
	  }
	}
	if ( r == 0 && i>=l) continue;
	diag[l] -= p;
	offdiag[l] = g;
	offdiag[m] = 0.0;
      }
    } while (m != l);
  }
}

double pythag(double a, double b){
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa>absb) return absa*sqrt(1.0+absb*absb/absa/absa);
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+absa*absa/absb/absb));
}

double sign(double a, double b){
  return (b >= 0) ? ( a >= 0 ? a : -a) : ( a>=0 ? -a : a);
}

#endif
