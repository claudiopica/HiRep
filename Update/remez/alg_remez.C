/*
	Modified by C. Pica - 10 May 2007

	In function search: increase delta if spread is decreasing to 
	                    accelerate convergence
	Using new initial approximation (or a user supplied guess)

#####################################################################

  Mike Clark - 25th May 2005

  alg_remez.C

  AlgRemez is an implementation of the Remez algorithm, which in this
  case is used for generating the optimal nth root rational
  approximation.

  Note this class requires the gnu multiprecision (GNU MP) library.

*/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"alg_remez.h"


// Constructor
AlgRemez::AlgRemez(double lower, double upper, long precision) 
{
  prec = precision;
  bigfloat::setDefaultPrecision(prec);

  apstrt = lower;
  apend = upper;
  apwidt = apend - apstrt;

	/*
  printf("Approximation bounds are [%e,%e]\n", (double)apstrt,(double)apend);
  printf("Precision of arithmetic is %d\n", precision);
	*/

  alloc = 0;
  n = 0;
  d = 0;

  foundRoots = 0;

  // Only require the approximation spread to be less than 1 ulp
  tolerance = 1e-16;

}

// Destructor
AlgRemez::~AlgRemez()
{
  if (alloc) {
    delete [] param;
    delete [] roots;
    delete [] poles;
    delete [] xx;
    delete [] mm;
    delete [] a_power;
    delete [] a;
  }
}

// Free memory and reallocate as necessary
void AlgRemez::allocate(int num_degree, int den_degree)
{
  // Arrays have previously been allocated, deallocate first, then allocate
  if (alloc) {
    delete [] param;
    delete [] roots;
    delete [] poles;
    delete [] xx;
    delete [] mm;
  }

  // Note use of new and delete in memory allocation - cannot run on qcdsp
  param = new bigfloat[num_degree+den_degree+1];
  roots = new bigfloat[num_degree];
  poles = new bigfloat[den_degree];
  xx = new bigfloat[num_degree+den_degree+2];
  mm = new bigfloat[num_degree+den_degree+2];

  if (!alloc) {
    // The coefficients of the sum in the exponential
    a = new bigfloat[SUM_MAX];
    a_power = new int[SUM_MAX];
  }

  alloc = 1;
}

// Reset the bounds of the approximation
void AlgRemez::setBounds(double lower, double upper)
{
  apstrt = lower;
  apend = upper;
  apwidt = apend - apstrt;
}

// Generate the rational approximation x^(pnum/pden)
double AlgRemez::generateApprox(int degree, unsigned long pnum, 
				unsigned long pden, bigfloat *dmm)
{
  return generateApprox(degree, degree, pnum, pden, dmm);
}

double AlgRemez::generateApprox(int num_degree, int den_degree, 
				unsigned long pnum, unsigned long pden, bigfloat *dmm)
{
  double *a_param = 0;
  int *a_pow = 0;
  return generateApprox(num_degree, den_degree, pnum, pden, dmm, 0, a_param, a_pow);
}

// Generate the rational approximation x^(pnum/pden)
double AlgRemez::generateApprox(int num_degree, int den_degree, 
				unsigned long pnum, unsigned long pden,
				bigfloat *dmm, int a_len, double *a_param, int *a_pow)
{
  char *fname = "generateApprox(int, unsigned long, unsigned long)";

	/*
  printf("Degree of the approximation is (%d,%d)\n", num_degree, den_degree);
  printf("Approximating the function x^(%d/%d)\n", pnum, pden);
	*/

  // Reallocate arrays, since degree has changed
  if (num_degree != n || den_degree != d) allocate(num_degree,den_degree);

  if (a_len > SUM_MAX) {
    printf("Error: a_length > SUM_MAX");
    exit(0);
  }

  step = new bigfloat[num_degree+den_degree+2];

  a_length = a_len;
  for (int j=0; j<a_len; j++) {
    a[j]= a_param[j];
    a_power[j] = a_pow[j];
  }

  power_num = pnum;
  power_den = pden;
  spread = 1.0e37;
  iter = 0;

  n = num_degree;
  d = den_degree;
  neq = n + d + 1;

	if(dmm!=0){
		/*
		for (int j=0;j<neq+1;++j)
			mm[j]=apstrt+bigfloat(dmm[j])*(apend-apstrt);
		for (int j=0;j<neq+2;++j)
			xx[j]=apstrt+bigfloat(dmm[neq+1+j])*(apend-apstrt);
		*/
		for (int j=0;j<neq+1;++j)
			mm[j]=exp_bf(log_bf(apstrt)+dmm[j]*(log_bf(apend)-log_bf(apstrt)));
		for (int j=0;j<neq+1;++j)
			xx[j]=exp_bf(log_bf(apstrt)+dmm[neq+1+j]*(log_bf(apend)-log_bf(apstrt)));
	} else {
		initialGuess();
	}
  stpini(step);

  while (spread > tolerance) { //iterate until convergance

		++iter;
		/*
    if (iter%100==0) 
      printf("Iteration %d, spread %e delta %e\n", 
	     iter-1,(double)spread,(double)delta);
		*/

    equations();
    if (delta < tolerance) {
      fprintf(stderr,"Delta too small, try increasing precision\n");
      exit(0);
    }

    search(step);

  }

  int sign;
  double error = (double)getErr(mm[0],&sign);
  fprintf(stderr,"Converged at %d iterations, error = %e\n",iter,error);

  // Once the approximation has been generated, calculate the roots
  if(!root()) {
    fprintf(stderr,"Root finding failed\n");
  } else {
    foundRoots = 1;
  }
  
  delete [] step;

  // Return the maximum error in the approximation
  return error;
}

void AlgRemez::getMM(bigfloat *dmm) {
	/*
	for (int i=0;i<neq+1;++i)
		dmm[i]=(double)((mm[i]-apstrt)/(apend-apstrt));
	for (int i=0;i<neq+2;++i)
		dmm[neq+1+i]=(double)((xx[i]-apstrt)/(apend-apstrt));
	*/
	for (int i=0;i<neq+1;++i)
		dmm[i]=(log_bf(mm[i])-log_bf(apstrt))/(log_bf(apend)-log_bf(apstrt));
	for (int i=0;i<neq+1;++i)
		dmm[neq+1+i]=(log_bf(xx[i])-log_bf(apstrt))/(log_bf(apend)-log_bf(apstrt));
/*
#define SCALE(x) ((x)/(log(apend)-log(apstrt)))
#define SCALEO(x) (((x)-log(apstrt))/(log(apend)-log(apstrt)))

	printf("ll+ %e %e\n",double(apstrt),double(apend));
	printf("mm+0 %e %e\n",SCALEO(log(apstrt)),SCALE(log(mm[0])-log(apstrt)));
	for (int i=0;i<neq;++i)
		printf("mm+ %e %e\n",SCALEO(log(mm[i])),SCALE(log(mm[i+1])-log(mm[i])));
	printf("mm+ %e %e\n",SCALEO(log(mm[neq])),SCALE(log(apend)-log(mm[neq])));
	printf("xx+ %e %e\n",SCALEO(log(apstrt)),SCALE(log(xx[0])-log(apstrt)));
	for (int i=0;i<neq+2;++i)
		printf("xx+ %e %e\n",SCALEO(log(xx[i])),SCALE(log(xx[i+1])-log(xx[i])));
	printf("xx+ %e %e\n",SCALEO(log(xx[neq+1])),SCALE(log(apend)-log(xx[neq+1])));
	*/
}

// Return the partial fraction expansion of the approximation x^(pnum/pden)
int AlgRemez::getPFE(double *Res, double *Pole, double *Norm) {

  if (n!=d) {
    printf("Cannot handle case: Numerator degree neq Denominator degree\n");
    return 0;
  }

  if (!alloc) {
    printf("Approximation not yet generated\n");
    return 0;
  }

  if (!foundRoots) {
    printf("Roots not found, so PFE cannot be taken\n");
    return 0;
  }

  bigfloat *r = new bigfloat[n];
  bigfloat *p = new bigfloat[d];
  
  for (int i=0; i<n; i++) r[i] = roots[i];
  for (int i=0; i<d; i++) p[i] = poles[i];
  
  // Perform a partial fraction expansion
  pfe(r, p, norm);

  // Convert to double and return
  *Norm = (double)norm;
  for (int i=0; i<n; i++) Res[i] = (double)r[i];
  for (int i=0; i<d; i++) Pole[i] = (double)p[i];

  delete [] r;
  delete [] p;

  // Where the smallest shift is located
  return 0;
}

// Return the partial fraction expansion of the approximation x^(pnum/pden)
int AlgRemez::getPR(double *Root, double *Pole, double *Norm) {

  if (n!=d) {
    printf("Cannot handle case: Numerator degree neq Denominator degree\n");
    return 0;
  }

  if (!alloc) {
    printf("Approximation not yet generated\n");
    return 0;
  }

  if (!foundRoots) {
    printf("Roots not found, so PFE cannot be taken\n");
    return 0;
  }

  // Convert to double and return
  *Norm = (double)norm;
  for (int i=0; i<n; i++) Root[i] = (double)roots[i];
  for (int i=0; i<d; i++) Pole[i] = (double)poles[i];

  // Where the smallest shift is located
  return 0;
}

// Return the partial fraction expansion of the approximation x^(-pnum/pden)
int AlgRemez::getIPFE(double *Res, double *Pole, double *Norm) {

  if (n!=d) {
    printf("Cannot handle case: Numerator degree neq Denominator degree\n");
    return 0;
  }

  if (!alloc) {
    printf("Approximation not yet generated\n");
    return 0;
  }

  if (!foundRoots) {
    printf("Roots not found, so PFE cannot be taken\n");
    return 0;
  }

  bigfloat *r = new bigfloat[d];
  bigfloat *p = new bigfloat[n];
  
  // Want the inverse function
  for (int i=0; i<n; i++) {
    r[i] = poles[i];
    p[i] = roots[i];
  }

  // Perform a partial fraction expansion
  pfe(r, p, (bigfloat)1l/norm);

  // Convert to double and return
  *Norm = (double)((bigfloat)1l/(norm));
  for (int i=0; i<n; i++) {
    Res[i] = (double)r[i];
    Pole[i] = (double)p[i];
  }

  delete [] r;
  delete [] p;

  // Where the smallest shift is located
  return 0;
}

// Initial values of maximal and minimal errors
void AlgRemez::initialGuess() {

  // Supply initial guesses for solution points
  long ncheb = neq;			// Degree of Chebyshev error estimate
  bigfloat a, r;

  // Find ncheb+1 extrema of Chebyshev polynomial

  a = ncheb;
  mm[0] = apstrt;
  double st=(double)apstrt;
  double ld=(log((double)apend)-log((double)apstrt))/double(ncheb);
  for (long i = 1; i < ncheb; i++) {
    /* comment this for the original guess
    r = 0.5 * (1 - cos((M_PI * i)/(double) a));
    r = (exp((double)r)-1.0)/(exp(1.0)-1.0);
    mm[i] = apstrt + r * apwidt;
    /*/
    st = exp(log(st)+ld);
    mm[i]=st;
    //*/
  }
  mm[ncheb] = apend;

  a = 2.0 * ncheb;
  st=(double)apstrt;
  st = exp(log(st)+ld/2.);
  for (long i = 0; i < ncheb; i++) {
    /* comment this for the original guess
    r = 0.5 * (1 - cos(M_PI * (2*i+1)/(double) a));
    r = (exp((double)r)-1.0)/(exp(1.0)-1.0);
    xx[i] = apstrt + r * apwidt;
    /*/
    xx[i]=st;
    st = exp(log(st)+ld);
    //*/
  }
	xx[ncheb]=apend;

}

// Initialise step sizes
void AlgRemez::stpini(bigfloat *step) {
  //xx[neq+1] = apend;
  delta = 0.25;
  step[0] = xx[0] - apstrt;
  for (int i = 1; i < neq; i++) step[i] = xx[i] - xx[i-1];
  step[neq] = step[neq-1];
}

// Search for error maxima and minima
void AlgRemez::search(bigfloat *step) {
  bigfloat a, q, xm, ym, xn, yn, xx0, xx1;
  int i, j, meq, emsign, ensign, steps;

  meq = neq + 1;
  bigfloat *yy = new bigfloat[meq];

  bigfloat eclose = 1.0e30;
  bigfloat farther = 0l;

  j = 1;
  xx0 = apstrt;

  for (i = 0; i < meq; i++) {
    steps = 0;
		//xx1 = xx[i]; // Next zero
    if (i == meq-1) xx1 = apend;
		else xx1 = xx[i];
    xm = mm[i];
    ym = getErr(xm,&emsign);
    q = step[i];
    xn = xm + q;
    if (xn < xx0 || xn >= xx1) {	// Cannot skip over adjacent boundaries
      q = -q;
      xn = xm;
      yn = ym;
      ensign = emsign;
    } else {
      yn = getErr(xn,&ensign);
      if (yn < ym) {
				q = -q;
				xn = xm;
				yn = ym;
				ensign = emsign;
      }
    }
  
    while(yn >= ym) {		// March until error becomes smaller.
      if (++steps > 10) break;
      ym = yn;
      xm = xn;
      emsign = ensign;
      a = xm + q;
      if (a == xm || a <= xx0 || a >= xx1) break;// Must not skip over the zeros either side.
      xn = a;
      yn = getErr(xn,&ensign);
    }

    mm[i] = xm;			// Position of maximum
    yy[i] = ym;			// Value of maximum

    if (eclose > ym) eclose = ym;
    if (farther < ym) farther = ym;

    xx0 = xx1; // Walk to next zero.
  } // end of search loop

  q = (farther - eclose);	// Decrease step size if error spread increased
  if (eclose != 0.0) q /= eclose; // Relative error spread
  /* comment this for the original code
  if (q >= spread) delta *= 0.2; 
  /*/
  if (q >= spread) delta *= 0.2; 
  else delta*=1.5;
  //*/
  spread = q;

  for (i = 0; i < neq; i++) {
    q = yy[i+1];
    if (q != 0.0) q = yy[i] / q  - (bigfloat)1l;
    else q = 0.0625;
    if (q > (bigfloat)0.25) q = 0.25;
    q *= mm[i+1] - mm[i];
    step[i] = q * delta;
  }
  step[neq] = step[neq-1];
  
  for (i = 0; i < neq; i++) {	// Insert new locations for the zeros.
    xm = xx[i] - step[i];
    if (xm <= apstrt) continue;
    if (xm >= apend) continue;
    if (xm <= mm[i]) xm = (bigfloat)0.5 * (mm[i] + xx[i]);
    if (xm >= mm[i+1]) xm = (bigfloat)0.5 * (mm[i+1] + xx[i]);
    xx[i] = xm;
  }

  delete [] yy;
}

// Solve the equations
void AlgRemez::equations(void) {
  bigfloat x, y, z;
  int i, j, ip;
  bigfloat *aa;

  bigfloat *AA = new bigfloat[(neq)*(neq)];
  bigfloat *BB = new bigfloat[neq];
  
  for (i = 0; i < neq; i++) {	// set up the equations for solution by simq()
    ip = neq * i;		// offset to 1st element of this row of matrix
    x = xx[i];			// the guess for this row
    y = func(x);		// right-hand-side vector

    z = (bigfloat)1l;
    aa = AA+ip;
    for (j = 0; j <= n; j++) {
      *aa++ = z;
      z *= x;
    }

    z = (bigfloat)1l;
    for (j = 0; j < d; j++) {
      *aa++ = -y * z;
      z *= x;
    }
    BB[i] = y * z;		// Right hand side vector
  }

  // Solve the simultaneous linear equations.
  if (simq(AA, BB, param, neq)) {
    printf("simq failed\n");
    exit(0);
  }

  delete [] AA;
  delete [] BB;

}

// Evaluate the rational form P(x)/Q(x) using coefficients
// from the solution vector param
bigfloat AlgRemez::approx(const bigfloat x) {
  bigfloat yn, yd;
  int i;

  // Work backwards toward the constant term.
  yn = param[n];		// Highest order numerator coefficient
  for (i = n-1; i >= 0; i--) yn = x * yn  +  param[i]; 
  yd = x + param[n+d];	// Highest degree coefficient = 1.0
  for (i = n+d-1; i > n; i--) yd = x * yd  +  param[i];

  return(yn/yd);
}

// Compute size and sign of the approximation error at x
bigfloat AlgRemez::getErr(bigfloat x, int *sign) {
  bigfloat e, f;

  f = func(x);
  e = approx(x) - f;
  if (f != 0) e /= f;
  if (e < (bigfloat)0.0) {
    *sign = -1;
    e = -e;
  }
  else *sign = 1;
  
  return(e);
}

// Calculate function required for the approximation.
bigfloat AlgRemez::func(const bigfloat x) {

  bigfloat z = (bigfloat)power_num / (bigfloat)power_den;
  bigfloat y;

  if (x == (bigfloat)1.0) y = (bigfloat)1.0;
  else y = pow_bf(x,z);

  if (a_length > 0) {
    bigfloat sum = 0l;
    for (int j=0; j<a_length; j++) sum += a[j]*pow_bf(x,a_power[j]);
    return y * exp_bf(sum);
  } else {
    return y;
  }

}

// Solve the system AX=B
int AlgRemez::simq(bigfloat A[], bigfloat B[], bigfloat X[], int n) {

  int i, j, ij, ip, ipj, ipk, ipn;
  int idxpiv, iback;
  int k, kp, kp1, kpk, kpn;
  int nip, nkp, nm1;
  bigfloat em, q, rownrm, big, size, pivot, sum;
  bigfloat *aa;

  // simq() work vector
  int *IPS = new int[(neq) * sizeof(int)];

  nm1 = n - 1;
  // Initialize IPS and X
  
  ij = 0;
  for (i = 0; i < n; i++) {
    IPS[i] = i;
    rownrm = 0.0;
    for(j = 0; j < n; j++) {
      q = abs_bf(A[ij]);
      if(rownrm < q) rownrm = q;
      ++ij;
    }
    if (rownrm == (bigfloat)0l) {
      printf("simq rownrm=0\n");
      delete [] IPS;
      return(1);
    }
    X[i] = (bigfloat)1.0 / rownrm;
  }
  
  for (k = 0; k < nm1; k++) {
    big = 0.0;
    idxpiv = 0;
    
    for (i = k; i < n; i++) {
      ip = IPS[i];
      ipk = n*ip + k;
      size = abs_bf(A[ipk]) * X[ip];
      if (size > big) {
	big = size;
	idxpiv = i;
      }
    }
    
    if (big == (bigfloat)0l) {
      printf("simq big=0\n");
      delete [] IPS;
      return(2);
    }
    if (idxpiv != k) {
      j = IPS[k];
      IPS[k] = IPS[idxpiv];
      IPS[idxpiv] = j;
    }
    kp = IPS[k];
    kpk = n*kp + k;
    pivot = A[kpk];
    kp1 = k+1;
    for (i = kp1; i < n; i++) {
      ip = IPS[i];
      ipk = n*ip + k;
      em = -A[ipk] / pivot;
      A[ipk] = -em;
      nip = n*ip;
      nkp = n*kp;
      aa = A+nkp+kp1;
      for (j = kp1; j < n; j++) {
	ipj = nip + j;
	A[ipj] = A[ipj] + em * *aa++;
      }
    }
  }
  kpn = n * IPS[n-1] + n - 1;	// last element of IPS[n] th row
  if (A[kpn] == (bigfloat)0l) {
    printf("simq A[kpn]=0\n");
    delete [] IPS;
    return(3);
  }

  
  ip = IPS[0];
  X[0] = B[ip];
  for (i = 1; i < n; i++) {
    ip = IPS[i];
    ipj = n * ip;
    sum = 0.0;
    for (j = 0; j < i; j++) {
      sum += A[ipj] * X[j];
      ++ipj;
    }
    X[i] = B[ip] - sum;
  }
  
  ipn = n * IPS[n-1] + n - 1;
  X[n-1] = X[n-1] / A[ipn];
  
  for (iback = 1; iback < n; iback++) {
    //i goes (n-1),...,1
    i = nm1 - iback;
    ip = IPS[i];
    nip = n*ip;
    sum = 0.0;
    aa = A+nip+i+1;
    for (j= i + 1; j < n; j++) 
      sum += *aa++ * X[j];
    X[i] = (X[i] - sum) / A[nip+i];
  }
  
  delete [] IPS;
  return(0);
}

// Calculate the roots of the approximation
int AlgRemez::root() {

  long i,j;
  bigfloat x,dx=0.05;
  bigfloat upper=1, lower=-100000;
  bigfloat tol = 1e-20;

  bigfloat *poly = new bigfloat[neq+1];

  // First find the numerator roots
  for (i=0; i<=n; i++) poly[i] = param[i];
  //for (i=0; i<=n; i++) printf("%d %e\n", i, (double)poly[i]);

  for (i=n-1; i>=0; i--) {
    roots[i] = rtnewt(poly,i+1,lower,upper,tol);
    //printf("root[%d] = %e\n", i, (double)roots[i]);
    if (roots[i] == 0.0) {
      printf("Failure to converge on root %d/%d\n", i+1, n);
      return 0;
    }
    poly[0] = -poly[0]/roots[i];
    for (j=1; j<=i; j++) poly[j] = (poly[j-1] - poly[j])/roots[i];
  }
  
 // Now find the denominator roots
  poly[d] = 1l;
  for (i=0; i<d; i++) poly[i] = param[n+1+i];
  //for (i=0; i<=d; i++) printf("%d %e\n", i, (double)poly[i]);

  for (i=d-1; i>=0; i--) {
    poles[i]=rtnewt(poly,i+1,lower,upper,tol);
    //printf("pole[%d] = %e\n", i, (double)poles[i]);
    if (poles[i] == 0.0) {
      printf("Failure to converge on pole %d/%d\n", i+1, d);
      return 0;
    }
    poly[0] = -poly[0]/poles[i];
    for (j=1; j<=i; j++) poly[j] = (poly[j-1] - poly[j])/poles[i];
  }

  norm = param[n];
  //printf("Normalisation constant is %e\n",(double)norm);
  //for (i=0; i<n; i++) printf("%ld root = %e\n",i,(double)roots[i]);
  //for (i=0; i<d; i++) printf("%ld pole = %e\n",i,(double)poles[i]);

  delete [] poly;

  return 1;
}

// Evaluate the polynomial
bigfloat AlgRemez::polyEval(bigfloat x, bigfloat *poly, long size) {
  bigfloat f = poly[size];
  for (int i=size-1; i>=0; i--) f = f*x + poly[i];
  return f;
}

/*
// Evaluate the polynomial using complex numbers
complex_bf AlgRemez::polyEval(complex_bf x,complex_bf *poly, long size) {
  complex_bf f = poly[size];
  for (long i=size-1; i>=0; i--) f = f*x + poly[i];
  return f;
}
*/
// Evaluate the differential of the polynomial
bigfloat AlgRemez::polyDiff(bigfloat x, bigfloat *poly, long size) {
  bigfloat df = (bigfloat)size*poly[size];
  for (int i=size-1; i>0; i--) df = df*x + (bigfloat)i*poly[i];
  return df;
}


// Newton's method to calculate roots
bigfloat AlgRemez::rtnewt(bigfloat *poly, long i, bigfloat x1, 
			  bigfloat x2, bigfloat xacc) {
  int j;
  bigfloat df, dx, f, rtn;

  rtn=(bigfloat)0.5*(x1+x2);
  for (j=1; j<=JMAX;j++) {
    f = polyEval(rtn, poly, i);
    df = polyDiff(rtn, poly, i);
    dx = f/df;
    rtn -= dx;
    //if ((x1-rtn)*(rtn-x2) < (bigfloat)0.0)
    //printf("Jumped out of brackets in rtnewt\n");
    if (abs_bf(dx) < xacc) return rtn;
  }
  printf("Maximum number of iterations exceeded in rtnewt\n");
  return 0.0;
}

// Evaluate the partial fraction expansion of the rational function
// with res roots and poles poles.  Result is overwritten on input
// arrays.
void AlgRemez::pfe(bigfloat *res, bigfloat *poles, bigfloat norm) {
  int i,j,small;
  bigfloat temp;
  bigfloat *numerator = new bigfloat[n];
  bigfloat *denominator = new bigfloat[d];

  // Construct the polynomials explicitly 
  for (i=1; i<n; i++) {
    numerator[i] = 0l;
    denominator[i] = 0l;
  }
  numerator[0]=1l;
  denominator[0]=1l;

  for (j=0; j<n; j++) {
    for (i=n-1; i>=0; i--) {
      numerator[i] *= -res[j];
      denominator[i] *= -poles[j];
      if (i>0) {
	numerator[i] += numerator[i-1];
	denominator[i] += denominator[i-1];
      }
    }
  }

  // Convert to proper fraction form.
  // Fraction is now in the form 1 + n/d, where O(n)+1=O(d)
  for (i=0; i<n; i++) numerator[i] -= denominator[i];

  // Find the residues of the partial fraction expansion and absorb the
  // coefficients.
  for (i=0; i<n; i++) {
    res[i] = 0l;
    for (j=n-1; j>=0; j--) {
      res[i] = poles[i]*res[i]+numerator[j];
    }
    for (j=n-1; j>=0; j--) {
      if (i!=j) res[i] /= poles[i]-poles[j];
    }
    res[i] *= norm;
  }  

  // res now holds the residues
  j = 0;
  for (i=0; i<n; i++) poles[i] = -poles[i];

  // Move the ordering of the poles from smallest to largest
  for (j=0; j<n; j++) {
    small = j;
    for (i=j+1; i<n; i++) {
      if (poles[i] < poles[small]) small = i;
    }
    if (small != j) {
      temp = poles[small];
      poles[small] = poles[j];
      poles[j] = temp;
      temp = res[small];
      res[small] = res[j];
      res[j] = temp;
    }
    //printf("Residue = %e, Pole = %e\n", (double)res[j], (double)poles[j]);
  }

  delete [] numerator;
  delete [] denominator;
}

double AlgRemez::evaluateApprox(double x) {
  return (double)approx((bigfloat)x);
}

double AlgRemez::evaluateInverseApprox(double x) {
  return 1.0/(double)approx((bigfloat)x);
}

double AlgRemez::evaluateFunc(double x) {
  return (double)func((bigfloat)x);
}

double AlgRemez::evaluateInverseFunc(double x) {
  return 1.0/(double)func((bigfloat)x);
}
