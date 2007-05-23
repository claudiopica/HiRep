/*
  Mike Clark - 25th May 2005

  Quick test code for calculating optimal rational approximations for
  the functions x^(y/z) with appropriate bounds.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include"alg_remez.h"

int main (int argc, char* argv[]) {

  int i=0;
  int n; // The degree of the numerator polynomial
  int d; // The degree of the denominator polynomial
  int y; // The numerator of the exponent
  int z; // The denominator of the exponent
  int precision; // The precision that gmp uses
  double lambda_low, lambda_high; // The bounds of the approximation

  // Set the exponent
  sscanf(argv[++i],"%d",&y);
  sscanf(argv[++i],"%d",&z);  

  // Set the required degree of approximation
  sscanf(argv[++i],"%d",&n);
  sscanf(argv[++i],"%d",&d);

  // Set the approximation bounds
  sscanf(argv[++i],"%le",&lambda_low);
  sscanf(argv[++i],"%le",&lambda_high);

  // Set the precision of the arithmetic
  sscanf(argv[++i],"%d",&precision);

  // The error from the approximation (the relative error is minimised
  // - if another error minimisation is requried, then line 398 in
  // alg_remez.C is where to change it)
  double error;

  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;
  double *res = new double[n];
  double *pole = new double[d];

  double bulk = exp(0.5*(log(lambda_low)+log(lambda_high)));
  char FORCE_FILE[100], ENERGY_FILE[100];
  sprintf(FORCE_FILE, "force_%d_%d_%d_%d_%f.dat", y, z, d, n, bulk);
  sprintf(ENERGY_FILE, "energy_%d_%d_%d_%d_%f.dat", y, z, d, n, bulk);

  // Instantiate the Remez class
  AlgRemez remez(lambda_low,lambda_high,precision);

  // Generate the required approximation
  error = remez.generateApprox(n,d,y,z);

  FILE *output = fopen("approx.dat", "w");

  fprintf(output, "Approximation Range: [%18.16e,%18.16e]\nError: %18.16e\n", lambda_low, lambda_high, error);

  remez.getPR(res,pole,&norm);
  
  fprintf(output, "Approximation to f(x) = x^(-%d/%d)=a0*(x-a1)/(x-b1)*(x-a2)/(x-b2)*...\n\n", y, z);
  fprintf(output, "app->a[0] = %18.16e;\n", 1./norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "app->a[%d] = %18.16e ; app->b[%d] = %18.16e ;\n", 
	    n-i, pole[i], n-i-1, res[i]);
  }


  // Find the partial fraction expansion of the approximation 
  // to the function x^{y/z} (this only works currently for 
  // the special case that n = d)
  fprintf(output, "Approximation to f(x) = x^(%d/%d)\n\n", y, z);

  remez.getPFE(res,pole,&norm);
  
  fprintf(output, "constant min_epsilon = %18.16e\n\n", lambda_low);
  fprintf(output, "Ra_a0 = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "Ra_a[%d] = %18.16e ; Ra_b[%d] = %18.16e ;\n", 
	    n-i-1, res[i], n-i-1, -pole[i]);
  }
/*
  fprintf(output, "alpha[0] = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	    i+1, res[i], i+1, pole[i]);
  }
*/
  // Find pfe of inverse function
  remez.getIPFE(res,pole,&norm);
  fprintf(output, "\nApproximation to f(x) = x^(-%d/%d)\n\n", y, z);
  fprintf(output, "constant min_epsilon = %18.16e\n\n", lambda_low);
  fprintf(output, "Ra_a0 = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "Ra_a[%d] = %18.16e ; Ra_b[%d] = %18.16e ;\n", 
	    n-i-1, res[i], n-i-1, -pole[i]);
  }

  fclose(output);

  FILE *error_file = fopen("error.dat", "w");
  for (double x=lambda_low; x<lambda_high; x*=1.01) {
    double f = remez.evaluateFunc(x);
    double r = remez.evaluateApprox(x);
    fprintf(error_file,"%e %e\n", x,  (r - f)/f);
  }
  fclose(error_file);

  delete res;
  delete pole;

  exit(0);

}     
