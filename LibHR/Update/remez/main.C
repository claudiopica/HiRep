/*
	Modified by C. Pica - 10 May 2007
  Mike Clark - 25th May 2005

  Quick test code for calculating optimal rational approximations for
  the functions x^(y/z) with appropriate bounds.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include"alg_remez.h"

typedef struct _guess_inf {
	int c;
	int n1,n2;
	double l1,l2;
} guess_inf;

int main (int argc, char* argv[]) {

	int i=0;
	int n; // The degree of the numerator polynomial
	int d; // The degree of the denominator polynomial
	int y; // The numerator of the exponent
	int z; // The denominator of the exponent
	int precision; // The precision that gmp uses
	double lambda_low, lambda_high; // The bounds of the approximation
	double reqprec, old_low, lowlimit, uplimit;
	double min_lambda;

	guess_inf low_inf;
	low_inf.c=0; /* init struct */

	// Set the exponent
	sscanf(argv[++i],"%d",&y);
	sscanf(argv[++i],"%d",&z);  

	//precision required by the approximations
	sscanf(argv[++i],"%le",&reqprec);

	// Set the required degree of approximation
	sscanf(argv[++i],"%d",&n);
	//sscanf(argv[++i],"%d",&d);

	// Set the approximation bounds
	sscanf(argv[++i],"%le",&min_lambda);
	sscanf(argv[++i],"%le",&lambda_high);

	// Set the precision of the arithmetic
	sscanf(argv[++i],"%d",&precision);

	// The partial fraction expansion takes the form 
	// r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
	double norm,oldnorm, olderror;

	//lambda_low=(min_lambda+lambda_high)*.5;
	old_low=lambda_low=2.*min_lambda;
	lowlimit=0;
	uplimit=lambda_high*.999;
	while (old_low>min_lambda) {
		//reqprec=1.e-6;
		d=n; //numerator and den equal
		//lowlimit=0.;
		//lambda_low*=.25*10.;
		//uplimit=lambda_low;

		double *res = new double[n];
		double *pole = new double[d];

		// The error from the approximation (the relative error is minimised
		// - if another error minimisation is requried, then line 398 in
		// alg_remez.C is where to change it)
		double error=0.;
		olderror=0.;
		bigfloat *dmm=0;
		while ((fabs(reqprec-error)/reqprec)>.01 || reqprec<error) {

			
			// Instantiate the Remez class
			AlgRemez remez(lambda_low,lambda_high,precision);

			// Generate the required approximation
			error = remez.generateApprox(n,d,y,z,dmm);
			if (dmm==0) dmm=new bigfloat[2*(n+d)+4];
			if (dmm!=0) remez.getMM(dmm);
			//printf("Limits: %e, %e - Error: %e\n",lowlimit,uplimit,error);


			if (error>reqprec){
				lowlimit=lambda_low;
			} else {
				uplimit=lambda_low;
			}

			if (olderror==0.){
				old_low=lambda_low;
				olderror=error;
				lambda_low=(lowlimit==0.)?0.5*(lowlimit+uplimit):sqrt(lowlimit*uplimit);
			} else {
				double tmp=lambda_low;
				lambda_low+=(old_low-lambda_low)*(log(reqprec*.999/error)/log(olderror/error));
				if (lambda_low<lowlimit || lambda_low>uplimit)
					lambda_low=(lowlimit==0.)?0.5*(lowlimit+uplimit):sqrt(lowlimit*uplimit);
				old_low=tmp;
				olderror=error;
			}

		}
		lambda_low=old_low;

		AlgRemez remez(lambda_low,lambda_high,precision);
		error = remez.generateApprox(n,d,y,z,dmm);
		remez.getPR(res,pole,&norm);

		/* print out solution */
		printf("%d,%d,%d",y,z,n); /* num,den,degree */
		printf(",%18.16e", error); /* error */
		printf(",%18.16e,%18.16e", lambda_low, lambda_high); /* range */
		printf(",%18.16e", 1./norm); /* a[0] */
		for (int i = 0; i < n; i++)
			printf(",%18.16e",pole[n-1-i]);
		for (int i = 0; i < n; i++)
			printf(",%18.16e",res[n-1-i]);
		printf("\n");
		fflush(stdout);


		FILE *output = fopen("approx.dat", "w");

		fprintf(output, "Approximation Range: [%18.16e,%18.16e]\nError: %18.16e\n", lambda_low, lambda_high, error);

		//remez.getPR(res,pole,&norm);


		fprintf(output, "Approximation to f(x) = x^(-%d/%d)=a0*(x-a1)/(x-b1)*(x-a2)/(x-b2)*...\n\n", y, z);
		fprintf(output, "app->a[0] = %18.16e;\n", 1./norm);
		for (int i = 0; i < n; i++) {
			fprintf(output, "app->a[%d] = %18.16e ; app->b[%d] = %18.16e ;\n", 
					n-i, pole[i], n-i-1, res[i]);
		}


		/*
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
		*/
		/* 
		// Find pfe of inverse function
		remez.getIPFE(res,pole,&norm);
		fprintf(output, "\nApproximation to f(x) = x^(-%d/%d)\n\n", y, z);
		fprintf(output, "constant min_epsilon = %18.16e\n\n", lambda_low);
		fprintf(output, "Ra_a0 = %18.16e\n", norm);
		for (int i = 0; i < n; i++) {
		fprintf(output, "Ra_a[%d] = %18.16e ; Ra_b[%d] = %18.16e ;\n", 
		n-i-1, res[i], n-i-1, -pole[i]);
		}
		*/

		fclose(output);

		/*
		// error check file
		FILE *error_file = fopen("error.dat", "w");
		for (double x=lambda_low; x<lambda_high; x*=1.01) {
		double f = remez.evaluateFunc(x);
		double r = remez.evaluateApprox(x);
		fprintf(error_file,"%e %e\n", x,  (r - f)/f);
		}
		fclose(error_file);
		*/

		if (dmm!=0) delete[] dmm;
		delete[] res;
		delete[] pole;

		/* guess new limits for next n*/
		if(low_inf.c==0){
			/* first point: put it in n1,l1 */
			low_inf.n1=n;
			low_inf.l1=log(lambda_low);
			/* leave lambda_low and uplimit the same */
			/* but lower lowlimit */
			lowlimit=0.;
		} else {
			low_inf.n2=n;
			low_inf.l2=log(lambda_low);
			/* guess new interval */
			/* the new guess is based on a linear interpolation in log scale */
			double delta=(low_inf.l1-low_inf.l2)/(double)(low_inf.n2-low_inf.n1); /* this is positive */
			double newlw = low_inf.l2-delta;
			lambda_low=exp(newlw);
			lowlimit=exp(newlw-0.05*delta);
			uplimit=exp(newlw+0.05*delta);
			//printf("Guess: %e %e %e\n",lowlimit,lambda_low,uplimit);
			//printf("Guess: %e %e %e %d %d\n",newlw,low_inf.l1,low_inf.l2,low_inf.n1,low_inf.n2);
			low_inf.n1=low_inf.n2;
			low_inf.l1=low_inf.l2;
		}

		++low_inf.c;
		++n;
	}

	exit(0);

}     
