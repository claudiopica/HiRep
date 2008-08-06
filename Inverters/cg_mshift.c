#include "inverters.h"
#include "linear_algebra.h"
#include "malloc.h"
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "memory.h"
#include "logger.h"

/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
static int cg_mshift_core(short int *sflags, mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out){

   suNf_spinor *k,*r,*Mk;
   suNf_spinor **p;
   double omega, oldomega, gamma;
   double alpha, lambda, delta;
   double innorm2;
   double *z1, *z2, *z3;

   int i;
   int cgiter;
   unsigned short notconverged;
	 unsigned int spinorlen;
   
   /* fare qualche check sugli input */
   /*
   printf("numero vettori n=%d\n",par->n);
   for (i=0; i<(par->n); ++i) {
      printf("shift[%d]=%f\n",i,par->shift[i]);
      printf("out[%d]=%p\n",i,out[i]);      
   }
   */
   
   /* allocate spinors fields and aux real variables */
	 get_spinor_len(&spinorlen);
   p = malloc(sizeof(suNf_spinor*)*(par->n));
   p[0] = alloc_spinor_field_f(3+par->n);
   for (i=1; i<(par->n); ++i) {
      p[i] = p[i-1]+spinorlen;
   }
   k=p[par->n-1]+spinorlen;
   r=k+spinorlen;
   Mk=r+spinorlen;

   z1 = malloc(sizeof(*z1)*(par->n));
   z2 = malloc(sizeof(*z2)*(par->n));
   z3 = malloc(sizeof(*z3)*(par->n));
   
   /* init recursion */
   cgiter = 0;
   omega = 1.;
   gamma = 0.;
   innorm2=spinor_field_sqnorm_f(in);
	 if(par->n==1) { /* non multishift case */
		 /* use out[0] as initial guess */
		 M(Mk,out[0]);
		 ++cgiter;
		 spinor_field_mul_add_assign_f(Mk,-par->shift[0],out[0]);
		 spinor_field_sub_f(r,in,Mk);
			
	 } else { /* initial guess = 0 for multishift */
	   spinor_field_copy_f(r, in);
	 }
   spinor_field_copy_f(k, r);
	 delta=spinor_field_sqnorm_f(r);
   for (i=0; i<(par->n); ++i) {
      z1[i]=z2[i]=1.;
      spinor_field_copy_f(p[i], r);
      if(par->n!=1) spinor_field_zero_f(out[i]);
      sflags[i]=1;
   }
   
   /* cg recursion */
   do {
      M(Mk,k);
      alpha = spinor_field_prod_re_f(k,Mk);
      oldomega = omega;
      omega = - delta/alpha;
      for (i=0; i<(par->n); ++i) {
         if(sflags[i]) {
            z3[i] = oldomega*z1[i]*z2[i]/(omega*gamma*(z1[i]-z2[i])+z1[i]*oldomega*(1.+par->shift[i]*omega));
            spinor_field_mul_add_assign_f(out[i],-omega*z3[i]/z2[i],p[i]);
         }
      }
      spinor_field_mul_add_assign_f(r,omega,Mk);
      lambda=spinor_field_sqnorm_f(r);
      gamma=lambda/delta;
      delta=lambda;
      
      spinor_field_mul_f(k,gamma,k);
      spinor_field_add_assign_f(k,r);
      notconverged=0; /* assume that all vectors have converged */
      for (i=0; i<(par->n); ++i) {
         /* check convergence of vectors */
         if(delta*z3[i]<par->err2*innorm2) sflags[i]=0;
         if(sflags[i]){
            notconverged++;
            spinor_field_mul_f(p[i],gamma*z3[i]*z3[i]/(z2[i]*z2[i]),p[i]);
            spinor_field_mul_add_assign_f(p[i],z3[i],r);
            z1[i]=z2[i];
            z2[i]=z3[i];
         }
      }

      /* Uncomment this to print cg recursion parameters
      printf("[ %d ] alpha=%e\n",cgiter,alpha);
      printf("[ %d ] omega=%e\n",cgiter,omega);
      printf("[ %d ] still runnning=%d\n",cgiter,notconverged);
      for (i=0;i<par->n;++i) printf("z3[%d]=%e; ",i,z3[i]);
      printf("\n[ %d ] gamma=%e\n",cgiter,gamma);
      printf("[ %d ] delta=%e\n",cgiter,delta);
      */
      
      ++cgiter;
   } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

   /* test results */
   for(i=0;i<par->n;++i){
     double norm;
     M(Mk,out[i]);
		 ++cgiter;
     spinor_field_mul_add_assign_f(Mk,-par->shift[i],out[i]);
     spinor_field_sub_f(Mk,Mk,in);
     norm=spinor_field_sqnorm_f(Mk)/spinor_field_sqnorm_f(in);
		 sflags[i]=1;
     if (fabs(norm)>par->err2){
			 sflags[i]=0;
       lprintf("INVERTER",30,"CG failed on vect %d: err2 = %1.8e > %1.8e\n",i,norm,par->err2);
		 } else {
       lprintf("INVERTER",20,"CG inversion: err2 = %1.8e < %1.8e\n",norm,par->err2);
		 }
   }
   
   /* free memory */
   free_field(p[0]);
   free(p);
   free(z1); free(z2); free(z3);

   /* return number of cg iter */
   return cgiter;
}

int cg_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out){ 
	int cgiter,msiter;
	int i;
	mshift_par par_save=*par;
	short int *sflags = malloc(sizeof(*sflags)*(par->n));

	if (par->n==1) spinor_field_zero_f(out[0]);
	cgiter=cg_mshift_core(sflags, par, M, in, out);
	msiter=cgiter; /* save multishift iterations for logging */

	par->n=1;
	for(i=0;i<par_save.n;++i) {
		int rep=0;
		while (sflags[i]==0) {
			par->shift=par_save.shift+i;
			cgiter+=cg_mshift_core(sflags+i, par, M, in, out+i);
			if((++rep)%5==0)
				lprintf("INVERTER",-10,"CG_mshift recursion = %d (precision too high?)\n",rep);
			if((rep)>10) {
				lprintf("INVERTER",-10,"CG_mshift failed on vector %d (precision too high?)\n",i);
				sflags[i]=1;
			}
		}
	}

	*par=par_save;

	free(sflags);

	lprintf("INVERTER",10,"CG_mshift: MVM = %d/%d\n",msiter,cgiter);

	return cgiter;
}
