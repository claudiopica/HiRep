#include "global.h"
#include "inverters.h"
#include "linear_algebra.h"
#include "rational_functions.h"
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>

/* Functions for coefficients computation */

/* converts the coef from the root/poles form to the partial fraction exp
 * before : r(x)=a[0]*(x-a[1])/(x-b[0])*...*(x-a[n+1])/(x-b[n])
 * after  : r(x)=a[0]+a[1]/(x-b[0])+a[2]/(x-b[1])+...+a[n+1)/(x-b[n])
 */
void r_app_rp2pfe(rational_app *app) {
  int i,j;
  double *ctmp;

  ctmp=(double*)malloc(sizeof(double)*(app->order));
  for(i=0;i<app->order;++i)
    ctmp[i]=app->a[i+1];

  /* only a[1]...a[n+1] changes */
  for(i=0;i<app->order;++i) {
    app->a[i+1]=app->a[0];
    for(j=0;j<app->order;++j) {
      app->a[i+1]*=(i==j)?(app->b[i]-ctmp[j]):((app->b[i]-ctmp[j])/(app->b[i]-app->b[j]));
    }
  }

  free(ctmp);
  
}

/* converts the coef from the root/poles form to the partial fraction exp
 * as r_app_rp3pfe but for the inverse function
 */
void r_app_rp2pfe_inv(rational_app *app) {
  int i;
  double ctmp;
  
  app->a[0]=1./app->a[0];
  /* swap poles and roots: a[i+1] <-> b[i] */
  for(i=0;i<app->order;++i) {
    ctmp=app->a[i+1];
    app->a[i+1]=app->b[i];
    app->b[i]=ctmp;
  }
  
  r_app_rp2pfe(app);

}

/* approximation for f(x)=x^(-1/2) */
/* the order to use is written into app */
void inv_sqrt_rp(rational_app *app, double min, double max) {
  if (app->order==7) {
    /*
     * Approximation Range: [9.5000000000000000e-04,6.4000000000000000e+01]
     * Error: 9.4054878526558268e-05
     * Approximation to f(x) = x^(-1/2)=a0*(x-a1)/(x-b1)*(x-a2)/(x-b2)*...
    */
    if(min<9.5e-04 || max>64.)
      printf("ERRORE: limiti dell'approssimazione rationale o7 non validi!\n");

    app->error = 9.4054878526558268e-05;

    app->a[0] = 3.6845866964412707e-02;
    app->a[7] = -1.0758682836072590e-03 ; app->b[6] = -2.1864431996069619e-04 ;
    app->a[6] = -9.1776297443641811e-03 ; app->b[5] = -3.3608858903202148e-03 ;
    app->a[5] = -6.1029648667047832e-02 ; app->b[4] = -2.3883202443314405e-02 ;
    app->a[4] = -3.9250133337802012e-01 ; app->b[3] = -1.5490393236815636e-01 ;
    app->a[3] = -2.5457222558116230e+00 ; app->b[2] = -9.9623709668884586e-01 ;
    app->a[2] = -1.8090468401534210e+01 ; app->b[1] = -6.6248041916635634e+00 ;
    app->a[1] = -2.7807719867101741e+02 ; app->b[0] = -5.6512494072364333e+01 ;
    
    return;
    
  }
  if (app->order==15) {
    /*
     * Approximation Range: [6.0000000000000002e-05,6.4000000000000000e+01]
     * Error: 4.1956376084685596e-08
     * Approximation to f(x) = x^(-1/2)=a0*(x-a1)/(x-b1)*(x-a2)/(x-b2)*...
    */
    if(min<6.0000000000000002e-05 || max>64.)
      printf("ERRORE: limiti dell'approssimazione rationale o15 non validi!\n");

    app->error = 4.1956376084685596e-08;

    app->a[0] = 2.1373791678827175e-02; 
    app->a[15] = -1.9044673685788960e-05 ; app->b[14] = -4.4335604109217891e-06 ;
    app->a[14] = -1.0035871661876729e-04 ; app->b[13] = -4.8151963824264967e-05 ;
    app->a[13] = -3.4718264377796325e-04 ; app->b[12] = -1.9109579329923458e-04 ;
    app->a[12] = -1.0729003104855976e-03 ; app->b[11] = -6.1475465420512477e-04 ;
    app->a[11] = -3.1989603817244527e-03 ; app->b[10] = -1.8570400672186852e-03 ;
    app->a[10] = -9.4250925987057216e-03 ; app->b[9] = -5.4953456178144025e-03 ;
    app->a[9] = -2.7659614018772424e-02 ; app->b[8] = -1.6150166785026270e-02 ;
    app->a[8] = -8.1082569645470920e-02 ; app->b[7] = -4.7359130535578615e-02 ;
    app->a[7] = -2.3776844233957270e-01 ; app->b[6] = -1.3883057071562224e-01 ;
    app->a[6] = -6.9877315587790778e-01 ; app->b[5] = -4.0742305285439001e-01 ;
    app->a[5] = -2.0678067575306667e+00 ; app->b[4] = -1.2003899835514640e+00 ;
    app->a[4] = -6.2463943521746979e+00 ; app->b[3] = -3.5790836878982781e+00 ;
    app->a[3] = -2.0094633867669664e+01 ; app->b[2] = -1.1060460736786801e+01 ;
    app->a[2] = -7.9747526269425578e+01 ; app->b[1] = -3.8262745174263344e+01 ;
    app->a[1] = -8.6612105037306117e+02 ; app->b[0] = -2.0163117853079254e+02 ;
 
    return;
 
  }
  
  printf("ERRORE: limiti ordine dell'approssimazione non valido!\n");
  
}

/* approximation for f(x)=x^(-1/4) */
/* the order to use is written into app */
void inv_fourrt_rp(rational_app *app, double min, double max) {
  if (app->order==7) {
    /*
     * Approximation Range: [9.5000000000000000e-04,6.4000000000000000e+01]
     * Error: 6.4864095405253665e-05
     * Approximation to f(x) = x^(-1/4)=a0*(x-a1)/(x-b1)*(x-a2)/(x-b2)*...
     */
    if(min<9.5e-04 || max>64.)
      printf("ERRORE: limiti dell'approssimazione rationale o7 non validi!\n");

    app->error = 6.4864095405253665e-05 ;

    app->a[0] = 2.0083197112875492e-01;
    app->a[7] = -7.9528771002868569e-04 ; app->b[6] = -3.7311997022139724e-04 ;
    app->a[6] = -7.2550594001712194e-03 ; app->b[5] = -4.4051077626231167e-03 ;
    app->a[5] = -4.8506794567469737e-02 ; app->b[4] = -3.0386139354632986e-02 ;
    app->a[4] = -3.1091647021465657e-01 ; app->b[3] = -1.9555091423115584e-01 ;
    app->a[3] = -2.0009123005200009e+00 ; app->b[2] = -1.2534326488102867e+00 ;
    app->a[2] = -1.3802159510348808e+01 ; app->b[1] = -8.3803586774996131e+00 ;
    app->a[1] = -1.6295027029489540e+02 ; app->b[0] = -7.6450320095864399e+01 ;
    
    return;
    
  }
  if (app->order==15) {
    /*
     * Approximation Range: [6.0000000000000002e-05,6.4000000000000000e+01]
     * Error: 2.9310824640342691e-08
     * Approximation to f(x) = x^(-1/4)=a0*(x-a1)/(x-b1)*(x-a2)/(x-b2)*...
     */
    if(min<6.0000000000000002e-05 || max>64.)
      printf("ERRORE: limiti dell'approssimazione rationale o15 non validi!\n");

    app->error = 2.9310824640342691e-08;

    app->a[0] = 1.5294989062985581e-01;
    app->a[15] = -1.4627780103698447e-05 ; app->b[14] = -7.3465504148420134e-06 ;
    app->a[14] = -8.5086123908935052e-05 ; app->b[13] = -5.9123781517097357e-05 ;
    app->a[13] = -3.0124428938194797e-04 ; app->b[12] = -2.2367734978218326e-04 ;
    app->a[12] = -9.3716106996227994e-04 ; app->b[11] = -7.0964663174485118e-04 ;
    app->a[11] = -2.7990970575664780e-03 ; app->b[10] = -2.1331723166207149e-03 ;
    app->a[10] = -8.2479417100448787e-03 ; app->b[9] = -6.2991968440610138e-03 ;
    app->a[9] = -2.4194372260459375e-02 ; app->b[8] = -1.8490960383838913e-02 ;
    app->a[8] = -7.0877325535069272e-02 ; app->b[7] = -5.4178116499331137e-02 ;
    app->a[7] = -2.0766904045482454e-01 ; app->b[6] = -1.5871459522327327e-01 ;
    app->a[6] = -6.0960152461665273e-01 ; app->b[5] = -4.6557070054500987e-01 ;
    app->a[5] = -1.8001358681061330e+00 ; app->b[4] = -1.3718709716119952e+00 ;
    app->a[4] = -5.4111438400804639e+00 ; app->b[3] = -4.0974813434733877e+00 ;
    app->a[3] = -1.7167585380188868e+01 ; app->b[2] = -1.2747129606600641e+01 ;
    app->a[2] = -6.4948484374084785e+01 ; app->b[1] = -4.5130743105771622e+01 ;
    app->a[1] = -5.2269429639278928e+02 ; app->b[0] = -2.6251420056752875e+02 ;    
 
    return;
 
  }
  
  printf("ERRORE: limiti ordine dell'approssimazione non valido!\n");
  
}

void inv_sqrt_coef(rational_app *app, double min, double max) {
  inv_sqrt_rp(app, min, max);
  r_app_rp2pfe(app);
  /*
  {
    int i;
    printf("a[0]= %1.15e\n",app->a[0]);
    for (i=0;i<app->order;++i){
      printf("a[%d]= %1.15e; ",i+1,app->a[i+1]);
      printf("b[%d]= %1.15e;\n",i,app->b[i]);
    }
    }*/
}

void sqrt_coef(rational_app *app, double min, double max) {
  inv_sqrt_rp(app, min, max);
  r_app_rp2pfe_inv(app);  
}

void inv_fourrt_coef(rational_app *app, double min, double max) {
  inv_fourrt_rp(app, min, max);
  r_app_rp2pfe(app);  
}

void fourrt_coef(rational_app *app, double min, double max) {
  inv_fourrt_rp(app, min, max);
  r_app_rp2pfe_inv(app);  
}

/*
 * computes: out = (a[0]*I + \sum_i a[i]*(Q-b[i])^-1 ) in
 * this MUST work in the case: out==in 
 * where Q is a linear operator acting on spinor in
 * This implementation uses CG_mshift => Q must be hertian positive definite!
 */
void rational_func(rational_app *coef, spinor_operator Q, suNf_spinor *out, suNf_spinor *in) {
   static cg_mshift_par cg_par;
   suNf_spinor **cg_out;
   int i;

   /* allocate cg output vectors */
   cg_out = (suNf_spinor **)malloc(sizeof(suNf_spinor*)*(coef->order));
   for (i=0; i<(coef->order); ++i) {
      cg_out[i] = (suNf_spinor *)malloc(sizeof(suNf_spinor)*VOLUME);
   }

   /* set up cg parameters */
   cg_par.n = coef->order;
   cg_par.shift = coef->b;
   cg_par.err2=coef->error*0.1;    /* CAMBIARE: METTERE PARAMETRI COMUNI ALL'UPDATE */
   cg_par.max_iter=1000;
   
   /* compute inverse vectors */
   cg_mshift(&cg_par, Q, in, cg_out);

   /* sum all the contributions */
   spinor_field_mul_f(out,coef->a[0],in);
   for (i=1; i<(coef->order)+1; ++i) {
      spinor_field_mul_add_assign_f(out,coef->a[i],cg_out[i-1]);
   }

   /* free memory */
   for (i=0; i<(coef->order); ++i) {
      free(cg_out[i]);
   }
   free(cg_out);

}
