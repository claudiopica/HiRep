#include "inverters.h"
#include "linear_algebra.h"
#include "global.h"
#include "utils.h"
#include "memory.h"
#include "logger.h"
#include "update.h"
#include "inverters.h"
#include <stdlib.h>
#include <math.h>

static spinor_field *ev=NULL;
static double *eigval=NULL;
static spinor_operator H_loc=NULL;

static void free_loc_mem() {
  if (ev!=NULL) {free_spinor_field_f(ev); ev=NULL; }
  if (eigval!=NULL) {free(eigval); eigval=NULL; }
}


static eva_prec loc_par={0};

static int alloc_loc_mem(unsigned int nevt, geometry_descriptor *type) {
  static int init=1;
  static unsigned int oldnevt=0;
  static geometry_descriptor *oldtype=NULL;
  if (init) { /* register memory clean-up function */
    atexit(&free_loc_mem);
    init=0;
  }
  
  if (nevt<1) { return 1; } 
  if (nevt!=oldnevt || type!=oldtype) {
    free_loc_mem();
    if (ev==NULL) {
      ev=alloc_spinor_field_f(nevt,type);
    }
    if (eigval==NULL) {
      eigval=malloc(sizeof(*eigval)*nevt);
    }
    /* update old nevt */
    oldnevt=nevt;
    oldtype=type;

    return 1;
  }

  return 0;
}

static void find_low_eig_H2(const eva_prec *e_par, geometry_descriptor *type, int start) {

  int n;

  /* EVA parameters */
  double max;
  int status,ie;
  /* END of EVA parameters */
  int MVM=0; /* counter for matrix-vector multiplications */

  lprintf("EVA_PREC",0,"Starting EVA\n");

  max_H(H_loc, type, &max);
  max*=1.1;

  ie=eva(e_par->nev,e_par->nevt,start,e_par->kmax,e_par->maxiter,max,e_par->omega1,e_par->omega2,H_loc,ev,eigval,&status);
  MVM+=status;
  while (ie!=0) { /* if failed restart EVA */
    lprintf("EVA_PREC",0,"Restarting EVA!\n");
    ie=eva(e_par->nev,e_par->nevt,2,e_par->kmax,e_par->maxiter,max,e_par->omega1,e_par->omega2,H_loc,ev,eigval,&status);
    MVM+=status;
  }

  lprintf("EVA_PREC",0,"EVA MVM = %d\n",MVM);
  for (n=0;n<e_par->nev;++n) {
    lprintf("EVA_PREC",0,"Eig %d = %1.15e\n",n,eigval[n]);
  }

}

static void check_ortho(spinor_field *in, int n) {
  int i,j;
  complex p = {0,0};
  for (i=0;i<n;++i) {
    for (j=i;j<n;++j) {
      p=spinor_field_prod_f(&in[i],&in[j]);
      lprintf("TESTORTHO",0,"(%e,%e) ",p.re,p.im);
    }
    lprintf("TESTORTHO",0,"\n",p.re,p.im);
  }
}

static void orthogonalize(spinor_field *out, spinor_field *in, int n) {
  complex p;
  while(n>0) {
    --n;
    p=spinor_field_prod_f(&in[n],out);
    _complex_minus(p,p);
    spinor_field_mulc_add_assign_f(out,p,&in[n]);
  }
}

/*static void GS_vect(spinor_field *b, int n) {
  int i;
  double inorm;
  for (i=0;i<n;++i) {
    orthogonalize(&b[i], &b[i+1], n-i-1); 
    inorm=1./sqrt(spinor_field_sqnorm_f(&b[i]));
    spinor_field_mul_f(&b[i],inorm,&b[i]);
  }
}*/

void set_def_matrix(eva_prec *e_par, spinor_operator H, geometry_descriptor *type) {
  int changed;

  lprintf("EVA_PREC",0,"Setting new preconditioning matrix.\n");

  loc_par=*e_par;
  H_loc=H;

  changed=alloc_loc_mem(e_par->nevt, type);
  if (changed) 
    find_low_eig_H2(e_par,type,0);
  else 
    find_low_eig_H2(e_par,type,2);

  check_ortho(ev,e_par->nev);
}


void eva_def(spinor_field *out, spinor_field *in){
  spinor_field_copy_f(out,in);
  orthogonalize(out,ev,loc_par.nev);
  /*
  for (i=0;i<loc_par.nev;++i) {
    complex p=spinor_field_prod_f(&ev[i],out);
    _complex_mulr(p,(1./eigval[i]-1.),p);
    spinor_field_mulc_add_assign_f(out,p,&ev[i]);
  }
  */
}

void eva_def_inv(spinor_field *out, spinor_field *in, double m){
  int i;
  spinor_field_zero_f(out);
  for (i=0;i<loc_par.nev;++i) {
    complex p=spinor_field_prod_f(&ev[i],in);
    _complex_mulr(p,1./(eigval[i]-m),p);
    spinor_field_mulc_add_assign_f(out,p,&ev[i]);
  }
}










