#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "dirac.h"
#include "inverters.h"
#include "rational_functions.h"
#include "representation.h"

#include <malloc.h>
#include <stdio.h>

/* declared in update_rhmc.c */
extern rhmc_par _update_par;
extern suNf_spinor *pf;
extern rational_app r_MD; /* used in the action MD evolution */

#define _print_avect(a) printf("(%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e)\n",(a).c1,(a).c2,(a).c3,(a).c4,(a).c5,(a).c6,(a).c7,(a).c8)

#define _print_mat(a) printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.re,(a).c1_2.re,(a).c1_3.re,(a).c2_1.re,(a).c2_2.re,(a).c2_3.re,(a).c3_1.re,(a).c3_2.re,(a).c3_3.re);printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.im,(a).c1_2.im,(a).c1_3.im,(a).c2_1.im,(a).c2_2.im,(a).c2_3.im,(a).c3_1.im,(a).c3_2.im,(a).c3_3.im)

/* we need to compute  Tr  U(x,mu) g_5*(1-g_mu) chi2 # chi1^+
 * where # denotes the tensor product and Tr is the trace on Lorentz space
 * the strategy is the following:
 * given the form of g_5(1-g_mu) one can compute only the first two lorentz
 * components of the spinor; so we first apply g_5(1-g_mu) to chi2 to find the first
 * two components; then we multiply these two vectors by U(x,mu) and
 * store the result in p.c1, p.c2; when computing the trace we can factorize p.c1 and p.c2
 * as they both multiply two components of chi1^+; we store these factors in p.c3 and p.c4.
 * the tensor product is perform by the macro 
 * _suNf_FMAT(u,p): u = p.c1 # p.c3^+ + p.c2 # p.c4^+
 */

/* these macros use the variables ptmp, p */
#define _F_DIR0(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c1,(chi2)->c3);		      \
  _suNf_multiply(p.c1,*(pu_gauge_f(x,0)),ptmp);		      \
  _vector_add_f(ptmp,(chi2)->c2,(chi2)->c4);		      \
  _suNf_multiply(p.c2,*(pu_gauge_f(x,0)),ptmp);		      \
  _vector_sub_f(p.c3,(chi1)->c1,(chi1)->c3);	      \
  _vector_sub_f(p.c4,(chi1)->c2,(chi1)->c4);	      \
  _suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c1,(chi2)->c4);		      \
  _suNf_multiply(p.c1,*(pu_gauge_f(x,1)),ptmp);		      \
  _vector_i_add_f(ptmp,(chi2)->c2,(chi2)->c3);		      \
  _suNf_multiply(p.c2,*(pu_gauge_f(x,1)),ptmp);		      \
  _vector_i_sub_f(p.c3,(chi1)->c1,(chi1)->c4);	      \
  _vector_i_sub_f(p.c4,(chi1)->c2,(chi1)->c3);	      \
  _suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c1,(chi2)->c4);		      \
  _suNf_multiply(p.c1,*(pu_gauge_f(x,2)),ptmp);		      \
  _vector_sub_f(ptmp,(chi2)->c2,(chi2)->c3);		      \
  _suNf_multiply(p.c2,*(pu_gauge_f(x,2)),ptmp);		      \
  _vector_sub_f(p.c3,(chi1)->c1,(chi1)->c4);	      \
  _vector_add_f(p.c4,(chi1)->c2,(chi1)->c3);	      \
  _suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c1,(chi2)->c3);		      \
  _suNf_multiply(p.c1,*(pu_gauge_f(x,3)),ptmp);		      \
  _vector_i_sub_f(ptmp,(chi2)->c2,(chi2)->c4);		      \
  _suNf_multiply(p.c2,*(pu_gauge_f(x,3)),ptmp);		      \
  _vector_i_sub_f(p.c3,(chi1)->c1,(chi1)->c3);	      \
  _vector_i_add_f(p.c4,(chi1)->c2,(chi1)->c4);	      \
  _suNf_FMAT((u),p)



void Force_rhmc_f(float dt, suNg_algebra_vector *force){
  int i, n;
  static suNg_algebra_vector f;
  static suNf_vector ptmp;
  static suNf_spinor p;
  static suNf s1;
  static cg_mshift_par cg_par;
  suNf_spinor **chi;
  suNf_spinor *Hchi;

  /* allocate spinors */
  chi = (suNf_spinor **)malloc(sizeof(suNf_spinor*)*(r_MD.order));
  chi[0] = (suNf_spinor *)malloc(sizeof(suNf_spinor)*((r_MD.order+1)*VOLUME));
  for (i=1; i<(r_MD.order); ++i) {
    chi[i]=chi[i-1]+VOLUME;
  } 
  Hchi = chi[r_MD.order-1]+VOLUME;

  /* Compute (H^2-b[n])^-1 * pf */
  /* set up cg parameters */
  cg_par.n = r_MD.order;
  cg_par.shift = r_MD.b;
  cg_par.err2= r_MD.error;
  cg_par.max_iter=1000;
    
  /* compute inverse vectors chi[i] = (H^2 - b[i])^1 * pf */
  cg_mshift(&cg_par, &H2, pf, chi);

  for (n=0; n<r_MD.order; ++n) {

    g5Dphi(_update_par.mass, Hchi, chi[n]);

    for (i=0;i<4*VOLUME;++i) {
      int x,y, mu;
      suNf_spinor *chi1, *chi2;
      index_to_coord(i,x,mu);
      _suNf_zero(s1);
      switch (mu) {
      case 0:
	y=iup[x][0];
	chi1=Hchi+x;
	chi2=chi[n]+y;
	_F_DIR0(s1,chi1,chi2);
	chi1=chi[n]+x;
	chi2=Hchi+y;
	_F_DIR0(s1,chi1,chi2);
	break;
      case 1:
	y=iup[x][1];
	chi1=Hchi+x;
	chi2=chi[n]+y;
	_F_DIR1(s1,chi1,chi2);
	chi1=chi[n]+x;
	chi2=Hchi+y;
	_F_DIR1(s1,chi1,chi2);
	break;
      case 2:
	y=iup[x][2];
	chi1=Hchi+x;
	chi2=chi[n]+y;
	_F_DIR2(s1,chi1,chi2);
	chi1=chi[n]+x;
	chi2=Hchi+y;
	_F_DIR2(s1,chi1,chi2);
	break;
      default: /* DIR 3 */
	y=iup[x][3];
	chi1=Hchi+x;
	chi2=chi[n]+y;
	_F_DIR3(s1,chi1,chi2);
	chi1=chi[n]+x;
	chi2=Hchi+y;
	_F_DIR3(s1,chi1,chi2);
      }

      _algebra_project(f,s1);
      /*_print_avect(f); */
      _algebra_vector_mul_add_assign_g(force[i],dt*r_MD.a[n+1]*(_REPR_NORM2/_FUND_NORM2),f);	
      
    }
    

  }  
  
  free(chi[0]);
  free(chi);
  
   

}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3
