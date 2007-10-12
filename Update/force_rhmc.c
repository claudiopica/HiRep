#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "dirac.h"
#include "inverters.h"
#include "rational_functions.h"
#include "representation.h"
#include "logger.h"
#include "linear_algebra.h"
#include "memory.h"

#include <malloc.h>
#include <stdio.h>
#include <math.h>

/* declared in update_rhmc.c */
extern rhmc_par _update_par;
extern suNf_spinor **pf;
extern rational_app r_MD; /* used in the action MD evolution */

#define _print_avect(a) printf("(%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e)\n",(a).c1,(a).c2,(a).c3,(a).c4,(a).c5,(a).c6,(a).c7,(a).c8)

#define _print_mat(a) printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.re,(a).c1_2.re,(a).c1_3.re,(a).c2_1.re,(a).c2_2.re,(a).c2_3.re,(a).c3_1.re,(a).c3_2.re,(a).c3_3.re);printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.im,(a).c1_2.im,(a).c1_3.im,(a).c2_1.im,(a).c2_2.im,(a).c2_3.im,(a).c3_1.im,(a).c3_2.im,(a).c3_3.im)

/* we need to compute  Tr  U(x,mu) g_5*(1-g_mu) chi2 # chi1^+
 * where # denotes the tensor product and Tr is the trace on Lorentz space.
 * the strategy is the following:
 * given the form of g_5(1-g_mu) one can compute only the first two lorentz
 * components of the spinor; so we first apply g_5(1-g_mu) to chi2 to find the first
 * two components; then we multiply these two vectors by U(x,mu) and
 * store the result in p.c[0], p.c[1]; when computing the trace we can factorize p.c[0] and p.c[1]
 * as they both multiply two components of chi1^+; we store these factors in p.c[2] and p.c[3].
 * the tensor product is performed by the macro 
 * _suNf_FMAT(u,p): u = p.c[0] # p.c[2]^+ + p.c[1] # p.c[3]^+
 */

/* these macros use the variables ptmp, p */
#define _F_DIR0(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,0)),ptmp);		      \
  _vector_add_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,0)),ptmp);		      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,1)),ptmp);		      \
  _vector_i_add_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,1)),ptmp);		      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_i_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,2)),ptmp);		      \
  _vector_sub_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,2)),ptmp);		      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_add_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,3)),ptmp);		      \
  _vector_i_sub_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,3)),ptmp);		      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_i_add_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)



void Force_rhmc_f(double dt, suNg_algebra_vector *force){
	int i, n, k;
	static suNg_algebra_vector f;
	static suNf_vector ptmp;
	static suNf_spinor p;
	static suNf s1;
	static mshift_par inv_par;
	suNf_spinor **chi;
	suNf_spinor *Hchi;
	double avrforce,maxforce;
	double nsq;
	unsigned int len;

	get_spinor_len(&len);
	/* allocate spinors */
	chi = (suNf_spinor **)malloc(sizeof(*chi)*(r_MD.order));
	chi[0] = alloc_spinor_field_f(r_MD.order+1);
	for (i=1; i<(r_MD.order); ++i) {
		chi[i]=chi[i-1]+len;
	} 
	Hchi = chi[r_MD.order-1]+len;

	/* Compute (H^2-b[n])^-1 * pf */
	/* set up cg parameters */
	inv_par.n = r_MD.order;
	inv_par.shift = r_MD.b;
	inv_par.err2= _update_par.force_prec; /* this should be high for reversibility */
	inv_par.max_iter=0; /* no limit */

	for (k=0; k<_update_par.n_pf; ++k) {
		/* compute inverse vectors chi[i] = (H^2 - b[i])^1 * pf */
		cg_mshift(&inv_par, &H2, pf[k], chi);

		for (n=0; n<r_MD.order; ++n) {

			g5Dphi(_update_par.mass, Hchi, chi[n]);

			lprintf("FORCE_RHMC",50,"[%d] |chi| = %1.8e |Hchi| = %1.8e\n",n,
					sqrt(spinor_field_sqnorm_f(chi[n])),
					sqrt(spinor_field_sqnorm_f(Hchi))
					);

			avrforce=0.;
			maxforce=0.;

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

				_algebra_vector_sqnorm_g(nsq,f);
				avrforce+=sqrt(nsq);
				for(x=0;x<NG*NG-1;++x){
					if(maxforce<fabs(*(((double*)&f)+x))) maxforce=fabs(*(((double*)&f)+x));
				}
			}

			avrforce*=dt*r_MD.a[n+1]*(_REPR_NORM2/_FUND_NORM2)/((double)(4*VOLUME));
			maxforce*=dt*r_MD.a[n+1]*(_REPR_NORM2/_FUND_NORM2);
			lprintf("FORCE_RHMC",50,"[%d] avr |force| = %1.8e maxforce = %1.8e a = %1.8e b = %1.8e\n",n,avrforce,maxforce,r_MD.a[n+1],r_MD.b[n]);
		}  
	}

	free_field(chi[0]);
	free(chi);

}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3
