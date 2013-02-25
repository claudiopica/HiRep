/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/


#include "global.h"
#include "linear_algebra.h"
#include "suN.h"
#include "observables.h"
#include "communications.h"

#define _INDEX_(i,s) ( (s)*NF+(i) )

#define _spinor_c_(r,i) (*((suNf_vector*)(&r)+(i)))

#define _spinor_perm_prod_re(k,r,s) \
  { double _ptmp;\
		_vector_prod_re_f(_ptmp,_spinor_c_(r,_C0_),(s).c[0]); (k)=(_S0_)*_ptmp; \
		_vector_prod_re_f(_ptmp,_spinor_c_(r,_C1_),(s).c[1]); (k)+=(_S1_)*_ptmp; \
		_vector_prod_re_f(_ptmp,_spinor_c_(r,_C2_),(s).c[2]); (k)+=(_S2_)*_ptmp; \
		_vector_prod_re_f(_ptmp,_spinor_c_(r,_C3_),(s).c[3]); (k)+=(_S3_)*_ptmp; \
	} while(0)

#define _spinor_perm_prod_im(k,r,s) \
  { double _ptmp;\
		_vector_prod_im_f(_ptmp,_spinor_c_(r,_C0_),(s).c[0]); (k)=(_S0_)*_ptmp; \
		_vector_prod_im_f(_ptmp,_spinor_c_(r,_C1_),(s).c[1]); (k)+=(_S1_)*_ptmp; \
		_vector_prod_im_f(_ptmp,_spinor_c_(r,_C2_),(s).c[2]); (k)+=(_S2_)*_ptmp; \
		_vector_prod_im_f(_ptmp,_spinor_c_(r,_C3_),(s).c[3]); (k)+=(_S3_)*_ptmp; \
	} while(0)



/***************************************************************************\

qp[_INDEX_(i,a)] is supposed to be \psi^{i,a}=(g5 D)^{-1}\eta^{i,a} where
\eta^{i,a} is a point source:
\eta^{i,a}_{j,b}(t,x) = \delta_{ij} \delta_{ab} \delta_{t,t0} \delta_{x,x0}

out[t] = 1/L^3 \sum_x \sum_{a,b,i}
         \psi^{a,i}(x,t0+t)^\dag g5 \bar{\Gamma} \psi^{b,i}(x,t0+t)
         (g5 \Gamma)_{ba}

\***************************************************************************/

#define MESON_DEFINITION(name) \
void name##_correlator(double *out, int t0, spinor_field *qp) { \
  int t,x,y,z, i; \
  suNf_spinor *s1; \
  suNf_spinor *s2; \
  for (t=0; t<GLB_T; t++) out[t] = 0.; \
  for (t=0; t<T; t++) { \
    double _tmp,hc=0.; \
    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { \
      for (i=0; i<NF; ++i) { \
        s1 = _FIELD_AT(&qp[_INDEX_(i,0)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_C0_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_re(_tmp,*s1,*s2); \
        hc += (_S0_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,1)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_C1_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_re(_tmp,*s1,*s2); \
        hc += (_S1_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,2)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_C2_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_re(_tmp,*s1,*s2); \
        hc += (_S2_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,3)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_C3_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_re(_tmp,*s1,*s2); \
        hc += (_S3_)*_tmp; \
      } \
    } \
    out[(zerocoord[0]+t+GLB_T-t0)%GLB_T] = (_SIGN_)*hc/GLB_VOL3; \
  } \
  global_sum(out,GLB_T); \
}


/***************************************************************************\

qp[_INDEX_(i,a)] is supposed to be \psi^{i,a}=(g5 D)^{-1}\eta^{i,a} where
\eta^{i,a} is a point source:
\eta^{i,a}_{j,b}(t,x) = \delta_{ij} \delta_{ab} \delta_{t,t0} \delta_{x,x0}

out[t] = Re{
         1/L^3 \sum_x \sum_{a,b,i}
         \psi^{a,i}(x,t0+t)^\dag g5 \bar{\Gamma}_1 \psi^{b,i}(x,t0+t)
         (g5 \Gamma_2)_{ba}
         }

\***************************************************************************/

#define MESON_DEFINITION_TWO_RE(name) \
void name##_re_correlator(double *out, int t0, spinor_field *qp) { \
  int t,x,y,z, i; \
  suNf_spinor *s1; \
  suNf_spinor *s2; \
  for (t=0; t<GLB_T; t++) out[t] = 0.; \
  for (t=0; t<T; t++) { \
    double _tmp,hc=0.; \
    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { \
      for (i=0; i<NF; ++i) { \
        s1 = _FIELD_AT(&qp[_INDEX_(i,0)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_D0_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_re(_tmp,*s1,*s2); \
        hc += (_T0_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,1)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_D1_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_re(_tmp,*s1,*s2); \
        hc += (_T1_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,2)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_D2_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_re(_tmp,*s1,*s2); \
        hc += (_T2_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,3)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_D3_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_re(_tmp,*s1,*s2); \
        hc += (_T3_)*_tmp; \
      } \
    } \
    out[(zerocoord[0]+t+GLB_T-t0)%GLB_T] = -(_SIGN_)*hc/GLB_VOL3; \
  } \
  global_sum(out,GLB_T); \
}


/***************************************************************************\

qp[_INDEX_(i,a)] is supposed to be \psi^{i,a}=(g5 D)^{-1}\eta^{i,a} where
\eta^{i,a} is a point source:
\eta^{i,a}_{j,b}(t,x) = \delta_{ij} \delta_{ab} \delta_{t,t0} \delta_{x,x0}

out[t] = Im{
         1/L^3 \sum_x \sum_{a,b,i}
         \psi^{a,i}(x,t0+t)^\dag g5 \bar{\Gamma}_1 \psi^{b,i}(x,t0+t)
         (g5 \Gamma_2)_{ba}
         }

\***************************************************************************/

#define MESON_DEFINITION_TWO_IM(name) \
void name##_im_correlator(double *out, int t0, spinor_field *qp) { \
  int t,x,y,z, i; \
  suNf_spinor *s1; \
  suNf_spinor *s2; \
  for (t=0; t<GLB_T; t++) out[t] = 0.; \
  for (t=0; t<T; t++) { \
    double _tmp,hc=0.; \
    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { \
      for (i=0; i<NF; ++i) { \
        s1 = _FIELD_AT(&qp[_INDEX_(i,0)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_D0_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_im(_tmp,*s1,*s2); \
        hc += (_T0_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,1)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_D1_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_im(_tmp,*s1,*s2); \
        hc += (_T1_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,2)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_D2_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_im(_tmp,*s1,*s2); \
        hc += (_T2_)*_tmp; \
        s1 = _FIELD_AT(&qp[_INDEX_(i,3)],ipt(t,x,y,z)); \
        s2 = _FIELD_AT(&qp[_INDEX_(i,_D3_)],ipt(t,x,y,z)); \
				_spinor_perm_prod_im(_tmp,*s1,*s2); \
        hc += (_T3_)*_tmp; \
      } \
    } \
    out[(zerocoord[0]+t+GLB_T-t0)%GLB_T] = -(_SIGN_)*hc/GLB_VOL3; \
  } \
  global_sum(out,GLB_T); \
}


/***************************************************************************\

Build the Gamma matrix (in Gamma[4][4]) from the MACROS, and compute the
sign defined as
g0 Gamma^\dag g0 = sign Gamma

\***************************************************************************/

#define SINGLE_TRACE_DEBUG(name) \
void name##_debug(complex Gamma[4][4], int* sign) { \
  int i,j; \
  for(i=0;i<4;i++) \
  for(j=0;j<4;j++) { \
    Gamma[i][j].re = Gamma[i][j].im = 0.; \
  } \
  if(_REAL_ == 1) { \
    *sign = - _SIGN_; \
    Gamma[0][_C0_].re = _S0_; \
    Gamma[1][_C1_].re = _S1_; \
    Gamma[2][_C2_].re = _S2_; \
    Gamma[3][_C3_].re = _S3_; \
  } else { \
    *sign = _SIGN_; \
    Gamma[0][_C0_].im = _S0_; \
    Gamma[1][_C1_].im = _S1_; \
    Gamma[2][_C2_].im = _S2_; \
    Gamma[3][_C3_].im = _S3_; \
  } \
  for(i=0;i<4;i++) { \
    Gamma[2][i].re = -Gamma[2][i].re; \
    Gamma[2][i].im = -Gamma[2][i].im; \
    Gamma[3][i].re = -Gamma[3][i].re; \
    Gamma[3][i].im = -Gamma[3][i].im; \
  } \
}


/***************************************************************************\

out = Re tr (g5 Gamma smat)

\***************************************************************************/

#define GAMMA_TRACE_RE_DEFINITION(name) \
void name##_trace_H(complex* out, complex* smat) { \
        out->re = _S0_*smat[SPIN_2D_INDEX(_C0_,0)].re \
                + _S1_*smat[SPIN_2D_INDEX(_C1_,1)].re \
                + _S2_*smat[SPIN_2D_INDEX(_C2_,2)].re \
                + _S3_*smat[SPIN_2D_INDEX(_C3_,3)].re; \
        out->im = _S0_*smat[SPIN_2D_INDEX(_C0_,0)].im \
                + _S1_*smat[SPIN_2D_INDEX(_C1_,1)].im \
                + _S2_*smat[SPIN_2D_INDEX(_C2_,2)].im \
                + _S3_*smat[SPIN_2D_INDEX(_C3_,3)].im; \
}


/***************************************************************************\

out = Im tr (g5 Gamma smat)

\***************************************************************************/

#define GAMMA_TRACE_IM_DEFINITION(name) \
void name##_trace_H(complex* out, complex* smat) { \
        out->im = _S0_*smat[SPIN_2D_INDEX(_C0_,0)].re \
                + _S1_*smat[SPIN_2D_INDEX(_C1_,1)].re \
                + _S2_*smat[SPIN_2D_INDEX(_C2_,2)].re \
                + _S3_*smat[SPIN_2D_INDEX(_C3_,3)].re; \
        out->re = - _S0_*smat[SPIN_2D_INDEX(_C0_,0)].im \
                - _S1_*smat[SPIN_2D_INDEX(_C1_,1)].im \
                - _S2_*smat[SPIN_2D_INDEX(_C2_,2)].im \
                - _S3_*smat[SPIN_2D_INDEX(_C3_,3)].im; \
}


/***************************************************************************\

out = g5 Gamma^\dag in

\***************************************************************************/

#define GAMMA_G5GAMMADAG_RE_DEFINITION(name) \
void name##_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in) { \
  int a; \
  for(a=0; a<NF; a++) { \
    out->c[0].c[a].re = _SIGN_DAG_ * _S0_ * in->c[_C0_].c[a].re; \
    out->c[1].c[a].re = _SIGN_DAG_ * _S1_ * in->c[_C1_].c[a].re; \
    out->c[2].c[a].re = _SIGN_DAG_ * _S2_ * in->c[_C2_].c[a].re; \
    out->c[3].c[a].re = _SIGN_DAG_ * _S3_ * in->c[_C3_].c[a].re; \
    out->c[0].c[a].im = _SIGN_DAG_ * _S0_ * in->c[_C0_].c[a].im; \
    out->c[1].c[a].im = _SIGN_DAG_ * _S1_ * in->c[_C1_].c[a].im; \
    out->c[2].c[a].im = _SIGN_DAG_ * _S2_ * in->c[_C2_].c[a].im; \
    out->c[3].c[a].im = _SIGN_DAG_ * _S3_ * in->c[_C3_].c[a].im; \
  } \
}


/***************************************************************************\

out = i g5 Gamma^\dag in

\***************************************************************************/

#define GAMMA_G5GAMMADAG_IM_DEFINITION(name) \
void name##_eval_g5GammaDag_times_spinor(suNf_spinor* out, suNf_spinor* in) { \
  int a; \
  for(a=0; a<NF; a++) { \
    out->c[0].c[a].re = -_SIGN_DAG_ * _S0_ * in->c[_C0_].c[a].im; \
    out->c[1].c[a].re = -_SIGN_DAG_ * _S1_ * in->c[_C1_].c[a].im; \
    out->c[2].c[a].re = -_SIGN_DAG_ * _S2_ * in->c[_C2_].c[a].im; \
    out->c[3].c[a].re = -_SIGN_DAG_ * _S3_ * in->c[_C3_].c[a].im; \
    out->c[0].c[a].im = _SIGN_DAG_ * _S0_ * in->c[_C0_].c[a].re; \
    out->c[1].c[a].im = _SIGN_DAG_ * _S1_ * in->c[_C1_].c[a].re; \
    out->c[2].c[a].im = _SIGN_DAG_ * _S2_ * in->c[_C2_].c[a].re; \
    out->c[3].c[a].im = _SIGN_DAG_ * _S3_ * in->c[_C3_].c[a].re; \
  } \
}




/***************************************************************************\

Single Gamma matrix MACROS

If Gamma has real elements then:
  NAME = Gamma
  _REAL_ = 1
  - g0 Gamma^\dag g0 = _SIGN_ Gamma
  (g5 Gamma)_{ab} = _Sa_ \delta_{_Ca_, b}
  Gamma^\dag = _SIGN_DAG_ \Gamma

If Gamma has imaginary elements then:
  NAME = Gamma
  _REAL_ = 0
  g0 Gamma^\dag g0 = _SIGN_ Gamma
  (g5 Gamma)_{ab} = i _Sa_ \delta_{_Ca_, b}
  Gamma^\dag = _SIGN_DAG_ \Gamma

\***************************************************************************/


#define NAME id

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _S0_ 1
#define _S1_ 1
#define _S2_ -1
#define _S3_ -1

#define _SIGN_ -1
#define _SIGN_DAG_ 1

#define _REAL_ 1

MESON_DEFINITION(id)
GAMMA_TRACE_RE_DEFINITION(id)
SINGLE_TRACE_DEBUG(id)
GAMMA_G5GAMMADAG_RE_DEFINITION(id)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_


#define NAME g0

#define _C0_ 2
#define _C1_ 3
#define _C2_ 0
#define _C3_ 1

#define _S0_ -1
#define _S1_ -1
#define _S2_ 1
#define _S3_ 1

#define _SIGN_ -1
#define _SIGN_DAG_ 1

#define _REAL_ 1

MESON_DEFINITION(g0)
GAMMA_TRACE_RE_DEFINITION(g0)
SINGLE_TRACE_DEBUG(g0)
GAMMA_G5GAMMADAG_RE_DEFINITION(g0)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g5

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _S0_ 1
#define _S1_ 1
#define _S2_ 1
#define _S3_ 1

#define _SIGN_ 1
#define _SIGN_DAG_ 1

#define _REAL_ 1

MESON_DEFINITION(g5)
GAMMA_TRACE_RE_DEFINITION(g5)
SINGLE_TRACE_DEBUG(g5)
GAMMA_G5GAMMADAG_RE_DEFINITION(g5)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g0g5

#define _C0_ 2
#define _C1_ 3
#define _C2_ 0
#define _C3_ 1

#define _S0_ 1
#define _S1_ 1
#define _S2_ 1
#define _S3_ 1

#define _SIGN_ -1
#define _SIGN_DAG_ -1

#define _REAL_ 1

MESON_DEFINITION(g0g5)
GAMMA_TRACE_RE_DEFINITION(g0g5)
SINGLE_TRACE_DEBUG(g0g5)
GAMMA_G5GAMMADAG_RE_DEFINITION(g0g5)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g1

#define _C0_ 3
#define _C1_ 2
#define _C2_ 1
#define _C3_ 0

#define _S0_ -1
#define _S1_ -1
#define _S2_ -1
#define _S3_ -1

#define _SIGN_ -1
#define _SIGN_DAG_ 1

#define _REAL_ 0

MESON_DEFINITION(g1)
GAMMA_TRACE_IM_DEFINITION(g1)
SINGLE_TRACE_DEBUG(g1)
GAMMA_G5GAMMADAG_IM_DEFINITION(g1)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g2

#define _C0_ 3
#define _C1_ 2
#define _C2_ 1
#define _C3_ 0

#define _S0_ -1
#define _S1_ 1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ 1
#define _SIGN_DAG_ 1

#define _REAL_ 1

MESON_DEFINITION(g2)
GAMMA_TRACE_RE_DEFINITION(g2)
SINGLE_TRACE_DEBUG(g2)
GAMMA_G5GAMMADAG_RE_DEFINITION(g2)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g3

#define _C0_ 2
#define _C1_ 3
#define _C2_ 0
#define _C3_ 1

#define _S0_ -1
#define _S1_ 1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ -1
#define _SIGN_DAG_ 1

#define _REAL_ 0

MESON_DEFINITION(g3)
GAMMA_TRACE_IM_DEFINITION(g3)
SINGLE_TRACE_DEBUG(g3)
GAMMA_G5GAMMADAG_IM_DEFINITION(g3)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g0g1

#define _C0_ 1
#define _C1_ 0
#define _C2_ 3
#define _C3_ 2

#define _S0_ -1
#define _S1_ -1
#define _S2_ -1
#define _S3_ -1

#define _SIGN_ 1
#define _SIGN_DAG_ -1

#define _REAL_ 0

MESON_DEFINITION(g0g1)
GAMMA_TRACE_IM_DEFINITION(g0g1)
SINGLE_TRACE_DEBUG(g0g1)
GAMMA_G5GAMMADAG_IM_DEFINITION(g0g1)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g0g2

#define _C0_ 1
#define _C1_ 0
#define _C2_ 3
#define _C3_ 2

#define _S0_ -1
#define _S1_ 1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ -1
#define _SIGN_DAG_ -1

#define _REAL_ 1

MESON_DEFINITION(g0g2)
GAMMA_TRACE_RE_DEFINITION(g0g2)
SINGLE_TRACE_DEBUG(g0g2)
GAMMA_G5GAMMADAG_RE_DEFINITION(g0g2)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g0g3

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _S0_ -1
#define _S1_ 1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ 1
#define _SIGN_DAG_ -1

#define _REAL_ 0

MESON_DEFINITION(g0g3)
GAMMA_TRACE_IM_DEFINITION(g0g3)
SINGLE_TRACE_DEBUG(g0g3)
GAMMA_G5GAMMADAG_IM_DEFINITION(g0g3)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g5g1

#define _C0_ 3
#define _C1_ 2
#define _C2_ 1
#define _C3_ 0

#define _S0_ -1
#define _S1_ -1
#define _S2_ 1
#define _S3_ 1

#define _SIGN_ -1
#define _SIGN_DAG_ -1

#define _REAL_ 0

MESON_DEFINITION(g5g1)
GAMMA_TRACE_IM_DEFINITION(g5g1)
SINGLE_TRACE_DEBUG(g5g1)
GAMMA_G5GAMMADAG_IM_DEFINITION(g5g1)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g5g2

#define _C0_ 3
#define _C1_ 2
#define _C2_ 1
#define _C3_ 0

#define _S0_ -1
#define _S1_ 1
#define _S2_ 1
#define _S3_ -1

#define _SIGN_ 1
#define _SIGN_DAG_ -1

#define _REAL_ 1

MESON_DEFINITION(g5g2)
GAMMA_TRACE_RE_DEFINITION(g5g2)
SINGLE_TRACE_DEBUG(g5g2)
GAMMA_G5GAMMADAG_RE_DEFINITION(g5g2)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g5g3

#define _C0_ 2
#define _C1_ 3
#define _C2_ 0
#define _C3_ 1

#define _S0_ -1
#define _S1_ 1
#define _S2_ 1
#define _S3_ -1

#define _SIGN_ -1
#define _SIGN_DAG_ -1

#define _REAL_ 0

MESON_DEFINITION(g5g3)
GAMMA_TRACE_IM_DEFINITION(g5g3)
SINGLE_TRACE_DEBUG(g5g3)
GAMMA_G5GAMMADAG_IM_DEFINITION(g5g3)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g0g5g1

#define _C0_ 1
#define _C1_ 0
#define _C2_ 3
#define _C3_ 2

#define _S0_ 1
#define _S1_ 1
#define _S2_ -1
#define _S3_ -1

#define _SIGN_ -1
#define _SIGN_DAG_ -1

#define _REAL_ 0

MESON_DEFINITION(g0g5g1)
GAMMA_TRACE_IM_DEFINITION(g0g5g1)
SINGLE_TRACE_DEBUG(g0g5g1)
GAMMA_G5GAMMADAG_IM_DEFINITION(g0g5g1)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g0g5g2

#define _C0_ 1
#define _C1_ 0
#define _C2_ 3
#define _C3_ 2

#define _S0_ 1
#define _S1_ -1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ 1
#define _SIGN_DAG_ -1

#define _REAL_ 1

MESON_DEFINITION(g0g5g2)
GAMMA_TRACE_RE_DEFINITION(g0g5g2)
SINGLE_TRACE_DEBUG(g0g5g2)
GAMMA_G5GAMMADAG_RE_DEFINITION(g0g5g2)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g0g5g3

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _S0_ 1
#define _S1_ -1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ -1
#define _SIGN_DAG_ -1

#define _REAL_ 0

MESON_DEFINITION(g0g5g3)
GAMMA_TRACE_IM_DEFINITION(g0g5g3)
SINGLE_TRACE_DEBUG(g0g5g3)
GAMMA_G5GAMMADAG_IM_DEFINITION(g0g5g3)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _SIGN_
#undef _SIGN_DAG_
#undef _REAL_



#define NAME g5_g0g5_re_correlator

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _D0_ 2
#define _D1_ 3
#define _D2_ 0
#define _D3_ 1

#define _S0_ 1
#define _S1_ 1
#define _S2_ 1
#define _S3_ 1

#define _T0_ 1
#define _T1_ 1
#define _T2_ 1
#define _T3_ 1

#define _SIGN_ -1

MESON_DEFINITION_TWO_RE(g5_g0g5)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _D0_
#undef _D1_
#undef _D2_
#undef _D3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _T0_
#undef _T1_
#undef _T2_
#undef _T3_
#undef _SIGN_



#define NAME g5_g0g5_im_correlator

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _D0_ 2
#define _D1_ 3
#define _D2_ 0
#define _D3_ 1

#define _S0_ 1
#define _S1_ 1
#define _S2_ 1
#define _S3_ 1

#define _T0_ 1
#define _T1_ 1
#define _T2_ 1
#define _T3_ 1

#define _SIGN_ -1

MESON_DEFINITION_TWO_IM(g5_g0g5)

#undef NAME
#undef _C0_
#undef _C1_
#undef _C2_
#undef _C3_
#undef _D0_
#undef _D1_
#undef _D2_
#undef _D3_
#undef _S0_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _T0_
#undef _T1_
#undef _T2_
#undef _T3_
#undef _SIGN_

