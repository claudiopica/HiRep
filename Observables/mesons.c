/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
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


#define MESON_DEFINITION \
void NAME(double *out, spinor_field *qp) { \
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
    out[COORD[0]*T+t] = (_SIGN_)*hc/(GLB_X*GLB_Y*GLB_Z); \
  } \
  global_sum(out,GLB_T); \
}


#define MESON_DEFINITION_TWO_RE \
void NAME(double *out, spinor_field *qp) { \
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
    out[COORD[0]*T+t] = -(_SIGN_)*hc/(GLB_X*GLB_Y*GLB_Z); \
  } \
  global_sum(out,GLB_T); \
}


#define MESON_DEFINITION_TWO_IM \
void NAME(double *out, spinor_field *qp) { \
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
    out[COORD[0]*T+t] = -(_SIGN_)*hc/(GLB_X*GLB_Y*GLB_Z); \
  } \
  global_sum(out,GLB_T); \
}


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



#define NAME id_correlator

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _S0_ 1
#define _S1_ 1
#define _S2_ -1
#define _S3_ -1

#define _SIGN_ -1

#define _REAL_ 1

MESON_DEFINITION
SINGLE_TRACE_DEBUG(id)

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
#undef _REAL_


#define NAME g0_correlator

#define _C0_ 2
#define _C1_ 3
#define _C2_ 0
#define _C3_ 1

#define _S0_ -1
#define _S1_ -1
#define _S2_ 1
#define _S3_ 1

#define _SIGN_ -1

#define _REAL_ 1

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g0)

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
#undef _REAL_



#define NAME g5_correlator

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _S0_ 1
#define _S1_ 1
#define _S2_ 1
#define _S3_ 1

#define _SIGN_ 1

#define _REAL_ 1

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g5)

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
#undef _REAL_



#define NAME g0g5_correlator

#define _C0_ 2
#define _C1_ 3
#define _C2_ 0
#define _C3_ 1

#define _S0_ 1
#define _S1_ 1
#define _S2_ 1
#define _S3_ 1

#define _SIGN_ -1

#define _REAL_ 1

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g0g5)

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
#undef _REAL_



#define NAME g1_correlator

#define _C0_ 3
#define _C1_ 2
#define _C2_ 1
#define _C3_ 0

#define _S0_ -1
#define _S1_ -1
#define _S2_ -1
#define _S3_ -1

#define _SIGN_ -1

#define _REAL_ 0

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g1)

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
#undef _REAL_



#define NAME g2_correlator

#define _C0_ 3
#define _C1_ 2
#define _C2_ 1
#define _C3_ 0

#define _S0_ -1
#define _S1_ 1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ 1

#define _REAL_ 1

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g2)

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
#undef _REAL_



#define NAME g3_correlator

#define _C0_ 2
#define _C1_ 3
#define _C2_ 0
#define _C3_ 1

#define _S0_ -1
#define _S1_ 1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ -1

#define _REAL_ 0

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g3)

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
#undef _REAL_



#define NAME g0g1_correlator

#define _C0_ 1
#define _C1_ 0
#define _C2_ 3
#define _C3_ 2

#define _S0_ -1
#define _S1_ -1
#define _S2_ -1
#define _S3_ -1

#define _SIGN_ 1

#define _REAL_ 0

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g0g1)

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
#undef _REAL_



#define NAME g0g2_correlator

#define _C0_ 1
#define _C1_ 0
#define _C2_ 3
#define _C3_ 2

#define _S0_ -1
#define _S1_ 1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ -1

#define _REAL_ 1

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g0g2)

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
#undef _REAL_



#define NAME g0g3_correlator

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _S0_ -1
#define _S1_ 1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ 1

#define _REAL_ 0

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g0g3)

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
#undef _REAL_



#define NAME g5g1_correlator

#define _C0_ 3
#define _C1_ 2
#define _C2_ 1
#define _C3_ 0

#define _S0_ -1
#define _S1_ -1
#define _S2_ 1
#define _S3_ 1

#define _SIGN_ -1

#define _REAL_ 0

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g5g1)

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
#undef _REAL_



#define NAME g5g2_correlator

#define _C0_ 3
#define _C1_ 2
#define _C2_ 1
#define _C3_ 0

#define _S0_ -1
#define _S1_ 1
#define _S2_ 1
#define _S3_ -1

#define _SIGN_ 1

#define _REAL_ 1

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g5g2)

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
#undef _REAL_



#define NAME g5g3_correlator

#define _C0_ 2
#define _C1_ 3
#define _C2_ 0
#define _C3_ 1

#define _S0_ -1
#define _S1_ 1
#define _S2_ 1
#define _S3_ -1

#define _SIGN_ -1

#define _REAL_ 0

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g5g3)

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
#undef _REAL_



#define NAME g0g5g1_correlator

#define _C0_ 1
#define _C1_ 0
#define _C2_ 3
#define _C3_ 2

#define _S0_ 1
#define _S1_ 1
#define _S2_ -1
#define _S3_ -1

#define _SIGN_ -1

#define _REAL_ 0

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g0g5g1)

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
#undef _REAL_



#define NAME g0g5g2_correlator

#define _C0_ 1
#define _C1_ 0
#define _C2_ 3
#define _C3_ 2

#define _S0_ 1
#define _S1_ -1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ 1

#define _REAL_ 1

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g0g5g2)

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
#undef _REAL_



#define NAME g0g5g3_correlator

#define _C0_ 0
#define _C1_ 1
#define _C2_ 2
#define _C3_ 3

#define _S0_ 1
#define _S1_ -1
#define _S2_ -1
#define _S3_ 1

#define _SIGN_ -1

#define _REAL_ 0

MESON_DEFINITION
SINGLE_TRACE_DEBUG(g0g5g3)

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

MESON_DEFINITION_TWO_RE

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

MESON_DEFINITION_TWO_IM

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

