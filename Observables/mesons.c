#include "global.h"
#include "linear_algebra.h"
#include "suN.h"
#include "observables.h"
#include "inverters.h"
#include <malloc.h>



#define _INDEX_(i,s) ( (s-1)*NF+(i) )

#define _spinor_c_(r,i) (*((suNf_vector*)(&r)+i-1))

#define _spinor_perm_prod_re(r,s) \
  ((_S1_)*_vector_prod_re_f(_spinor_c_(r,_C1_),(s).c1) \
  +(_S2_)*_vector_prod_re_f(_spinor_c_(r,_C2_),(s).c2) \
  +(_S3_)*_vector_prod_re_f(_spinor_c_(r,_C3_),(s).c3) \
  +(_S4_)*_vector_prod_re_f(_spinor_c_(r,_C4_),(s).c4))


#define MESON_DEFINITION \
void NAME(float *out, suNf_spinor **qp) { \
  int t,x,y,z, i; \
  suNf_spinor *s1; \
  suNf_spinor *s2; \
  for (t=0; t<T; ++t) { \
    double hc=0.; \
    for (x=0; x<L; ++x) for (y=0; y<L; ++y) for (z=0; z<L; ++z) { \
      for (i=0; i<NF; ++i) { \
        s1 = &(qp[_INDEX_(i,1)][ipt[t][x][y][z]]); \
        s2 = &(qp[_INDEX_(i,_C1_)][ipt[t][x][y][z]]); \
        hc += (_S1_)*_spinor_perm_prod_re(*s1,*s2); \
        s1 = &(qp[_INDEX_(i,2)][ipt[t][x][y][z]]); \
        s2 = &(qp[_INDEX_(i,_C2_)][ipt[t][x][y][z]]); \
        hc += (_S2_)*_spinor_perm_prod_re(*s1,*s2); \
        s1 = &(qp[_INDEX_(i,3)][ipt[t][x][y][z]]); \
        s2 = &(qp[_INDEX_(i,_C3_)][ipt[t][x][y][z]]); \
        hc += (_S3_)*_spinor_perm_prod_re(*s1,*s2); \
        s1 = &(qp[_INDEX_(i,4)][ipt[t][x][y][z]]); \
        s2 = &(qp[_INDEX_(i,_C4_)][ipt[t][x][y][z]]); \
        hc += (_S4_)*_spinor_perm_prod_re(*s1,*s2); \
      } \
    } \
    out[t] = (_S0_)*hc; \
  } \
}



#define NAME id_correlator

#define _C1_ 1
#define _C2_ 2
#define _C3_ 3
#define _C4_ 4

#define _S1_ 1
#define _S2_ 1
#define _S3_ -1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g0_correlator

#define _C1_ 3
#define _C2_ 4
#define _C3_ 1
#define _C4_ 2

#define _S1_ 1
#define _S2_ 1
#define _S3_ -1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g5_correlator

#define _C1_ 1
#define _C2_ 2
#define _C3_ 3
#define _C4_ 4

#define _S1_ 1
#define _S2_ 1
#define _S3_ 1
#define _S4_ 1

#define _S0_ 1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g0g5_correlator

#define _C1_ 3
#define _C2_ 4
#define _C3_ 1
#define _C4_ 2

#define _S1_ 1
#define _S2_ 1
#define _S3_ 1
#define _S4_ 1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g1_correlator

#define _C1_ 4
#define _C2_ 3
#define _C3_ 2
#define _C4_ 1

#define _S1_ -1
#define _S2_ -1
#define _S3_ -1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g2_correlator

#define _C1_ 4
#define _C2_ 3
#define _C3_ 2
#define _C4_ 1

#define _S1_ 1
#define _S2_ -1
#define _S3_ 1
#define _S4_ -1

#define _S0_ 1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g3_correlator

#define _C1_ 3
#define _C2_ 4
#define _C3_ 1
#define _C4_ 2

#define _S1_ 1
#define _S2_ -1
#define _S3_ 1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g0g1_correlator

#define _C1_ 2
#define _C2_ 1
#define _C3_ 4
#define _C4_ 3

#define _S1_ -1
#define _S2_ -1
#define _S3_ -1
#define _S4_ -1

#define _S0_ 1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g0g2_correlator

#define _C1_ 2
#define _C2_ 1
#define _C3_ 4
#define _C4_ 3

#define _S1_ 1
#define _S2_ -1
#define _S3_ 1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g0g3_correlator

#define _C1_ 1
#define _C2_ 2
#define _C3_ 3
#define _C4_ 4

#define _S1_ 1
#define _S2_ -1
#define _S3_ 1
#define _S4_ -1

#define _S0_ 1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g5g1_correlator

#define _C1_ 4
#define _C2_ 3
#define _C3_ 2
#define _C4_ 1

#define _S1_ 1
#define _S2_ 1
#define _S3_ -1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g5g2_correlator

#define _C1_ 4
#define _C2_ 3
#define _C3_ 2
#define _C4_ 1

#define _S1_ -1
#define _S2_ 1
#define _S3_ 1
#define _S4_ -1

#define _S0_ 1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g5g3_correlator

#define _C1_ 3
#define _C2_ 4
#define _C3_ 1
#define _C4_ 2

#define _S1_ -1
#define _S2_ 1
#define _S3_ 1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g0g5g1_correlator

#define _C1_ 2
#define _C2_ 1
#define _C3_ 4
#define _C4_ 3

#define _S1_ 1
#define _S2_ 1
#define _S3_ -1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g0g5g2_correlator

#define _C1_ 2
#define _C2_ 1
#define _C3_ 4
#define _C4_ 3

#define _S1_ -1
#define _S2_ 1
#define _S3_ 1
#define _S4_ -1

#define _S0_ 1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



#define NAME g0g5g3_correlator

#define _C1_ 1
#define _C2_ 2
#define _C3_ 3
#define _C4_ 4

#define _S1_ -1
#define _S2_ 1
#define _S3_ 1
#define _S4_ -1

#define _S0_ -1

MESON_DEFINITION

#undef NAME
#undef _C1_
#undef _C2_
#undef _C3_
#undef _C4_
#undef _S1_
#undef _S2_
#undef _S3_
#undef _S4_
#undef _S0_



