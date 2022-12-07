#ifndef GEOMETRY_CHECK_H
#define GEOMETRY_CHECK_H
#ifdef __cplusplus
    extern "C" {
#endif

#include "geometry.h"
#include "error.h"

#define _CHECK_GEOMETRY_MATCHING(s1, s2) \
    error((s1)->type!=(s2)->type, 1, __FILE__ ": ", "Geometries don't match!\n");

#define _CHECK_GEOMETRY_EO(s1, s2) \
    do {\
        int pass1 = (s1)->type->glattice && (s2)->type->glattice;\
        int pass2 = (s1)->type->glat_odd && (s2)->type->glat_even;\
        int pass3 = (s1)->type->glat_eve && (s2)->type->glat_odd;\
        error(!(pass1 || pass2 || pass3), 1, __FILE__ ": ", "Incorrect combination of geometries!\n");\
    } while (0)

#ifdef CHECK_SPINOR_MATCHING

#define _TWO_SPINORS_MATCHING(s1,s2) \
        error((s1)->type!=(s2)->type,1,__FILE__ ": ", "Spinors don't match!");

#define _ARRAY_SPINOR_MATCHING(s,n) \
        for(int _i=0; _i<n; _i++) \
                error((s)->type!=((s)+_i)->type,1,__FILE__ ": ", "Spinors don't match!");

#else /* CHECK_SPINOR_MATCHING */

#define _TWO_SPINORS_MATCHING(s1,s2)

#define _ARRAY_SPINOR_MATCHING(s,n)

#endif /* CHECK_SPINOR_MATCHING */

#ifdef __cplusplus
    }
#endif
#endif
