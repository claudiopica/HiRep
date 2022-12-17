/**
 * @file geometry_check.h
 * @brief Validity checks on the geometry of spinor fields that are passed as 
 *        function parameters. Use these macros to make the code more compact. These
 *        are compiled to empty lines, if CHECK_SPINOR_MATCHING is inactive.
 */

#ifndef GEOMETRY_CHECK_H
#define GEOMETRY_CHECK_H
#ifdef __cplusplus
    extern "C" {
#endif

// TODO: This needs adjustment to the new geometry (SAM)

#include "geometry.h"
#include "error.h"

#ifdef CHECK_SPINOR_MATCHING

    /**
     * @brief Check whether two spinors are defined on the same subset of the lattice.
     *        Usually this means to check whether the fields either both have odd parity, 
     *        even parity or are both defined on the full lattice.
     *
     * @param s1                        First input field
     * @param s2                        Second input field against which to check the first.
     */
    #define _TWO_SPINORS_MATCHING(s1,s2) \
            error((s1)->type != (s2)->type, 1, __FILE__ ": ", \
                    "Spinor geometries don't agree!");

    /**
     * @brief Check whether all elements of a spinor field array are allocated with the
     *        same parity (most of the time either glat_even, glat_odd or glattice).
     *
     * @param s                         Spinor field array to check
     * @param n                         Number of fields in the array.
     */
    #define _ARRAY_SPINOR_MATCHING(s,n) \
            for(int _i=0; _i<n; _i++) \
                    error((s)->type != ((s)+_i)->type, 1, __FILE__ ": ", \
                            "Spinors geometries don;t agree!");

    /**
     * @brief Check whether two fields are defined on the same subset of the lattice/with
     *        matching parity. Usually this means to check whether the fields either both
     *        have odd parity, even parity or are both defined on the full lattice.
     *
     * @param s1                First input field
     * @param s2                Second input field against which to check the first.
     */
    #define _CHECK_GEOMETRY_MATCHING(s1, s2) \
        error((s1)->type!=(s2)->type, 1, __FILE__ ": ", 
                    "Field geometries don't agree!\n");

    /**
     * @brief Check whether two fields are defined with even-odd parity, meaning if one
     *        field is odd the other one is even or the other way around. Allowed is also
     *        the case where both are allocated on the full lattice.
     *
     * @param s1                First input field
     * @param s2                Second input field against which to check the first.
     */
    #define _CHECK_GEOMETRY_EO(s1, s2) \
        do {\
            int pass1 = (s1)->type->glattice && (s2)->type->glattice;\
            int pass2 = (s1)->type->glat_odd && (s2)->type->glat_even;\
            int pass3 = (s1)->type->glat_even && (s2)->type->glat_odd;\
            error(!(pass1 || pass2 || pass3), 1, __FILE__ ": ", \
                    "Incorrect combination of geometries! " \
                    "Need even to odd, odd to even or both defined on the full lattice\n");\
        } while (0)

    #define _CHECK_GEOMETRY_FULL(s1) \
        do { \
            error(s1->type!=&glattice, 1, "Dphi_gpu [Dphi_gpu.c]", \
            "Spinor is not defined on all the lattice!"); \
        } while (0)

    #define _CHECK_GEOMETRY_EVEN(s1) \
        do { \
            error(s1->type!=&glat_even, 1, "Dphi_gpu [Dphi_gpu.c]", \
            "Spinor needs to be even!"); \
        } while (0)

    #define _CHECK_GEOMETRY_ODD(s1) \
        do { \
            error(s1->type!=&glat_odd, 1, "Dphi_gpu [Dphi_gpu.c]", \
            "Spinor needs to be odd!"); \
        } while (0)

#else 

    /**
     * @brief This macro does nothing unless CHECK_SPINOR_MATCHING is defined.
     */
    #define _TWO_SPINORS_MATCHING(s1,s2) 

    /**
     * @brief This macro does nothing unless CHECK_SPINOR_MATCHING is defined.
     */
    #define _ARRAY_SPINOR_MATCHING(s,n) 

    /**
     * @brief This macro does nothing unless CHECK_SPINOR_MATCHING is defined.
     */
    #define _CHECK_GEOMETRY_MATCHING(s1, s2)

    /**
     * @brief This macro does nothing unless CHECK_SPINOR_MATCHING is defined.
     */
    #define _CHECK_GEOMETRY_EO(s1, s2)

    #define _CHECK_GEOMETRY_FULL(s1)

    #define _CHECK_GEOMETRY_EVEN(s1)

    #define _CHECK_GEOMETRY_ODD(s1)

#endif

#ifdef __cplusplus
    }
#endif
#endif
