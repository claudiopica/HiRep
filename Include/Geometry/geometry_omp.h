/**
 * @file geometry_omp.h
 * @brief This file contains useful macros that perform OpenMP reduction operations
 *        and are necessary to define iterations over sites on the lattice.
 */

#ifndef GEOMETRY_OMP_H
#define GEOMETRY_OMP_H

#include "hr_omp.h"

/**
 * @brief Iterate over the pieces of the given geometry.
 *
 * @param type          geometry_descriptor that contains information on the geometry
 *                      of the local and global lattice
 * @param ip            local variable that denotes the piece index
 */
#define _PIECE_FOR(type, ip) for (int ip = 0; ip < (type)->local_master_pieces; ip++)

/**
 * @brief Reduced iteration over sites of a given piece. Variables given as redop parameters
 * 	  are reduced using an OpenMP reduction sum.
 *
 * @param type		geometry_descriptor that contains information on the geometry
 * 			of the local and global lattice
 * @param ip		Identifying index of the current piece
 * @param is		Local variable that runs over all site indices on the piece
 * @param redop1        Variable to reduce
 * @param redop2	Variable to reduce
 *
 */
#define _SITE_FOR_RED(type, ip, is, redop1, redop2) \
    _OMP_PRAGMA(_omp_parallel)                      \
    _OMP_PRAGMA(_omp_for redop1 redop2)             \
    for (size_t is = (type)->master_start[ip]; is <= (type)->master_end[ip]; is++)

/**
 * @brief Iterate over sites of a given piece.
 *
 * @param type		geometry_descriptor that contains information on the geometry of
 * 			the local lattice
 * @param ip		Identifying index of the current piece
 * @param is		Local variable that runs over all site indices on the given piece
 */
#define _SITE_FOR(type, ip, is) _SITE_FOR_RED(type, ip, is, nowait, )

/**
 * @brief Iterate over sites on a given piece and omp-reduce using a sum.
 *
 * @param type		geometry_descriptor that contains information on the geometry
 * 			of the local lattice.
 * @param ip		Identifying index of the current piece
 * @param is		local variable that runs over all site indices on the given piece
 */
#define _SITE_FOR_SUM(type, ip, is, ...) _SITE_FOR_RED(type, ip, is, _omp_sum(__VA_ARGS__), )

/**
 * @brief Iterate over sites on a given piece and omp-reduce using a max.
 * 
 * @param type          geometry_descriptor that contains information on the geometry
 *                      of the local lattice.
 * @param ip            Identifying index of the current piece
 * @param is            local variable that runs over all site indices on the given piece
 */
#define _SITE_FOR_MAX(type, ip, is, ...) _SITE_FOR_RED(type, ip, is, _omp_max(__VA_ARGS__), )

/**
 * @brief Iterate over sites on a given piece and omp-reduce using a min.
 *  *
 * @param type          geometry_descriptor that contains information on the geometry
 *                      of the local lattice.
 * @param ip            Identifying index of the current piece
 * @param is            local variable that runs over all site indices on the given piece
 */
#define _SITE_FOR_MIN(type, ip, is, ...) _SITE_FOR_RED(type, ip, is, _omp_min(__VA_ARGS__), )

/**
 * @brief Reduced iteration over all sites of the local lattice. Variables given as redop
 *        parameters will be reduced using an OpenMP reduction sum.
 *
 * @param type          geometry_descriptor that contains information on the geometry
 *                      of the local and global lattice
 * @param is            Local variable that runs over all site indices on the local lattice
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce
 */
#define _MASTER_FOR_RED(type, is, redop1, redop2) \
    _PIECE_FOR((type), _master_for_ip_##is)       \
        _SITE_FOR_RED((type), _master_for_ip_##is, is, redop1, redop2)

/**
 * @brief Iterate over all sites of the local lattice
 *
 * @param type		geometry_descriptor that contains information on the geometry
 * 			of the local lattice
 * @param is		Local variable that runs over all site indices
 */
#define _MASTER_FOR(type, is) _MASTER_FOR_RED(type, is, , )

/**
 * @brief Reduced iteration over all sites of the local lattice and omp-reduce using a sum.
 *
 * @param type 		geometry_descriptor that contains information on the geometry
 * 			of the local lattice.
 * @param is		Local variable that runs over all site indices on the given piece
 */
#define _MASTER_FOR_SUM(type, is, ...) _MASTER_FOR_RED(type, is, _omp_sum(__VA_ARGS__), )

/**
 * @brief Reduced iteration over all sites of the local lattice and omp-reduce using a max.
 * 
 * @param type          geometry_descriptor that contains information on the geometry
 *                      of the local lattice.
 * @param is            Local variable that runs over all site indices on the given piece
 */
#define _MASTER_FOR_MAX(type, is, ...) _MASTER_FOR_RED(type, is, _omp_max(__VA_ARGS__), )

/**
 * @brief Reduced iteration over all sites of the local lattice and omp-reduce using a min.
 * 
 * @param type          geometry_descriptor that contains information on the geometry
 *                      of the local lattice.
 * @param is            Local variable that runs over all site indices on the given piece
 */
#define _MASTER_FOR_MIN(type, is, ...) _MASTER_FOR_RED(type, is, _omp_min(__VA_ARGS__), )

/**
 * @brief Fuse reduced on the current piece FIXME: more desc
 *
 * @param type		geometry_descriptor that contains information on the geometry
 * 			of the local lattice.
 * @param ip 		Identifying index of the current piece
 * @param is		Local variable that runs over all site indices on the given piece.
 * @param redop1	Variable to reduce
 * @param redop2	Variable to reduce		
 */
#define _FUSE_FOR_RED(type, ip, is, redop1, redop2) \
    _OMP_PRAGMA(_omp_parallel)                      \
    _OMP_PRAGMA(_omp_for redop1 redop2)             \
    for (int ip = 0; ip < type->fuse_gauge_size; ip++)

/**
 * @brief Fuse reduce on the whole local lattice TODO: more desc
 *
 * @param type          geometry_descriptor that contains information on the geometry
 *                      of the local lattice.
 * @param is            Local variable that runs over all site indices on the given piece.
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce 
 */
#define _FUSE_MASTER_FOR_RED(type, is, redop1, redop2) _FUSE_FOR_RED((type), _fuse_master_for_ip_##is, is, redop1, redop2)

/**
 * @brief Iterate over all sites of the local lattice but not by index in memory
 *        but by spinor, applying a OpenMP reduction to the other variables.
 *        The current spinor can be found using _SPINOR_PTR(s). 
 *
 * @param s             Input spinor field
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce
 */
#define _ONE_SPINOR_FOR_RED(s, redop1, redop2) _MASTER_FOR_RED((s)->type, _spinor_for_is, redop1, redop2)

/**
 * @brief Iterate over all sites of the local lattice but not by index in memory
 *        but by site, applying a OpenMP reduction to the other variables.
 *        The current site can be found using _SITE_PTR(s). 
 *
 * @param s             Input field
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce
 */
#define _ONE_SITE_FOR_RED(s, redop1, redop2) _MASTER_FOR_RED((s)->type, _site_for_is, redop1, redop2)

/**
 * @brief Iterate over all sites of the local lattice but not by index in memory
 *        but by spinor, applying an OpenMP sum reduction to the other variables.
 *        The current spinor can be found using _SPINOR_PTR(s).
 *
 * @param s             Input spinor field
 * @param ...           Variables to reduce
 */
#define _ONE_SPINOR_FOR_SUM(s, ...) _ONE_SPINOR_FOR_RED(s, _omp_sum(__VA_ARGS__), )

/**
 * @brief Iterate over all sites of the local lattice but not by index in memory
 *        but by site, applying an OpenMP sum reduction to the other variables.
 *        The current site can be found using _SITE_PTR(s).
 *
 * @param s             Input field
 * @param ...           Variables to reduce
 */
#define _ONE_SITE_FOR_SUM(s, ...) _ONE_SITE_FOR_RED(s, _omp_sum(__VA_ARGS__), )

/**
 * @brief Iterate over two corresponding spinors on the given fields, applying
 *        an OpenMP reduction operation on the other given variables. The current 
 *        spinors can be found using _SPINOR_PTR(s1) and _SPINOR_PTR(s2).
 *
 * @param s1            First input spinor field
 * @param s2            Second input spinor field
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce
 */
#define _TWO_SPINORS_FOR_RED(s1, s2, redop1, redop2) \
    _TWO_SPINORS_MATCHING(s1, s2);                   \
    _ONE_SPINOR_FOR_RED(s1, redop1, redop2)

/**
 * @brief Iterate over two corresponding sites on the given fields, applying
 *        an OpenMP reduction operation on the other given variables. The current 
 *        sites can be found using _SITE_PTR(s1) and _SITE_PTR(s2).
 *
 * @param s1            First input field
 * @param s2            Second input field
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce
 */
#define _TWO_SITE_FOR_RED(s1, s2, redop1, redop2) \
    _CHECK_GEOMETRY_MATCHING(s1, s2);             \
    _ONE_SITE_FOR_RED(s1, redop1, redop2)

/**
 * @brief Iterate over two corresponding spinors on the given fields, applying
 *        an OpenMP sum reduction on the other given variables. The current
 *        spinors can be found using _SPINOR_PTR(s1) and _SPINOR_PTR(s2).
 *
 * @param s1            First input spinor field
 * @param s2            Second input spinor field
 * @param ...           Variables to reduce
 */
#define _TWO_SPINORS_FOR_SUM(s1, s2, ...) _TWO_SPINORS_FOR_RED(s1, s2, _omp_sum(__VA_ARGS__), )

/**
 * @brief Iterate over two corresponding sites on the given fields, applying
 *        an OpenMP sum reduction on the other given variables. The current
 *        sites can be found using _SITE_PTR(s1) and _SITE_PTR(s2).
 *
 * @param s1            First input field
 * @param s2            Second input field
 * @param ...           Variables to reduce
 */
#define _TWO_SITE_FOR_SUM(s1, s2, ...) _TWO_SITE_FOR_RED(s1, s2, _omp_sum(__VA_ARGS__), )

/**
 * @brief Iterate over three corresponding spinors on the given fields, applying
 *        an OpenMP reduction on the other given variables. The current
 *        spinors can be found using _SPINOR_PTR(s1) and _SPINOR_PTR(s2).
 * 
 * @param s1            First input spinor field
 * @param s2            Second input spinor field
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce
 */
#define _THREE_SPINORS_FOR_RED(s1, s2, s3, redop1, redop2) \
    _TWO_SPINORS_MATCHING(s1, s2);                         \
    _TWO_SPINORS_MATCHING(s1, s3);                         \
    _ONE_SPINOR_FOR_RED(s1, redop1, redop2)

/**
 * @brief Iterate over three corresponding sites on the given fields, applying
 *        an OpenMP reduction on the other given variables. The current
 *        sites can be found using _SITE_PTR(s1) and _SITE_PTR(s2).
 * 
 * @param s1            First input field
 * @param s2            Second input field
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce
 */
#define _THREE_SITE_FOR_RED(s1, s2, s3, redop1, redop2) \
    _CHECK_GEOMETRY_MATCHING(s1, s2);                   \
    _CHECK_GEOMETRY_MATCHING(s1, s3);                   \
    _ONE_SITE_FOR_RED(s1, redop1, redop2)

/**
 * @brief Iterate over all sites of the local lattice but not by index in memory but by
 *        spinor. The current spinor can be found using _SPINOR_PTR(s). 
 *
 * @param s		Input spinor field
 */
#define _ONE_SPINOR_FOR(s) _ONE_SPINOR_FOR_RED(s, , )

/**
 * @brief Iterate over all sites of the local lattice but not by index in memory but by
 *        site. The current site can be found using _SITE_PTR(s). 
 *
 * @param s		Input field
 */
#define _ONE_SITE_FOR(s) _ONE_SITE_FOR_RED(s, , )

/**
 * @brief Iterate over two corresponding spinors on the given fields. The current spinors 
 * 	  can be found using _SPINOR_PTR(s1) and _SPINOR_PTR(s2).
 *
 * @param s1		First input spinor field
 * @param s2		Second input spinor field
 */
#define _TWO_SPINORS_FOR(s1, s2) _TWO_SPINORS_FOR_RED(s1, s2, , )

/**
 * @brief Iterate over two corresponding sites on the given fields. The current sites 
 * 	  can be found using _SITE_PTR(s1) and _SITE_PTR(s2).
 *
 * @param s1		First input spinor field
 * @param s2		Second input spinor field
 */
#define _TWO_SITE_FOR(s1, s2) _TWO_SITE_FOR_RED(s1, s2, , )

/**
 * @brief Iterate over all three corresponding spinors on the given fields. The current 
 * 	  spinors can be found using _SPINOR_PTR(s1), _SPINOR_PTR(s2), _SPINOR_PTR(s3).
 *
 * @param s1		First input spinor field
 * @param s2		Second input spinor field
 * @param s3		Third input spinor field
 */
#define _THREE_SPINORS_FOR(s1, s2, s3) _THREE_SPINORS_FOR_RED(s1, s2, s3, , )

/**
 * @brief Iterate over all three corresponding sites on the given fields. The current 
 * 	  sites can be found using _SITE_PTR(s1), _SITE_PTR(s2), _SITE_PTR(s3).
 *
 * @param s1		First input spinor field
 * @param s2		Second input spinor field
 * @param s3		Third input spinor field
 */
#define _THREE_SITE_FOR(s1, s2, s3) _THREE_SITE_FOR_RED(s1, s2, s3, , )

/**
 * @brief Retrieve current spinor. This macro only works inside _SPINOR_FOR, 
 * 	  _TWO_SPINORS_FOR or _THREE_SPINORS_FOR.
 *
 * @param s 		Spinor field that is being iterated over.
 */
#define _SPINOR_PTR(s) _FIELD_AT(s, _spinor_for_is)

/**
 * @brief Retrieve current site. This macro only works inside _SITE_FOR, 
 * 	  _TWO_SITE_FOR or _THREE_SITE_FOR.
 *
 * @param s 		Field that is being iterated over.
 */
#define _SITE_PTR(__s, __mu, __dim) (_DFIELD_AT(__s, _site_for_is, __mu, __dim))

#endif
