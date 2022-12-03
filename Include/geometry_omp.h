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
#define _PIECE_FOR(type, ip) \
  for (int ip = 0; ip < (type)->local_master_pieces; ip++)

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
  _PIECE_FOR((type), _master_for_ip_##is)         \
  _SITE_FOR_RED((type), _master_for_ip_##is, is, redop1, redop2)

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
  _OMP_PRAGMA(_omp_parallel)                        \
  _OMP_PRAGMA(_omp_for redop1 redop2)               \
  for (int ip = 0; ip < type->fuse_gauge_size; ip++)

/**
 * @brief Fuse reduce on the whole local lattice FIXME: more desc
 *
 *  * @param type          geometry_descriptor that contains information on the geometry
 *                      of the local lattice.
 * @param is            Local variable that runs over all site indices on the given piece.
 * @param redop1        Variable to reduce
 * @param redop2        Variable to reduce 
 */
#define _FUSE_MASTER_FOR_RED(type, is, redop1, redop2) \
  _FUSE_FOR_RED((type), _fuse_master_for_ip_##is, is, redop1, redop2)

/**
 * @brief Have the spinor to iterate
 */
#define _ONE_SPINOR_FOR_RED(s,redop1,redop2) _MASTER_FOR_RED((s)->type,_spinor_for_is,redop1,redop2)

#define _ONE_SPINOR_FOR_SUM(s,...) _ONE_SPINOR_FOR_RED(s,_omp_sum(__VA_ARGS__),)

#define _TWO_SPINORS_FOR_RED(s1,s2,redop1,redop2) \
  _TWO_SPINORS_MATCHING(s1,s2); \
  _ONE_SPINOR_FOR_RED(s1,redop1,redop2)

#define _TWO_SPINORS_FOR_SUM(s1,s2,...) _TWO_SPINORS_FOR_RED(s1,s2,_omp_sum(__VA_ARGS__),)


#define _THREE_SPINORS_FOR_RED(s1,s2,s3,redop1,redop2) \
  _TWO_SPINORS_MATCHING(s1,s2); \
  _TWO_SPINORS_MATCHING(s1,s3); \
  _ONE_SPINOR_FOR_RED(s1,redop1,redop2)
#endif
