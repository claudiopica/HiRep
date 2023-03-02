/**
 * @file geometry_fuse.h
 * @brief FIXME: Add docs
 */

#ifndef GEOMETRY_FUSE_H
#define GEOMETRY_FUSE_H

#include "geometry_omp.h"

#define _FUSE_IDX(type, is) int is = (type)->fuse_mask[_fuse_master_for_ip_##is]

#define _FUSE_MASTER_FOR(type, is) _FUSE_MASTER_FOR_RED(type, is, nowait, )
#define _FUSE_MASTER_FOR_SUM(type, is, ...) _FUSE_MASTER_FOR_RED(type, is, _omp_sum(__VA_ARGS__), )
#define _FUSE_MASTER_FOR_MAX(type, is, ...) _FUSE_MASTER_FOR_RED(type, is, _omp_max(__VA_ARGS__), )
#define _FUSE_MASTER_FOR_MIN(type, is, ...) _FUSE_MASTER_FOR_RED(type, is, _omp_min(__VA_ARGS__), )

#endif
