/***************************************************************************\
 * Copyright (c) 2014, Claudio Pica                                          *
 * All rights reserved.                                                      * 
 \***************************************************************************/

 /**
  * @file hr_omp.h
  * @brief OpenMP reduction operations for HiRep
  */

#ifndef HR_OMP_H
#define HR_OMP_H

#define _DO_PRAGMA(s) _Pragma ( #s )

#ifdef _OPENMP
#define _OMP_PRAGMA(s) _DO_PRAGMA( omp s )

#include <omp.h>

//define OpenMP behavior
#define _omp_parallel parallel default(shared)
#define _omp_for for schedule(static)
#define _omp_sum(...) reduction(+:__VA_ARGS__)
#define _omp_max(...) reduction(max:__VA_ARGS__)
#define _omp_min(...) reduction(min:__VA_ARGS__)

#define hr_threadId() omp_get_thread_num()

#else //to avoid compilation warnings
#define _OMP_PRAGMA(s)
#define hr_threadId() ((int)(0))
#endif

#define _OMP_BARRIER _OMP_PRAGMA( barrier )

#define _OMP_PARALLEL_FOR \
_OMP_PRAGMA( _omp_parallel) \
_OMP_PRAGMA( _omp_for )

#endif
