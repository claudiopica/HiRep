/***************************************************************************\
 * Copyright (c) 2014, Claudio Pica                                          *
 * All rights reserved.                                                      * 
 \***************************************************************************/

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



#else //to avoid compilation warnings
#define _OMP_PRAGMA(s)
#endif
#define _OMP_BARRIER _OMP_PRAGMA( barrier )


#endif
