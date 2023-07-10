/***************************************************************************\
* Copyright (c) 2022 Sofie Martins                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef GPU_AFFINITY_H
#define GPU_AFFINITY_H

#ifdef WITH_GPU

#if defined(WITH_MPI) && defined(HWLOC)
#include <hwloc.h>
#include <hwloc/cudart.h>
int find_physically_close_CPU_core();
#endif

#endif
#endif
