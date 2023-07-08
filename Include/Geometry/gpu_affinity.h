/***************************************************************************\
* Copyright (c) 2022 Sofie Martins                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef GPU_AFFINITY_H
#define GPU_AFFINITY_H

#if defined(WITH_MPI) && defined(HWLOC)
int find_physically_close_CPU_core();
#endif

#endif