/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file dirac.h
 * @brief Implementation of the Dirac operator
 */

#ifndef DIRAC_H
#define DIRAC_H

#include "suN_types.h"
#include "utils.h"

// For cross-compilation we need to mark CPU functions as C code. We only cross-compile
// C&C++ for the CUDA version.
#ifdef __cplusplus
extern "C" {
#endif

   // Default declarations for Dirac Operator
   #include "dirac_default.h"

   // CPU declarations
   #include "dirac_cpu.h"

   // Clover Operations
   #include "dirac_clover.h"

   // Four-fermion interactions
   #include "dirac_4f.h"

   // Twisted-Mass
   #include "dirac_tm.h"

#ifdef __cplusplus
}
#endif

// GPU declarations, will be empty if not compiled WITH_GPU
#include "dirac_gpu.h" 

#endif
