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

#ifdef WITH_GPU
   // GPU Dirac operator implementations
   #include "dirac_gpu.h" 
#endif

#endif
