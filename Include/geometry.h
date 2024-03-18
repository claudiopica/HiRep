/***************************************************************************\
 * Copyright (c) 2008-2022, Claudio Pica, Sofie Martins                     *
 * All rights reserved.                                                     *
 \**************************************************************************/

/**
 * @file geometry.h
 * @brief This file contains information on the geometry of the local 
 * 	  lattice, block decomposed geometry, buffers, etc. All geometric
 * 	  properties that are needed to perform operations correctly
 * 	  have to be contained in the struct geometry_descriptor.
 * 	  In order to encapsulate the complicated geometry structure
 * 	  iteration over lattice sites is simplified through macros
 * 	  given in this file. This simplified performing operations
 * 	  on lattice data.
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Geometry/geometry_type.h"
#include "Geometry/geometry_indexing.h"
#include "Geometry/communications.h"
#include "Geometry/communications_reduced.h"
#include "Geometry/cpu_geometry.h"
#include "Geometry/geometry_check.h"
#include "Geometry/geometry_descriptor.h"
#include "Geometry/geometry_fuse.h"
#include "Geometry/geometry_init.h"
#include "Geometry/geometry_maskstate.h"
#include "Geometry/geometry_omp.h"
#include "Geometry/new_geometry.h"
#include "Geometry/setup.h"
#include "Geometry/strided_reads.h"
#include "Geometry/hr_sendrecv.h"

#ifdef WITH_GPU
#include "Geometry/gpu_affinity.h"
#include "Geometry/geometry_gpu_init.h"
#include "Geometry/gpu_geometry.h"
#include "Geometry/strided_reads_gpu.h"
#include "Geometry/read_clover.h"
#endif

/* this define the width of the borders for parallel dimensions
 * For different actions it must be modified
 * FIXME: Put into global variables, or put more geometry variables from global to
 * here.
 */
#define BORDERSIZE 1

#endif
