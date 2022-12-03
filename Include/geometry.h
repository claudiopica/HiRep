/***************************************************************************\
 * Copyright (c) 2008-2014, Claudio Pica                                    *
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

#include "hr_omp.h"
#include "geometry_omp.h"
#include "geometry_fuse.h"
#include "geometry_init.h"
#include "cpu_geometry.h"
#include "geometry_check.h"
#ifdef WITH_GPU
   #include "gpu_geometry.h"
#endif

/* this define the width of the borders for parallel dimensions
 * For different actions it must be modified
 * FIXME: Put into global variables, or put more geometry variables from global to
 * here.
 */

#define BORDERSIZE 1

/**
 * @struct geometry_descriptor
 * @brief This struct should contain all information necessary to perform operations
 *        on the local lattice including correct inter-core and inter-node communications.
 *        More information on GPU geometry with illustrations can be found in the 
 *        developer guide in the section on GPU geometry.
 *
 * @var geometry_descriptor::inner_master_pieces
 * 				Pieces that contain only bulk elements and can therefore
 * 				need no communications
 * @var geometry_descriptor::local_master_pieces
 * 				Pieces that contain only sites of the local lattice. 
 * @var geometry_descriptor::total_spinor_master_pieces
 * 				Total number of pieces for a spinor field. Since we need
 * 				more buffers for the gauge fields, this number is smaller 
 * 				than the number of gauge master pieces and the buffer pieces
 * 				are a subset of the buffer pieces for gauge field 
 * 				communications.
 * @var geometry_descriptor::total_gauge_master_pieces
 * 				Total number of pieces necessary for communications of a
 * 				gauge field. Since we need more buffers for the gauge fields
 * 				this number is bigger than the number of spinor master
 * 				pieces and the buffer pieces for the gauge fields are a 
 * 				superset of the buffer pieces for the spinor field.
 * @var geometry_descriptor::master_start
 * 				Array that takes the piece number as index and returns
 * 				the index of the first site in the field data that is 
 * 				element of the piece with the given index.
 * 				It is assumed that other elements on the piece are stored
 * 				on the following indices.
 * @var geometry_descriptor::master_end
 * 				Array that takes the piece number as index and returns
 * 				the index of the last site in field data that is element
 * 				of the piece with the given index.
 * 				It is assumed that other elements on the piece are stored
 * 				in preceding indices.
 * @var geometry_descriptor::master_shift
 * 				This gives the shift of spinor relative to a full lattice.
 * 				For example: A spinor with odd parity has a geometry that
 * 				is deplaced by half the number of lattice sites relative 
 * 				to a spinor that is defined on the full lattice.
 * @var geometry_descriptor::ncopies_spinor
 * 				This gives the number of pieces that need to be copied
 * 				from the local lattice to the local lattice in order
 * 				to ensure contingency of piece data in memory for the 
 * 				spinor field.
 * @var geometry_descriptor::ncopies_gauge
 * 				This gives the number of pieces that need to be copied
 * 				from the local lattice to the local lattice in order
 * 				to ensure contingency of piece data in memory for the 
 * 				gauge field.
 * @var geometry_descriptor::nbuffers_spinor
 * 				Number of buffers that need to be transferred for the 
 * 				spinor field.
 * @var geometry_descriptor::nbuffers_gauge
 * 				Number of buffers that need to be transferred for the gauge
 * 				field.
 * @var geometry_descriptor::rbuf_len
 * 				Array that takes a buffer id and then returns the length
 * 				of this buffer on the receive side.
 * @var geometry_descriptor::sbuf_len
 * 				Array that takes a buffer id and then returns the length
 * 				of this buffer on the send site.
 * @var geometry_descriptor::rbuf_from_proc
 * 				Array that takes a buffer id and then return the process id
 * 				of the target process.
 * @var geometry_descriptor::rbuf_start
 * 				Array that takes a buffer id and returns an index of the 
 * 				first site in memory that contains a site that is an element
 * 				of the buffer piece.
 * @var geometry_descriptor::sbuf_to_proc
 *				Array that takes a buffer id and returns the process id of
 *				the process that sends the corresponding buffer.
 * @var geometry_descriptor::sbuf_start
 * 				Array that takes a buffer id and returns an index of the 
 * 				first site in memory that contains a site that is an 
 * 				element of the send buffer piece.
 * FIXME: Complete this...
 *
 *
 */
typedef struct
{
  int inner_master_pieces; 
  int local_master_pieces; 
  int total_spinor_master_pieces;
  int total_gauge_master_pieces;
  int *master_start, *master_end; 
  int master_shift;
  int ncopies_spinor;
  int ncopies_gauge;
  int *copy_from, *copy_to, *copy_len;
  int copy_shift;
  int nbuffers_spinor;
  int nbuffers_gauge;
  int *rbuf_len, *sbuf_len;
  int *rbuf_from_proc, *rbuf_start;
  int *sbuf_to_proc, *sbuf_start;
  int gsize_spinor;
  int gsize_gauge;
  int *fuse_mask;
  int fuse_gauge_size;
  int fuse_inner_counter;
} geometry_descriptor;

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
  _OMP_PRAGMA(_omp_parallel)                        \
  _OMP_PRAGMA(_omp_for redop1 redop2)               \
  for (int is = (type)->master_start[ip]; is <= (type)->master_end[ip]; is++)

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
 * @brief Iterate over all sites of the local lattice
 *
 * @param type		geometry_descriptor that contains information on the geometry
 * 			of the local lattice
 * @param is		Local variable that runs over all site indices
 */
#define _MASTER_FOR(type, is) _MASTER_FOR_RED(type, is, , )

/**
 * @brief Iterate over all sites of the local lattice but not by index in memory but by
 *        spinor. The current spinor can be found using _SPINOR_PTR(s). 
 *
 * @param s		Input spinor field
 */
#define _ONE_SPINOR_FOR(s) _ONE_SPINOR_FOR_RED(s,,)

/**
 * @brief Iterate over two corresponding spinors on the given fields. The current spinors 
 * 	  can be found using _SPINOR_PTR(s1) and _SPINOR_PTR(s2).
 *
 * @param s1		First input spinor field
 * @param s2		Second input spinor field
 */
#define _TWO_SPINORS_FOR(s1,s2) _TWO_SPINORS_FOR_RED(s1,s2,,)

/**
 * @brief Iterate over all three corresponding spinors on the given fields. The current 
 * 	  spinors can be found using _SPINOR_PTR(s1), _SPINOR_PTR(s2), _SPINOR_PTR(s3).
 *
 * @param s1		First input spinor field
 * @param s2		Second input spinor field
 * @param s3		Third input spinor field
 */
#define _THREE_SPINORS_FOR(s1,s2,s3) _THREE_SPINORS_FOR_RED(s1,s2,s3,,)

/**
 * @brief Retrieve current spinor field. This macro only works insite _SPINOR_FOR, 
 * 	  _TWO_SPINORS_FOR or _THREE_SPINORS_FOR.
 *
 * @param s 		Spinor field that is being iterated over.
 */
#define _SPINOR_PTR(s) _FIELD_AT(s,_spinor_for_is)

#endif
