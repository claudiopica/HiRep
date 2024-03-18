/***************************************************************************\
* Copyright (c) 2022, Claudio Pica, Sofie Martins                           *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file geometry_descriptor.h
 * @brief Geometry descriptor struct, that contains all necessary information for multi-node/
 *        multi-GPU calculations.
 */

#ifndef GEOMETRY_DESCRIPTOR_H
#define GEOMETRY_DESCRIPTOR_H
#ifdef __cplusplus
extern "C" {
#endif

#include "Geometry/geometry_type.h"

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
 * @var geometry_descriptor::gsize_spinor
 *        Number of sites allocated for a spinor field
 * @var geometry_descriptor::gize_gauge
 *        Number of links allocated for a gauge field
 * @var geometry_descriptor::fuse_mask
 * @var geometry_descriptor::fuse_gauge_size
 * @var geometry_descriptor::fuse_inner_counter
 *
 */
typedef struct geometry_descriptor {
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
#ifdef WITH_NEW_GEOMETRY
    gd_type desc;
#endif
} geometry_descriptor;

#ifdef __cplusplus
}
#endif
#endif
