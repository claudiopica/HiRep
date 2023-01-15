/***************************************************************************\
* Copyright (c) 2022, Claudio Pica, Sofie Martins                           *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file new_geometry.h
 * @brief Functions needed for the new geometry implementation that
 *        will replace the current geometry in the future
 */

#ifndef NEW_GEOMETRY_H
#define NEW_GEOMETRY_H

#include "spinor_field.h"
#include "suN.h"
#include "suN_types.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
    extern "C" {
#endif

void define_geometry(void);
void* sendbuf_alloc(size_t bytes_per_site);
#ifdef WITH_GPU
void* sendbuf_alloc_gpu(size_t bytes_per_site);
#endif
void sync_field(geometry_descriptor *gd, int byte_per_site, int is_spinor_like, void *latticebuf, void *sb_ptr);
int test_define_geometry(void);
void sendbuf_report(void);
void sync_field_gpu(geometry_descriptor*, int, int, void*, void*);

typedef struct coord4 {
    uint8_t x[4];
} coord4;

enum box_type {
    // L0 = 0, //unused
    // L1 = 1, //unused
    L2 = 2,
    L3 = 3,
    INNER = 4,
    SENDBUF = 5
};

// TODO: this should be in geometry.h and geometry descriptor should contain it
// enum to define geometry type
// this is a simple bitmask with GLOBAL = EVEN | ODD
enum gd_type {
    EVEN   = 1, 
    ODD    = 2,
    GLOBAL = 3
};

//  ----h
//  ....|  NB: the h[4] is not in the box,   
//  ....|  i.e. coordinates needs to be l[]<= p[] <h[] 
//  l...|
typedef struct box_t {
    int l[4]; // the lower left corner, the box base point (e.g. in the extended lattice)
    int h[4]; // the upper right corner
    int base_index;
    int base_index_odd;
    int parity; //0 -> base point is even; 1 -> basepoint is odd
    char mask; //tells if the box is a border, e.g. if T_UP_MASK is set the box is in top T border of the extended lattice
    enum box_type type; // tell the type of the box, just a convenience for testing
    int *ipt_ext; //given the cordinate of a point in the box returns an index
    coord4 *icoord; //given an index in the box return the 4D coordinates of the point in the box relative to the l[4]
    int *icoord_idx; //given an index in the box return the index of the point in the box relative to the l[4]
    struct box_t *sendBox; //if this is a border corresponding to a Recv buffer, this is the box to copy data from, i.e. corresponding to the Send buffer
    struct box_t *next; // link to next box. NULL if last
} box_t;
//TODO: do we want to add vol, even_vol, odd_vol for avoid recomputing them every time?
//TODO: do we want to precompute ipt_ext for sendboxes?


int boxEvenVolume(box_t *B);
int boxOddVolume(box_t *B);
int boxVolume(box_t*);

#define _DECLARE_SYNC_TO_BUFFER(_name, _type) \
    void sync_box_to_buffer_gpu_##_name(geometry_descriptor*,box_t*,void*,void*); \
    void sync_field_to_buffer_##_name(geometry_descriptor*,void*,void*);

_DECLARE_SYNC_TO_BUFFER(spinor_field_f, suNf_spinor);
_DECLARE_SYNC_TO_BUFFER(spinor_field_f_flt, suNf_spinor_flt);
_DECLARE_SYNC_TO_BUFFER(sfield, double);

_DECLARE_SYNC_TO_BUFFER(gfield, suNg);
_DECLARE_SYNC_TO_BUFFER(gfield_f, suNf);
_DECLARE_SYNC_TO_BUFFER(gfield_flt, suNg_flt);
_DECLARE_SYNC_TO_BUFFER(gfield_f_flt, suNf_flt);
_DECLARE_SYNC_TO_BUFFER(suNg_scalar_field, suNg_vector);
_DECLARE_SYNC_TO_BUFFER(avfield, suNg_algebra_vector);
_DECLARE_SYNC_TO_BUFFER(gtransf, suNg);
_DECLARE_SYNC_TO_BUFFER(clover_ldl, ldl_t);
_DECLARE_SYNC_TO_BUFFER(clover_term, suNfc);
_DECLARE_SYNC_TO_BUFFER(clover_force, suNf);

#ifdef __cplusplus
    }
#endif

#endif