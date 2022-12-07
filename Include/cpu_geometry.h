/**
 * @file cpu_geometry.h
 * @brief This file contains macros to load elements of single sites of a field.
 */

#ifndef CPU_GEOMETRY_H
#define CPU_GEOMETRY_H

#include "field_ordering.h"

#ifdef __cplusplus
	extern "C" {
#endif

#define _FIELD_AT(s,i) (((s)->ptr) + i - (s)->type->master_shift)
#define _4FIELD_AT(s,i,mu) (((s)->ptr) + coord_to_index(i-(s)->type->master_shift,mu))
#define _6FIELD_AT(s,i,mu) (((s)->ptr) + (( i - (s)->type->master_shift)*6+mu))
#define _DFIELD_AT(s,i,mu,size) (size == 1) ? _FIELD_AT(s,i) : \
				((size == 4) ? _4FIELD_AT(s,i,mu) : \
				 _6FIELD_AT(s,i,mu))

#define START_SP_ADDRESS_GPU(sf) ((sf)->gpu_ptr + (sf)->type->master_start[0])
#define START_GF_ADDRESS_GPU(gf) ((gf)->gpu_ptr + (gf)->type->master_start[0])

#define _GPU_FIELD_BLK(s,i) (((s)->gpu_ptr) + (s)->type->master_start[(i)])
#define _GPU_4FIELD_BLK(s,i) (((s)->gpu_ptr) + 4*(s)->type->master_start[(i)])
#define _GPU_DFIELD_BLK(s,i,size) (((s)->gpu_ptr) + size*(s)->type->master_start[(i)])

#define _FIELD_BLK(s,i) (((s)->ptr) + (s)->type->master_start[(i)])
#define _4FIELD_BLK(s,i) (((s)->ptr) + 4*(s)->type->master_start[(i)])
#define _DFIELD_BLK(s,i,size) (((s)->ptr) + size*(s)->type->master_start[(i)]) 

#define _BUF_FIELD_BLK(s,i) (((s)->ptr) + (s)->type->sbuf_start[(i)])
#define _BUF_4FIELD_BLK(s,i) (((s)->ptr) + 4*(s)->type->sbuf_start[(i)])
#define _BUF_DFIELD_BLK(s,i,_size) (((s)->ptr) + (_size)*(s)->type->sbuf_start[(i)]) 

#define _BUF_GPU_FIELD_BLK(s,i) (((s)->gpu_ptr) + (s)->type->sbuf_start[(i)])
#define _BUF_GPU_4FIELD_BLK(s,i) (((s)->gpu_ptr) + 4*(s)->type->sbuf_start[(i)])
#define _BUF_GPU_DFIELD_BLK(s,i,size) (((s)->gpu_ptr) + size*(s)->type->sbuf_start[(i)])

#define _GPU_IDX_TO_LOCAL(in, ix, ixp) ix - in->type->master_start[(ixp)];
#define _SITE_IDX_GPU(ix, ixp, stride) (ix) + stride*(ixp)

#ifdef __cplusplus
	}
#endif
#endif
