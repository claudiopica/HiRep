/**
 * @file cpu_geometry.h
 * @brief This file contains macros to load elements of single sites of a field.
 */

#ifndef CPU_GEOMETRY_H
#define CPU_GEOMETRY_H

#include <stdlib.h>

#ifdef __cplusplus
	extern "C" {
#endif

/// Returns x mod y, where y is assumed positive and x can have any
/// sign. The returned value is in the interval [0,y)
inline static int safe_mod(int x, int y) {
  if (x >= 0) return (x % y);
  else return ((y - (abs(x) % y)) % y);
}

inline static int safe_mod_alt(int x, int y) {
  while (x<0) x+=y;
  while (x>=y) x-=y;
  return x;
}

/// Returns x mod y, where y is assumed positive and x can have any
/// sign. The returned value is in the interval [0,y)
/// This assumes -y <= x < 2y
inline static int safe_mod_fast(int x, int y) {
  if (x < 0) return x+y;
  else return (x < y) ? x : x-y;
}

/* NB: it is assumed in the code that different directions are contiguous in memory */
#define coord_to_index(ix,mu) (((ix)<<2)|(mu))
#define index_to_coord(i,ix,mu) (mu)=((i)&3);(ix)=((i)>>2)


#define _FIELD_AT(s,i) (((s)->ptr) + i - (s)->type->master_shift)
#define _4FIELD_AT(s,i,mu) (((s)->ptr) + coord_to_index(i-(s)->type->master_shift,mu))
#define _6FIELD_AT(s,i,mu) (((s)->ptr) + (( i - (s)->type->master_shift)*6+mu))
#define _DFIELD_AT(s,i,mu,size) (size == 1) ? _FIELD_AT(s,i) : \
				((size == 4) ? _4FIELD_AT(s,i,mu) : \
				 _6FIELD_AT(s,i,mu))

#define _FIELD_AT_PTR(s,i,_master_shift) (s + i - _master_shift)
#define _4FIELD_AT_PTR(s,i,mu,_master_shift) (s + coord_to_index(i-_master_shift,mu))
#define _6FIELD_AT_PTR(s,i,mu,_master_shift) (s + ((i - _master_shift)*6+mu))
#define _DFIELD_AT_PTR(s,i,mu,_master_shift,__size) (__size == 1) ? _FIELD_AT_PTR(s,i,_master_shift) : \
				((__size == 4) ? _4FIELD_AT_PTR(s,i,mu,_master_shift) : \
				 _6FIELD_AT_PTR(s,i,mu,_master_shift))

#define START_SP_ADDRESS_GPU(sf) ((sf)->gpu_ptr + (sf)->type->master_start[0])
#define START_GF_ADDRESS_GPU(gf) ((gf)->gpu_ptr + (gf)->type->master_start[0])

#define _GPU_FIELD_BLK(s,i) (((s)->gpu_ptr) + (s)->type->master_start[(i)])
#define _GPU_4FIELD_BLK(s,i) (((s)->gpu_ptr) + 4*(s)->type->master_start[(i)])
#define _GPU_DFIELD_BLK(s,i,size) (((s)->gpu_ptr) + size*(s)->type->master_start[(i)])

//#define _GPU_FIELD_BLK_BUF(s,i) (((s)->gpu_ptr) + (s)->type->master_end[(i)])
//#define _GPU_4FIELD_BLK_BUF(s,i) (((s)->gpu_ptr) + 4*(s)->type->master_end[(i)])
//#define _GPU_DFIELD_BLK_BUF(s,i,size) (((s)->gpu_ptr) + size*(s)->type->master_end[(i)])

#define _FIELD_BLK(s,i) (((s)->ptr) + (s)->type->master_start[(i)])
#define _4FIELD_BLK(s,i) (((s)->ptr) + 4*(s)->type->master_start[(i)])
#define _DFIELD_BLK(s,i,size) (((s)->ptr) + size*(s)->type->master_start[(i)]) 

#define _BUF_FIELD_BLK(s,i) (((s)->ptr) + (s)->type->rbuf_start[(i)])
#define _BUF_4FIELD_BLK(s,i) (((s)->ptr) + 4*(s)->type->rbuf_start[(i)])
#define _BUF_DFIELD_BLK(s,i,_size) (((s)->ptr) + (_size)*(s)->type->rbuf_start[(i)]) 

#define _BUF_GPU_FIELD_BLK(s,i) (((s)->gpu_ptr) + (s)->type->rbuf_start[(i)])
#define _BUF_GPU_4FIELD_BLK(s,i) (((s)->gpu_ptr) + 4*(s)->type->rbuf_start[(i)])
#define _BUF_GPU_DFIELD_BLK(s,i,size) (((s)->gpu_ptr) + size*(s)->type->rbuf_start[(i)])

#define _GPU_IDX_TO_LOCAL(in, ix, ixp) ix - in->type->master_start[(ixp)];
//#define _SITE_IDX_GPU(ix, ixp, stride) (ix) + stride*(ixp)

#ifdef __cplusplus
	}
#endif
#endif