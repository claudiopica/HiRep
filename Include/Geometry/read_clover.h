/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef READ_CLOVER_H
#define READ_CLOVER_H

#ifdef __cplusplus

__device__ void read_clover(hr_complex *c, const ldl_t *in, int ix, int dn, int j);
__device__ void write_clover(hr_complex *c, ldl_t *out, int ix, int dn, int j);
__device__ void read_force(hr_complex *c, const suNf *in, int ix, int comp, int i, int j);
__device__ void write_force(hr_complex *c, suNf *out, int ix, int comp, int i, int j);
__device__ void read_clover_term_comp(hr_complex *c, const suNfc *in, int ix, int comp, int i, int j);
__device__ void write_clover_term_comp(hr_complex *c, suNfc *out, int ix, int comp, int i, int j);

#endif
#endif
