#ifndef EVA_GPU_HPP
#define EVA_GPU_HPP

#ifdef WITH_GPU

#include "libhr_core.h"

#ifdef __cplusplus
extern "C" {
#endif

void rotate_gpu(int, spinor_field *, hr_complex *);

#ifdef __cplusplus
}
#endif

#endif
#endif