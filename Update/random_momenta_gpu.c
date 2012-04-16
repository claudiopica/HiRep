/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/
#ifdef WITH_GPU
#include "spinor_field.h"
#include "suN_repr_func.h"
#include "random.h"
#include "update.h"
#include "utils.h"
#include "memory.h"
#include <math.h>

void gaussian_momenta(suNg_av_field *momenta) {
  gaussian_momenta_cpu(momenta);
  suNg_av_field_copy_to_gpu(momenta);
}

void zero_momenta(suNg_av_field *momenta) {
  zero_momenta_cpu(momenta);
  suNg_av_field_copy_to_gpu(momenta);
}
#endif //WITH_GPU
