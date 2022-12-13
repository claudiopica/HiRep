#ifndef NEW_GEOM_GPU_H
#define NEW_GEOM_GPU_H

#include "new_geometry.h"

#ifdef WITH_GPU
void sync_box_to_buffer_gpu_spinor_field_f(geometry_descriptor*, box_t*, suNf_spinor*, void*);
void sync_box_to_buffer_gpu_gfield_f(geometry_descriptor*, box_t*, suNf*, void*);

void sync_buffer_to_box_gpu_spinor_field_f(geometry_descriptor*, box_t*, suNf_spinor*, void*);
void sync_buffer_to_box_gpu_gfield_f(geometry_descriptor*, box_t*, suNf*, void*);
#endif

#endif