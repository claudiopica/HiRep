#ifndef GLOBAL_GPU_H
#define GLOBAL_GPU_H

#ifdef MAIN_PROGRAM
#  define GLB_VAR_DEV(type,name,...) __device__ type name __VA_ARGS__
#else
#  define GLB_VAR_DEV(type,name,...) extern __device__ type name
#endif

#ifdef WITH_GPU
GLB_VAR_DEV(int,T_EXT_GPU,=0);
GLB_VAR_DEV(int,X_EXT_GPU,=0);
GLB_VAR_DEV(int,Y_EXT_GPU,=0);
GLB_VAR_DEV(int,Z_EXT_GPU,=0);

#define ipt_ext_gpu(t, x, y, z) ipt_gpu[_lexi(T_EXT_GPU, X_EXT_GPU, Y_EXT_GPU, Z_EXT_GPU, t, x, y, z)]
#endif

#endif