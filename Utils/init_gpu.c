/***************************************************************************\
* Copyright (c) 2012, Ari Hietanen                                          * 
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File init_gpu.c
*
* Functions to initialize and setup the global variables for GPU
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gpu.h"
#include "io.h"
#include "global.h"
#include "logger.h"
#include "utils.h"
#include "error.h"


void init_gpu(input_gpu gpu_var){
  int device_count = 0;
  cudaError_t error_id;
  cudaDeviceProp device_prop;
  int mem_clock,mem_bus_width, l2_cache_size,driver_version,runtime_version; /*Device properties*/

  const char *sComputeMode[] = {
    "Default (multiple host threads can use ::cudaSetDevice() with device simultaneously)",
    "Exclusive (only one host thread in one process is able to use ::cudaSetDevice() with this device)",
    "Prohibited (no host thread can use ::cudaSetDevice() with this device)",
    "Exclusive Process (many threads in one process is able to use ::cudaSetDevice() with this device)",
    "Unknown",
    NULL
  };

  lprintf("GPU_INIT",0,"Initializing GPU\n");

  gpu_id = gpu_var.gpuID;
  /*Select the right GPU */
  error_id = cudaGetDeviceCount(&device_count);
  error(error_id != cudaSuccess,1,"init_gpu","Error getting device count");
  error(gpu_id>=device_count,1,"init_gpu","Illegal device ID");
  error_id = cudaSetDevice(gpu_id);
  error(error_id != cudaSuccess,1,"init_gpu","Error selecting GPU");
  lprintf("GPU_INIT",0,"Using GPU #%d\n",gpu_id);

  /*Print out the device info assume CUDART >= 4000*/
  cudaGetDeviceProperties(&device_prop, gpu_id);
  lprintf("GPU_INIT",10,"Device: %d\n",device_prop.name);
  cudaDriverGetVersion(&driver_version);
  cudaRuntimeGetVersion(&runtime_version);
  lprintf("GPU_INIT",10,"CUDA Driver Version / Runtime Version: %d.%d / %d.%d\n", driver_version/1000, driver_version%100, runtime_version/1000, runtime_version%100);
  lprintf("GPU_INIT",10,"CUDA Capability Major/Minor version number: %d.%d\n", device_prop.major, device_prop.minor);

  lprintf("GPU_INIT",10,"Total amount of global memory: %.0f MB (%lluB)\n", 
	  (float)device_prop.totalGlobalMem/((float) (1<<20)), (unsigned long long) device_prop.totalGlobalMem);

  lprintf("GPU_INIT",10,"Multiprocessors: %d\n", device_prop.multiProcessorCount);
  /*,
	 ConvertSMVer2Cores(device_prop.major, device_prop.minor),
	 ConvertSMVer2Cores(device_prop.major, device_prop.minor) * device_prop.multiProcessorCount);*/

  lprintf("GPU_INT",10,"GPU Clock Speed: %.2f GHz\n", device_prop.clockRate * 1e-6f);
  cuDeviceGetAttribute( &mem_clock,CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE , gpu_id);
  lprintf("GPU_INIT",10,"Memory Clock rate: %.3f Mhz\n",mem_clock*1e-3f); 
  cuDeviceGetAttribute(&mem_bus_width, CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH,gpu_id);
  lprintf("GPU_INIT",10,"Memory Bus Widht: %d-bit\n",mem_bus_width); 
  cuDeviceGetAttribute(&l2_cache_size, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, gpu_id);
  lprintf("GPU_INIT",10,"L2 Cache Size: %dB\n",l2_cache_size); 
  
  lprintf("GPU_INIT",10,"Max Texture dimension size (x,y,z): 1D=(%d), 2D=(%d,%d), 3D=(%d,%d,%d)\n",
	  device_prop.maxTexture1D,device_prop.maxTexture2D[0],device_prop.maxTexture2D[1],
	  device_prop.maxTexture3D[0],device_prop.maxTexture3D[1],device_prop.maxTexture3D[2]);

  lprintf("GPU_INIT",10,"Max Layered Texture dimension size (dim) x layers: 1D=(%d) x %d, 2D=(%d,%d) x %d\n",
	  device_prop.maxTexture1DLayered[0],device_prop.maxTexture1DLayered[1],
	  device_prop.maxTexture2DLayered[0],device_prop.maxTexture2DLayered[1],device_prop.maxTexture2DLayered[2]);

  lprintf("GPU_INIT",10,"Total amount of constant memory: %dB\n",device_prop.totalConstMem);
  lprintf("GPU_INIT",10,"Total amount of shared memory per block: %dB\n",device_prop.sharedMemPerBlock); 
  lprintf("GPU_INIT",10,"Total number of register per block: %dB\n",device_prop.regsPerBlock); 
  lprintf("GPU_INIT",10,"Warp size: %dB\n",device_prop.warpSize); 
  error(device_prop.warpSize!=32,1,"init_gpu","Error: warp size 32 assumed in global sum\n");  
  lprintf("GPU_INIT",10,"Maximum number of threds per block: %d\n",device_prop.maxThreadsPerBlock);   
  lprintf("GPU_INIT",10,"Maximum size of each dimension of a block (x,y,z): (%d,%d,%d)\n",
	  device_prop.maxThreadsDim[0],device_prop.maxThreadsDim[1],device_prop.maxThreadsDim[2]);
  lprintf("GPU_INIT",10,"Maximum size of each dimension of a grid (x,y,z): (%d,%d,%d)\n",
	  device_prop.maxGridSize[0],device_prop.maxGridSize[1],device_prop.maxGridSize[2]);

  lprintf("GPU_INIT",10,"Maximum memory pitch: %uB\n",device_prop.memPitch);
  lprintf("GPU_INIT",10,"Texture alignment: %uB\n",device_prop.textureAlignment);

  lprintf("GPU_INIT",10,"Concurrent copy and execution: %s with %d copy engine(s)\n", (device_prop.deviceOverlap ? "Yes" : "No"), device_prop.asyncEngineCount);

  lprintf("GPU_INIT",10,"Run time limit on kernels:                     %s\n", device_prop.kernelExecTimeoutEnabled ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Integrated GPU sharing Host Memory:            %s\n", device_prop.integrated ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Support host page-locked memory mapping:       %s\n", device_prop.canMapHostMemory ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Concurrent kernel execution:                   %s\n", device_prop.concurrentKernels ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Alignment requirement for Surfaces:            %s\n", device_prop.surfaceAlignment ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Device has ECC support enabled:                %s\n", device_prop.ECCEnabled ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Device is using TCC driver mode:               %s\n", device_prop.tccDriver ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Device supports Unified Addressing (UVA):      %s\n", device_prop.unifiedAddressing ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Device PCI Bus ID / PCI location ID:           %d / %d\n", device_prop.pciBusID, device_prop.pciDeviceID );
  lprintf("GPU_INIT",10,"Compute Mode:\n");
  lprintf("GPU_INIT",10,"  < %s >\n", sComputeMode[device_prop.computeMode]);

  /* Set up the global variables */
  grid_size_max_gpu = device_prop.maxGridSize[0];
 

}
