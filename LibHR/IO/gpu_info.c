/***************************************************************************\
* Copyright (c) 2012, Ari Hietanen                                          * 
* All rights reserved.                                                      * 
\***************************************************************************/
 
/**
 * @file gpu_info.c
 * @brief Prints information on available hardware and software on the 
 *        cluster relating to GPUs.
 */
#ifdef WITH_GPU
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gpu.h"
#include "io.h"
#include "global.h"
#include "logger.h"
#include "utils.h"
#include "error.h"
#include "geometry.h"
#include "gpu.h"
#ifdef WITH_MPI
  #include "mpi.h"
#endif

#ifdef __cplusplus

const char *sComputeMode[] = {
    "Default (multiple host threads can use ::cudaSetDevice() with device simultaneously)",
    "Exclusive (only one host thread in one process is able to use ::cudaSetDevice() with this device)",
    "Prohibited (no host thread can use ::cudaSetDevice() with this device)",
    "Exclusive Process (many threads in one process is able to use ::cudaSetDevice() with this device)",
    "Unknown",
    NULL
  };

/**
 * @brief Query number of GPUs and print related information
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_device_count_info(input_gpu gpu_var) 
{
  int device_count;
  CHECK_CUDA(cudaGetDeviceCount(&device_count));
  lprintf("GPU_INIT", 0, "GPU_ID = %d\n", gpu_var.gpuID);
  //error(gpu_id > device_count, 1, "init_gpu", "Illegal device ID"); 
  // I don't see what we need this for (SAM)
}

/**
 * @brief Print CUDA driver version information
 */
void print_driver_info() 
{
  int driver_version;
  CHECK_CUDA(cudaDriverGetVersion(&driver_version));
  lprintf("GPU_INIT", 10, "CUDA Driver Version / Runtime Version: %d.%d / %d.%d\n", 
                           driver_version/1000, 
                           driver_version%100, 
                           runtime_version/1000, 
                           runtime_version%100);
}

/**
 * @brief Print CUDA runtime version information
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_runtime_info(cudaDeviceProp device_prop) 
{
  int runtime_version;
  CHECK_CUDA(cudaRuntimeGetVersion(&runtime_version));
  lprintf("GPU_INIT", 10, "CUDA Capability Major/Minor version number: %d.%d\n", 
                          device_prop.major, device_prop.minor);
}

/**
 * @brief Print Global memory information including bandwidth 
 *        paramters and supported features
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_global_memory_info(cudaDeviceProp device_prop, input_gpu gpu_var) 
{
  int mem_bus_width;

  // Query properties
  CHECK_CUDA(cuDeviceGetAttribute(&mem_bus_width, CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH,gpu_var.gpuID));

  // Print formatted
  lprintf("GPU_INIT",10,"Total amount of global memory: %.0f MB (%lluB)\n", 
	                      (float)device_prop.totalGlobalMem/((float) (1<<20)), 
                        (unsigned long long) device_prop.totalGlobalMem);
  lprintf("GPU_INIT",10,"Memory Clock rate: %.3f Mhz\n",mem_clock*1e-3f); 
  lprintf("GPU_INIT",10,"Memory Bus Width: %d-bit\n",mem_bus_width); 
  lprintf("GPU_INIT",10,"Maximum memory pitch: %uB\n",device_prop.memPitch);
  lprintf("GPU_INIT",10,"Integrated GPU sharing Host Memory: %s\n", 
                          device_prop.integrated ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Support host page-locked memory mapping: %s\n", 
                          device_prop.canMapHostMemory ? "Yes" : "No");
}

/**
 * @brief Print information on shared memory
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_shared_memory_info(cudaDeviceProp device_prop) 
{
  lprintf("GPU_INIT",10,"Total amount of shared memory per block: %dB\n",device_prop.sharedMemPerBlock); 
}

/**
 * @brief Print (L2) cache information
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_cache_info(cudaDeviceProp device_prop, input_gpu gpu_var) 
{
  int l2_cache_size;

  // Query properties
  CHECK_CUDA(cuDeviceGetAttribute(&l2_cache_size, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, gpu_var.gpuID));

  // Print formatted
  lprintf("GPU_INIT",10,"L2 Cache Size: %dB\n",l2_cache_size); 
  // TODO: L1 cache??
}

/**
 * @brief Print information on constant memory, in particular
 *        amount available, alignment, texture and layered texture 
 *        memory paramters
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_constant_memory_info(cudaDeviceProp device_prop) 
{
  lprintf("GPU_INIT",10,"Max Texture dimension size (x,y,z): 1D=(%d), 2D=(%d,%d), 3D=(%d,%d,%d)\n",
	  device_prop.maxTexture1D, device_prop.maxTexture2D[0], device_prop.maxTexture2D[1],
	  device_prop.maxTexture3D[0], device_prop.maxTexture3D[1], device_prop.maxTexture3D[2]);

  lprintf("GPU_INIT",10,"Max Layered Texture dimension size (dim) x layers: 1D=(%d) x %d, 2D=(%d,%d) x %d\n",
	  device_prop.maxTexture1DLayered[0], device_prop.maxTexture1DLayered[1],
	  device_prop.maxTexture2DLayered[0], device_prop.maxTexture2DLayered[1], 
    device_prop.maxTexture2DLayered[2]);
  
  lprintf("GPU_INIT",10,"Texture alignment: %uB\n",device_prop.textureAlignment);
  lprintf("GPU_INIT",10,"Total amount of constant memory: %dB\n",device_prop.totalConstMem);
  lprintf("GPU_INIT",10,"Alignment requirement for Surfaces: %s\n", 
                            device_prop.surfaceAlignment ? "Yes" : "No");
}
 
/**
 * @brief Prints all memory related info, here on global, shared
 *        memory, cache and constant memory, see resp. functions
 *        for more information.
 * 
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_memory_info(cudaDeviceProp device_prop, input_gpu gpu_var) 
{
  print_global_memory_info(device_prop, gpu_var);
  print_shared_memory_info(device_prop, gpu_var);
  print_cache_info(device_prop, gpu_var);
  print_constant_memory_info(device_prop, gpu_var);
}

/**
 * @brief Prints information on compute parameters
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_compute_info(cudaDeviceProp device_prop, input_gpu gpu_var) 
{
  int mem_clock;

  // Query properties
  CHECK_CUDA(cuDeviceGetAttribute(&mem_clock, CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE, gpu_var.gpuID));

  // Print formatted
  lprintf("GPU_INIT", 10, "Multiprocessors: %d\n", device_prop.multiProcessorCount);
  lprintf("GPU_INT", 10, "GPU Clock Speed: %.2f GHz\n", device_prop.clockRate * 1e-6f);
  lprintf("GPU_INIT", 10, "Total number of register per block: %dB\n", device_prop.regsPerBlock); 
  lprintf("GPU_INIT",10, "Warp size: %dB\n", device_prop.warpSize); 
  error(device_prop.warpSize!=32, 1, "init_gpu", "Error: warp size 32 assumed in global sum\n");  
  lprintf("GPU_INIT", 10, "Maximum number of threads per block: %d\n", device_prop.maxThreadsPerBlock);   
  lprintf("GPU_INIT", 10, "Maximum size of each dimension of a block (x,y,z): (%d,%d,%d)\n",
	                        device_prop.maxThreadsDim[0], 
                          device_prop.maxThreadsDim[1],
                          device_prop.maxThreadsDim[2]);
  lprintf("GPU_INIT", 10, "Maximum size of each dimension of a grid (x,y,z): (%d,%d,%d)\n",
	                        device_prop.maxGridSize[0],
                          device_prop.maxGridSize[1],
                          device_prop.maxGridSize[2]);
  lprintf("GPU_INIT",10,"Concurrent copy and execution: %s with %d copy engine(s)\n", 
                          device_prop.deviceOverlap ? "Yes" : "No", 
                          device_prop.asyncEngineCount);
  lprintf("GPU_INIT",10,"Run time limit on kernels: %s\n", 
                          device_prop.kernelExecTimeoutEnabled ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Concurrent kernel execution: %s\n", 
                          device_prop.concurrentKernels ? "Yes" : "No");
}

/**
 * @brief Checks for a number of other supported features and
 *        prints information on them.
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_supported_features(cudaDeviceProp device_prop) 
{
  lprintf("GPU_INIT",10,"Device has ECC support enabled:                %s\n", device_prop.ECCEnabled ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Device is using TCC driver mode:               %s\n", device_prop.tccDriver ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Device supports Unified Addressing (UVA):      %s\n", device_prop.unifiedAddressing ? "Yes" : "No");
  lprintf("GPU_INIT",10,"Device PCI Bus ID / PCI location ID:           %d / %d\n", device_prop.pciBusID, device_prop.pciDeviceID );
  lprintf("GPU_INIT",10,"Compute Mode:\n");
  lprintf("GPU_INIT",10,"  < %s >\n", sComputeMode[device_prop.computeMode]);
}

/*Print out the device info assume CUDART >= 4000*/
/**
 * @brief Prints all information on hardware, meaning memory, 
 *        compute and features (no driver and runtime info).
 *        This assumes CUDART >= 4000.
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_hardware_info(cudaDeviceProp device_prop, input_gpu gpu_var) 
{
  lprintf("GPU_INIT", 10, "Device: %s\n", device_prop.name);

  print_memory_info(device_prop, gpu_var);
  print_compute_info(device_prop);
}

#endif
#endif
