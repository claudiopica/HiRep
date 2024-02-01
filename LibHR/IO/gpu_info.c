/***************************************************************************\
* Copyright (c) 2012, 2022, Ari Hietanen, Sofie Martins                     * 
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file gpu_info.c
 * @brief Prints information on available hardware and software on the 
 *        cluster relating to GPUs.
 */
#ifdef WITH_GPU

#include "io.h"
#include "error.h"
#include "gpu.h"
#ifdef WITH_MPI
#include "mpi.h"
//#include <mpi-ext.h> // Needed for CUDA-awareness check, see https://docs.open-mpi.org/en/v5.0.x/man-openmpi/man3/MPIX_Query_cuda_support.3.html
#endif

const char *sComputeMode[] = {
    "Default (multiple host threads can use ::cudaSetDevice() with device simultaneously)",
    "Exclusive (only one host thread in one process is able to use ::cudaSetDevice() with this device)",
    "Prohibited (no host thread can use ::cudaSetDevice() with this device)",
    "Exclusive Process (many threads in one process is able to use ::cudaSetDevice() with this device)",
    "Unknown",
    NULL
};

#ifndef HIP /* ----- CUDA ----- */
/**
 * @brief Query number of GPUs and print related information
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_device_count_info(input_gpu gpu_var_init) {
    int device_count;
    cudaGetDeviceCount(&device_count);
    lprintf("GPU_INIT", 0, "GPU_ID = %d\n", gpu_var_init.gpuID);
    //error(gpu_id > device_count, 1, "init_gpu", "Illegal device ID");
    // I don't see what we need this for (SAM)
}

/**
 * @brief Print CUDA driver and runtime version information
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_software_info(cudaDeviceProp device_prop) {
    int driver_version, runtime_version;
    cudaDriverGetVersion(&driver_version);
    cudaRuntimeGetVersion(&runtime_version);
    lprintf("GPU_INIT", 10, "CUDA Capability Major/Minor version number: %d.%d\n", device_prop.major, device_prop.minor);
    lprintf("GPU_INIT", 10, "CUDA Driver Version / Runtime Version: %d.%d / %d.%d\n", driver_version / 1000,
            driver_version % 100, runtime_version / 1000, runtime_version % 100);
}

/**
 * @brief Print Global memory information including bandwidth 
 *        paramters and supported features
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_global_memory_info(cudaDeviceProp device_prop, input_gpu gpu_var_init) {
    int mem_bus_width;
    int mem_clock;

    // Query properties
    cuDeviceGetAttribute(&mem_bus_width, CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH, gpu_var_init.gpuID);
    cuDeviceGetAttribute(&mem_clock, CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE, gpu_var_init.gpuID);

    // Print formatted
    lprintf("GPU_INIT", 10, "Total amount of global memory: %.0f MB (%lluB)\n",
            (float)device_prop.totalGlobalMem / ((float)(1 << 20)), (unsigned long long)device_prop.totalGlobalMem);
    lprintf("GPU_INIT", 10, "Memory Clock rate: %.3f Mhz\n", mem_clock * 1e-3f);
    lprintf("GPU_INIT", 10, "Memory Bus Width: %d-bit\n", mem_bus_width);
    lprintf("GPU_INIT", 10, "Maximum memory pitch: %uB\n", device_prop.memPitch);
    lprintf("GPU_INIT", 10, "Integrated GPU sharing Host Memory: %s\n", device_prop.integrated ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Support host page-locked memory mapping: %s\n", device_prop.canMapHostMemory ? "Yes" : "No");
}

/**
 * @brief Print information on shared memory
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_shared_memory_info(cudaDeviceProp device_prop) {
    lprintf("GPU_INIT", 10, "Total amount of shared memory per block: %dB\n", device_prop.sharedMemPerBlock);
}

/**
 * @brief Print (L2) cache information
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_cache_info(cudaDeviceProp device_prop, input_gpu gpu_var_init) {
    int l2_cache_size;

    // Query properties
    cuDeviceGetAttribute(&l2_cache_size, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, gpu_var_init.gpuID);

    // Print formatted
    lprintf("GPU_INIT", 10, "L2 Cache Size: %dB\n", l2_cache_size);
    // TODO: L1 cache??
}

/**
 * @brief Print information on constant memory, in particular
 *        amount available, alignment, texture and layered texture 
 *        memory paramters
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_constant_memory_info(cudaDeviceProp device_prop) {
    lprintf("GPU_INIT", 10, "Max Texture dimension size (x,y,z): 1D=(%d), 2D=(%d,%d), 3D=(%d,%d,%d)\n",
            device_prop.maxTexture1D, device_prop.maxTexture2D[0], device_prop.maxTexture2D[1], device_prop.maxTexture3D[0],
            device_prop.maxTexture3D[1], device_prop.maxTexture3D[2]);

    lprintf("GPU_INIT", 10, "Max Layered Texture dimension size (dim) x layers: 1D=(%d) x %d, 2D=(%d,%d) x %d\n",
            device_prop.maxTexture1DLayered[0], device_prop.maxTexture1DLayered[1], device_prop.maxTexture2DLayered[0],
            device_prop.maxTexture2DLayered[1], device_prop.maxTexture2DLayered[2]);

    lprintf("GPU_INIT", 10, "Texture alignment: %uB\n", device_prop.textureAlignment);
    lprintf("GPU_INIT", 10, "Total amount of constant memory: %dB\n", device_prop.totalConstMem);
    lprintf("GPU_INIT", 10, "Alignment requirement for Surfaces: %s\n", device_prop.surfaceAlignment ? "Yes" : "No");
}

/**
 * @brief Prints all memory related info, here on global, shared
 *        memory, cache and constant memory, see resp. functions
 *        for more information.
 * 
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_memory_info(cudaDeviceProp device_prop, input_gpu gpu_var_init) {
    print_global_memory_info(device_prop, gpu_var_init);
    print_shared_memory_info(device_prop);
    print_cache_info(device_prop, gpu_var_init);
    print_constant_memory_info(device_prop);
}

/**
 * @brief Prints information on compute parameters
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_compute_info(cudaDeviceProp device_prop, input_gpu gpu_var_init) {
    // Print formatted
    lprintf("GPU_INIT", 10, "Multiprocessors: %d\n", device_prop.multiProcessorCount);
    lprintf("GPU_INIT", 10, "GPU Clock Speed: %.2f GHz\n", device_prop.clockRate * 1e-6f);
    lprintf("GPU_INIT", 10, "Total number of register per block: %dB\n", device_prop.regsPerBlock);
    lprintf("GPU_INIT", 10, "Warp size: %dB\n", device_prop.warpSize);
    error(device_prop.warpSize != 32, 1, "init_gpu", "Error: warp size 32 assumed in global sum\n");
    lprintf("GPU_INIT", 10, "Maximum number of threads per block: %d\n", device_prop.maxThreadsPerBlock);
    lprintf("GPU_INIT", 10, "Maximum size of each dimension of a block (x,y,z): (%d,%d,%d)\n", device_prop.maxThreadsDim[0],
            device_prop.maxThreadsDim[1], device_prop.maxThreadsDim[2]);
    lprintf("GPU_INIT", 10, "Maximum size of each dimension of a grid (x,y,z): (%d,%d,%d)\n", device_prop.maxGridSize[0],
            device_prop.maxGridSize[1], device_prop.maxGridSize[2]);
    lprintf("GPU_INIT", 10, "Concurrent copy and execution: %s with %d copy engine(s)\n",
            device_prop.deviceOverlap ? "Yes" : "No", device_prop.asyncEngineCount);
    lprintf("GPU_INIT", 10, "Run time limit on kernels: %s\n", device_prop.kernelExecTimeoutEnabled ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Concurrent kernel execution: %s\n", device_prop.concurrentKernels ? "Yes" : "No");
}

/**
 * @brief Checks for a number of other supported features and
 *        prints information on them.
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_supported_features(cudaDeviceProp device_prop) {
    lprintf("GPU_INIT", 10, "Device has ECC support enabled:                %s\n", device_prop.ECCEnabled ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Device is using TCC driver mode:               %s\n", device_prop.tccDriver ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Device supports Unified Addressing (UVA):      %s\n",
            device_prop.unifiedAddressing ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Device PCI Bus ID / PCI location ID:           %d / %d\n", device_prop.pciBusID,
            device_prop.pciDeviceID);
    lprintf("GPU_INIT", 10, "Compute Mode:\n");
    lprintf("GPU_INIT", 10, "  < %s >\n", sComputeMode[device_prop.computeMode]);

// Multi-GPU calculations are not supported for the old geometry
#if defined(WITH_GPU) && defined(WITH_MPI) && !defined(WITH_NEW_GEOMETRY)
    error(
        1, 1, __func__,
        "Legacy geometry not supported for Multi-GPU compilation. Please enable WITH_NEW_GEOMETRY in compilation flags. Exiting. ");
#endif

#ifdef WITH_MPI
    int cuda_aware_support = 0;
// TODO: This possibly only works for OpenMPI (SAM)
#if defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
    cuda_aware_support = MPIX_Query_cuda_support();
#else
#error "Your MPI was not installed with CUDA support. This is unsupported in HiRep.\n"
#endif

    if (cuda_aware_support) {
        lprintf("GPU_INIT", 10, "MPI implementation CUDA-aware? yes.\n");
    } else {
        error(1, 1, __func__, "MPI implementation CUDA-aware? no. Exiting. \n");
    }
#endif
}

/**
 * @brief Prints peak performance metric estimates. 
 *        This allows to check, whether we are fully utilizing the 
 *        capabilities of the hardware
 * 
 * 
*/
void print_performance_metrics() {
    /**
   * Code snippet based on:
   * https://developer.nvidia.com/blog/how-query-device-properties-and-handle-errors-cuda-cc/
   * Access 2023-02-27
   * 
  */
    int n_devices;
    cudaGetDeviceCount(&n_devices);
    double peak_memory_bandwidth = 0;
    for (int i = 0; i < n_devices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        peak_memory_bandwidth += 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6;
    }
    lprintf("GPU_INIT", 10, "Peak Memory Bandwidth (GB/s): %1.6g\n", peak_memory_bandwidth);
}

/*Print out the device info assume CUDART >= 4000*/
/**
 * @brief Prints all information on hardware, meaning memory, 
 *        compute and features (no driver and runtime info).
 *        This assumes CUDART >= 4000.
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_hardware_info(cudaDeviceProp device_prop, input_gpu gpu_var_init) {
    lprintf("GPU_INIT", 10, "Device: %s\n", device_prop.name);

    print_memory_info(device_prop, gpu_var_init);
    print_compute_info(device_prop, gpu_var_init);
    print_performance_metrics();
}

#else /* ----- HIP -----*/

/**
 * @brief Query number of GPUs and print related information
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_device_count_info(input_gpu gpu_var_init) {
    int device_count;
    hipGetDeviceCount(&device_count);
    lprintf("GPU_INIT", 0, "GPU_ID = %d\n", gpu_var_init.gpuID);
    //error(gpu_id > device_count, 1, "init_gpu", "Illegal device ID");
    // I don't see what we need this for (SAM)
}

/**
 * @brief Print CUDA driver and runtime version information
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_software_info(hipDeviceProp_t device_prop) {
    int driver_version, runtime_version;
    hipDriverGetVersion(&driver_version);
    hipRuntimeGetVersion(&runtime_version);
    lprintf("GPU_INIT", 10, "CUDA Capability Major/Minor version number: %d.%d\n", device_prop.major, device_prop.minor);
    lprintf("GPU_INIT", 10, "CUDA Driver Version / Runtime Version: %d.%d / %d.%d\n", driver_version / 1000,
            driver_version % 100, runtime_version / 1000, runtime_version % 100);
}

/**
 * @brief Print Global memory information including bandwidth 
 *        paramters and supported features
 *
 * @param hipDeviceProp_t        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_global_memory_info(hipDeviceProp_t device_prop, input_gpu gpu_var_init) {
    int mem_bus_width;
    int mem_clock;

    // Query properties
    hipDeviceGetAttribute(&mem_bus_width, hipDeviceAttributeMemoryBusWidth, gpu_var_init.gpuID);
    hipDeviceGetAttribute(&mem_clock, hipDeviceAttributeMemoryClockRate, gpu_var_init.gpuID);

    // Print formatted
    lprintf("GPU_INIT", 10, "Total amount of global memory: %.0f MB (%lluB)\n",
            (float)device_prop.totalGlobalMem / ((float)(1 << 20)), (unsigned long long)device_prop.totalGlobalMem);
    lprintf("GPU_INIT", 10, "Memory Clock rate: %.3f Mhz\n", mem_clock * 1e-3f);
    lprintf("GPU_INIT", 10, "Memory Bus Width: %d-bit\n", mem_bus_width);
    lprintf("GPU_INIT", 10, "Maximum memory pitch: %uB\n", device_prop.memPitch);
    lprintf("GPU_INIT", 10, "Integrated GPU sharing Host Memory: %s\n", device_prop.integrated ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Support host page-locked memory mapping: %s\n", device_prop.canMapHostMemory ? "Yes" : "No");
}

/**
 * @brief Print information on shared memory
 *
 * @param hipDeviceProp_t        A CUDA class containing information on the device.
 */
void print_shared_memory_info(hipDeviceProp_t device_prop) {
    lprintf("GPU_INIT", 10, "Total amount of shared memory per block: %dB\n", device_prop.sharedMemPerBlock);
}

/**
 * @brief Print (L2) cache information
 *
 * @param hipDeviceProp_t        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_cache_info(hipDeviceProp_t device_prop, input_gpu gpu_var_init) {
    int l2_cache_size;

    // Query properties
    hipDeviceGetAttribute(&l2_cache_size, hipDeviceAttributeL2CacheSize, gpu_var_init.gpuID);

    // Print formatted
    lprintf("GPU_INIT", 10, "L2 Cache Size: %dB\n", l2_cache_size);
    // TODO: L1 cache??
}

/**
 * @brief Print information on constant memory, in particular
 *        amount available, alignment, texture and layered texture 
 *        memory paramters
 *
 * @param hipDeviceProp_t        A CUDA class containing information on the device.
 */
void print_constant_memory_info(hipDeviceProp_t device_prop) {
#if 0
    lprintf("GPU_INIT", 10, "Max Texture dimension size (x,y,z): 1D=(%d), 2D=(%d,%d), 3D=(%d,%d,%d)\n",
            device_prop.maxTexture1D, device_prop.maxTexture2D[0], device_prop.maxTexture2D[1], device_prop.maxTexture3D[0],
            device_prop.maxTexture3D[1], device_prop.maxTexture3D[2]);

    lprintf("GPU_INIT", 10, "Max Layered Texture dimension size (dim) x layers: 1D=(%d) x %d, 2D=(%d,%d) x %d\n",
            device_prop.maxTexture1DLayered[0], device_prop.maxTexture1DLayered[1], device_prop.maxTexture2DLayered[0],
            device_prop.maxTexture2DLayered[1], device_prop.maxTexture2DLayered[2]);

    lprintf("GPU_INIT", 10, "Texture alignment: %uB\n", device_prop.textureAlignment);
    lprintf("GPU_INIT", 10, "Total amount of constant memory: %dB\n", device_prop.totalConstMem);
    lprintf("GPU_INIT", 10, "Alignment requirement for Surfaces: %s\n", device_prop.surfaceAlignment ? "Yes" : "No");
#endif
}

/**
 * @brief Prints all memory related info, here on global, shared
 *        memory, cache and constant memory, see resp. functions
 *        for more information.
 * 
 * @param hipDeviceProp_t        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_memory_info(hipDeviceProp_t device_prop, input_gpu gpu_var_init) {
    print_global_memory_info(device_prop, gpu_var_init);
    print_shared_memory_info(device_prop);
    print_cache_info(device_prop, gpu_var_init);
    print_constant_memory_info(device_prop);
}

/**
 * @brief Prints information on compute parameters
 *
 * @param hipDeviceProp_t        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_compute_info(hipDeviceProp_t device_prop, input_gpu gpu_var_init) {
#if 0
    // Print formatted
    lprintf("GPU_INIT", 10, "Multiprocessors: %d\n", device_prop.multiProcessorCount);
    lprintf("GPU_INIT", 10, "GPU Clock Speed: %.2f GHz\n", device_prop.clockRate * 1e-6f);
    lprintf("GPU_INIT", 10, "Total number of register per block: %dB\n", device_prop.regsPerBlock);
    lprintf("GPU_INIT", 10, "Warp size: %dB\n", device_prop.warpSize);
    error(device_prop.warpSize != 32, 1, "init_gpu", "Error: warp size 32 assumed in global sum\n");
    lprintf("GPU_INIT", 10, "Maximum number of threads per block: %d\n", device_prop.maxThreadsPerBlock);
    lprintf("GPU_INIT", 10, "Maximum size of each dimension of a block (x,y,z): (%d,%d,%d)\n", device_prop.maxThreadsDim[0],
            device_prop.maxThreadsDim[1], device_prop.maxThreadsDim[2]);
    lprintf("GPU_INIT", 10, "Maximum size of each dimension of a grid (x,y,z): (%d,%d,%d)\n", device_prop.maxGridSize[0],
            device_prop.maxGridSize[1], device_prop.maxGridSize[2]);
    lprintf("GPU_INIT", 10, "Concurrent copy and execution: %s with %d copy engine(s)\n",
            device_prop.deviceOverlap ? "Yes" : "No", device_prop.asyncEngineCount);
    lprintf("GPU_INIT", 10, "Run time limit on kernels: %s\n", device_prop.kernelExecTimeoutEnabled ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Concurrent kernel execution: %s\n", device_prop.concurrentKernels ? "Yes" : "No");
#endif
}

/**
 * @brief Checks for a number of other supported features and
 *        prints information on them.
 *
 * @param hipDeviceProp_t        A CUDA class containing information on the device.
 */
void print_supported_features(hipDeviceProp_t device_prop) {
#if 0
    lprintf("GPU_INIT", 10, "Device has ECC support enabled:                %s\n", device_prop.ECCEnabled ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Device is using TCC driver mode:               %s\n", device_prop.tccDriver ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Device supports Unified Addressing (UVA):      %s\n",
            device_prop.unifiedAddressing ? "Yes" : "No");
    lprintf("GPU_INIT", 10, "Device PCI Bus ID / PCI location ID:           %d / %d\n", device_prop.pciBusID,
            device_prop.pciDeviceID);
    lprintf("GPU_INIT", 10, "Compute Mode:\n");
    lprintf("GPU_INIT", 10, "  < %s >\n", sComputeMode[device_prop.computeMode]);

// Multi-GPU calculations are not supported for the old geometry
#if defined(WITH_GPU) && defined(WITH_MPI) && !defined(WITH_NEW_GEOMETRY)
    error(
        1, 1, __func__,
        "Legacy geometry not supported for Multi-GPU compilation. Please enable WITH_NEW_GEOMETRY in compilation flags. Exiting. ");
#endif

#ifdef WITH_MPI
    int cuda_aware_support = 0;
// TODO: This possibly only works for OpenMPI (SAM)
//#if defined(OMPI_HAVE_MPI_EXT_CUDA) && OMPI_HAVE_MPI_EXT_CUDA
//    cuda_aware_support = MPIX_Query_cuda_support();
//#else
//#error "Your MPI was not installed with CUDA support. This is unsupported in HiRep.\n"
#endif

    //if (cuda_aware_support) {
    //    lprintf("GPU_INIT", 10, "MPI implementation CUDA-aware? yes.\n");
    //} else {
    //    error(1, 1, __func__, "MPI implementation CUDA-aware? no. Exiting. \n");
   // }
#endif
}

/**
 * @brief Prints peak performance metric estimates. 
 *        This allows to check, whether we are fully utilizing the 
 *        capabilities of the hardware
 * 
 * 
*/
void print_performance_metrics() {
    /**
   * Code snippet based on:
   * https://developer.nvidia.com/blog/how-query-device-properties-and-handle-errors-cuda-cc/
   * Access 2023-02-27
   * 
  */
    int n_devices;
    hipGetDeviceCount(&n_devices);
    double peak_memory_bandwidth = 0;
    for (int i = 0; i < n_devices; i++) {
        hipDeviceProp_t prop;
        hipGetDeviceProperties(&prop, i);
        peak_memory_bandwidth += 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6;
    }
    lprintf("GPU_INIT", 10, "Peak Memory Bandwidth (GB/s): %1.6g\n", peak_memory_bandwidth);
}

/*Print out the device info assume CUDART >= 4000*/
/**
 * @brief Prints all information on hardware, meaning memory, 
 *        compute and features (no driver and runtime info).
 *        This assumes CUDART >= 4000.
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_hardware_info(hipDeviceProp_t device_prop, input_gpu gpu_var_init) {
    lprintf("GPU_INIT", 10, "Device: %s\n", device_prop.name);

    print_memory_info(device_prop, gpu_var_init);
    print_compute_info(device_prop, gpu_var_init);
    print_performance_metrics();
}

#endif
#endif
