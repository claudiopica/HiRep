/***************************************************************************\
* Copyright (c) 2012, Ari Hietanen                                          * 
* All rights reserved.                                                      * 
\***************************************************************************/
 
/**
 * @file init_gpu.c
 * @brief Functions initialize the GPU for the simulation and print collective 
 *        information on cluster drivers, necessary software and available hardware.
 *        To be logged to simulation outfiles.
 */

#ifdef WITH_GPU
extern "C" {
  #include "logger.h"
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gpu.h"
#include "io.h"
#include "global.h"
#include "utils.h"
#include "error.h"
#include "geometry.h"
#include "gpu_info.h"
#ifdef WITH_MPI
  #include "mpi.h"
#endif

void select_GPU(input_gpu);
int enable_GPU_peer_to_peer_access();

/**
 * @brief Call this in an init function to setup available graphics card for
 *        use. This also logs information on available software and hardware.
 * 
 * @param input_gpu             A struct containing information on the current active
 *                              GPU
 */ 
void init_gpu(input_gpu gpu_var)
{
  lprintf("GPU_INIT", 0, "Initializing GPU\n");
  struct cudaDeviceProp device_prop;
  CHECK_CUDA(cudaGetDeviceProperties(&device_prop, gpu_var.gpuID));

  // Print GPU info
  print_device_count_info(gpu_var);
  print_driver_info(device_prop);
  print_runtime_info(device_prop);
  print_hardware_info(device_prop, gpu_var);

  // Select a card (no MPI) or bind cards to processes (MPI)
  //select_GPU(gpu_var);

  // Setup global variables necessary for optimal kernel execution
  grid_size_max_gpu = device_prop.maxGridSize[0];
}

/**
 * @brief Selects GPU(s) depending on whether MPI is enabled or not.
 *
 * @param input_gpu             A struct containing information on the current active
 *                              GPU
 */
void select_GPU(input_gpu gpu_var) 
{
  #ifndef WITH_MPI /* For Single GPU -> select device with ID=0 */
    CHECK_CUDA(cudaSetDevice(gpu_var.gpuID));
    lprintf("GPU_INIT", 0, "Using GPU #%d\n", gpu_var.gpuID);
  #else /* For Multi-GPU -> bind devices to local ranks using hwloc */
    /*
    TODO: The following code is wrong. What we need is
      * bind GPU to local rank instead of process ID
      * print info on local rank, process ID, node ID etc.
      * use hwloc to bind according to ideal CPU topology (SAM)
    */
    cudaSetDevice(PID);
    int current_device;
    CHECK_CUDA(cudaGetDevice(&current_device));
    lprintf("GPU_INIT", 0, "GPU Affinity: GPU Node %d has been bound to MPI Thread of Rank %d\n", current_device, PID);
    //enable_GPU_peer_to_peer_access();
  #endif
}

/**
 * @brief Enables peer to peer access for GPUs connected to the same node.
 *
 * @param input_gpu           A struct containing information on the current active
 *                            GPU
 */
int enable_GPU_peer_to_peer_access() 
{
  // TODO: For more than one node we need to use local MPI ranks instead of PIDs (SAM)
  #if defined(WITH_MPI) 
    int device_count = 0;
    CHECK_CUDA(cudaGetDeviceCount(&device_count));

    for (int i = 0; i < device_count; ++i) 
    {
      if (i > PID) 
      {
        int peer_access_available = 0;
        CHECK_CUDA(cudaDeviceCanAccessPeer(&peer_access_available, PID, i));
        lprintf("INFO", 0, "Peer-to-peer access: GPU Node %d can access node %d\n", PID, i);
        error(peer_access_available == 0, 1, "setup_GPU_peer_to_peer_access", "Unable to enable peer-to-peer access.\n");

        CHECK_CUDA(cudaDeviceEnablePeerAccess(PID, i));
        lprintf("INFO", 0, "Enabled peer-to-peer access from node %d to %d\n", PID, i);
      }
    }
  #endif
  return 0;
}

#endif