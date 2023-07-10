/***************************************************************************\
* Copyright (c) 2012, Ari Hietanen, Sofie Martins                           * 
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file init_gpu.c
 * @brief Functions initialize the GPU for the simulation and print collective 
 *        information on cluster drivers, necessary software and available hardware.
 *        To be logged to simulation outfiles.
 */

#ifdef WITH_GPU

#include "geometry.h"
#include "libhr_core.h"
#include "io.h"

//TODO: should these two be static?
//otherwise move to header file
void select_GPU(input_gpu);
int enable_GPU_peer_to_peer_access();

/**
 * @brief Call this in an init function to setup available graphics card for
 *        use. This also logs information on available software and hardware.
 * 
 * @param input_gpu             A struct containing information on the current active
 *                              GPU
 */
void init_gpu(input_gpu gpu_var_init) {
    lprintf("GPU_INIT", 0, "Initializing GPU\n");
    struct cudaDeviceProp device_prop;
    cudaGetDeviceProperties(&device_prop, gpu_var_init.gpuID);

    // Print GPU info
    print_device_count_info(gpu_var_init);
    print_software_info(device_prop);
    print_hardware_info(device_prop, gpu_var_init);
    print_supported_features(device_prop);

    // Select a card (no MPI) or bind cards to processes (MPI)
    select_GPU(gpu_var_init);

    // Setup global variables necessary for optimal kernel execution
    grid_size_max_gpu = device_prop.maxGridSize[0];
}

/**
 * @brief Selects GPU(s) depending on whether MPI is enabled or not.
 *
 * @param input_gpu             A struct containing information on the current active
 *                              GPU
 */
void select_GPU(input_gpu gpu_var_init) {
#ifndef WITH_MPI /* For Single GPU -> select device from input file */
    gpu_id = gpu_var_init.gpuID;
    cudaSetDevice(gpu_id);
    lprintf("GPU_INIT", 0, "Using GPU #%d\n", gpu_var_init.gpuID);
#else /* For Multi-GPU -> bind devices to local ranks using hwloc */
    /*
    TODO: The following code is wrong. What we need is
      * bind GPU to local rank instead of process ID
      * print info on local rank, process ID, node ID etc.
      * use hwloc to bind according to ideal CPU topology (SAM)
    */

    gpu_id = LID; // use local process id to select a GPU
    cudaSetDevice(gpu_id);
    int current_device;
    cudaGetDevice(&current_device);
    lprintf("GPU_INIT", 0, "GPU Affinity: GPU %d has been bound to MPI Rank %d (local %d)\n", current_device, PID, LID);
    enable_GPU_peer_to_peer_access();

#ifdef HWLOC
    find_physically_close_CPU_core();
#endif

#endif
}

/**
 * @brief Enables peer to peer access for GPUs connected to the same node.
 *
 * @param input_gpu           A struct containing information on the current active
 *                            GPU
 */
int enable_GPU_peer_to_peer_access() {
#if defined(WITH_MPI)
    int device_count = 0;
    cudaGetDeviceCount(&device_count);

    for (int i = 0; i < device_count; ++i) {
        if (i > PID) {
            int peer_access_available = 0;
            cudaDeviceCanAccessPeer(&peer_access_available, PID, i);
            lprintf("INFO", 0, "Peer-to-peer access: GPU Node %d can access node %d\n", PID, i);
            error(peer_access_available == 0, 1, "setup_GPU_peer_to_peer_access", "Unable to enable peer-to-peer access.\n");

            cudaDeviceEnablePeerAccess(PID, i);
            lprintf("INFO", 0, "Enabled peer-to-peer access from node %d to %d\n", PID, i);
        }
    }
#endif
    return 0;
}

#endif