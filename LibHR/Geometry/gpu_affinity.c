#ifdef WITH_MPI

#ifdef WITH_GPU

#include <mpi.h>
#include "global.h"
#include "io.h"
#include "Geometry/gpu_affinity.h"

/**
 * TODO:
 * =====
 * 
 * * Bind GPUs for multi-node calculation with hwloc
 * * Print detailed information on communication setup, warnings for possible performance problems in the hardware support
 * 
*/

#ifdef HWLOC
static hwloc_topology_t topology = NULL;

hwloc_topology_t init_topology() {
    hwloc_topology_t topology;
    hwloc_topology_init(&topology);
    hwloc_topology_set_flags(topology, HWLOC_TOPOLOGY_FLAG_IS_THISSYSTEM);
    hwloc_topology_load(topology);
    return topology;
}

hwloc_bitmap_t physically_close_cpuset(hwloc_topology_t top) {
    hwloc_bitmap_t cpuset = hwloc_bitmap_alloc();
    hwloc_cudart_get_device_cpuset(top, gpu_id, cpuset);
    return cpuset;
}

int find_physically_close_CPU_core() {
    // Adapted from professional CUDA C programming
    topology = init_topology();
    hwloc_bitmap_t cpuset = physically_close_cpuset(topology);

    hwloc_cpuset_t close_cpus = (hwloc_cpuset_t)malloc(sizeof(hwloc_cpuset_t));
    hwloc_set_cpubind(topology, close_cpus, HWLOC_CPUBIND_PROCESS);

    int match = 0;
    int i;
    int local_rank;
    int cpu = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);

    // Iterate over physically close cores
    // We need to make sure that we don't assign the same
    // CPU to two MPI processes because they both happen to
    // be close to a GPU (and first in the array or so)
    hwloc_bitmap_foreach_begin(i, cpuset) if (match == local_rank) {
        cpu = i;
        break;
    }
    ++match;
    hwloc_bitmap_foreach_end();

    // Single CPU we want to bind the process to
    // allocate and set to the bitmap of the cpu we
    // just found
    hwloc_bitmap_t onecpu = hwloc_bitmap_alloc();
    hwloc_bitmap_set(onecpu, cpu);

    // Bind to the correct process
    hwloc_set_cpubind(topology, onecpu, 0);

    // Free constructed
    hwloc_bitmap_free(onecpu);
    hwloc_bitmap_free(cpuset);
    hwloc_topology_destroy(topology);

    lprintf("HWLOC", 10, "Bound GPU %d at MPI local rank %d to physically close CPU %d.\n", PID, PID, cpu);
    return cpu;
}

#endif
#endif
#endif
