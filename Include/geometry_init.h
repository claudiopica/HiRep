/**
 * @file geometry_init.h
 * @brief Initialization functions, that determine all important parameters of the 
 *        geometry, so that communications and operations can be completed correctly.
 */

#ifndef GEOMETRY_INIT_H
#define GEOMETRY_INIT_H
#ifdef __cplusplus
   extern "C" {
#endif

void origin_coord(int *c);
void other_proc_origin_coord(int *proc_coord, int *c);
void glb_to_proc(int *g, int *p);
int proc_up(int id, int dir);
int proc_dn(int id, int dir);
int proc_id(int coords[4]);

int geometry_init(void);
void init_geometry_SAP(void);
void print_wdmatrix(char *filename);


// needed when WITH_NEW_GEOMETRY
#include <stddef.h>
void define_geometry();
void* sendbuf_alloc(size_t bytes_per_site);
void sync_field(geometry_descriptor *gd, int byte_per_site, int is_spinor_like, void *latticebuf, void *sb_ptr);
int test_define_geometry();
void sendbuf_report();

// needed for old gemetry
void geometry_mpi_eo(void);
void geometry_mem_alloc(void);
void test_geometry_mpi(void);
void test_geometry_mpi_eo(void);
void free_geometry_mpi_eo(void);

#ifdef WITH_GPU
#include "gpu.h"
   /**
    * @brief Call this in an init function to setup available graphics card for
    *        use. This also logs information on available software and hardware.
    * 
    * @param input_gpu             A struct containing information on the current active
    *                              GPU
    */
   void init_gpu(input_gpu gpu_var);
   void init_neighbors_gpu();
#endif

#ifdef __cplusplus
   }
#endif
#endif


