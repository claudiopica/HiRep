/**
 * @file geometry_init.h
 * @brief Initialization functions, that determine all important parameters of the 
 *        geometry, so that communications and operations can be completed correctly.
 */

#ifndef GEOMETRY_INIT_H
#define GEOMETRY_INIT_H

void origin_coord(int *c);
void other_proc_origin_coord(int *proc_coord, int *c);
void glb_to_proc(int *g, int *p);

int geometry_init(void);
void geometry_mpi_eo(void);
void geometry_mem_alloc(void);
int proc_up(int id, int dir);
int proc_dn(int id, int dir);

void free_geometry_mpi_eo(void);

void init_geometry_SAP(void);
void test_geometry_mpi(void);
void test_geometry_mpi_eo(void);
void print_wdmatrix(char *filename);

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

#endif


