/**
 * @file geometry_init.h
 * @brief Initialization functions, that determine all important parameters of the 
 *        geometry, so that communications and operations can be completed correctly.
 */
/// Headerfile for:
/// - geometry_init.c
/// - geometry_SAP.c
/// - print_wdmatrix.c
/// - geometry_mpi_eo.c

#ifndef GEOMETRY_INIT_H
#define GEOMETRY_INIT_H

#include "new_geometry.h"

#ifdef __cplusplus
   extern "C" {
#endif

// geometry_init.c
void origin_coord(int *c);
void other_proc_origin_coord(int *proc_coord, int *c);
void glb_to_proc(int *g, int *p);
int proc_up(int id, int dir);
int proc_dn(int id, int dir);
int proc_id(int coords[4]);
int geometry_init(void);
void print_gd(geometry_descriptor *gd);

// geometry_SAP.c
void init_geometry_SAP(void);
void empty_buffers(spinor_field *s); //in SAP

//print_wdmatrix.c
#ifdef WITH_UNTESTED
void print_wdmatrix(char *filename);
#endif

// geometry_mpi_eo.c
void geometry_mpi_eo(void);
void geometry_mem_alloc(void);
void test_geometry_mpi(void);
void test_geometry_mpi_eo(void);
void free_geometry_mpi_eo(void);

#ifdef __cplusplus
   }
#endif
#endif


