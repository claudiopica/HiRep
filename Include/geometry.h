/***************************************************************************\
 * Copyright (c) 2008-2014, Claudio Pica                                    *
 * All rights reserved.                                                     *
 \**************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

/* this define the width of the borders for parallel dimensions
 * For different actions it must be modified
 */

#define BORDERSIZE 1

typedef struct _geometry_descriptor {
  int inner_master_pieces; /* number of inner pieces (1 o 2 for even odd or no_eo)) */
  int local_master_pieces;  // 
  int total_spinor_master_pieces;
  int total_gauge_master_pieces;
  int *master_start, *master_end; // Beginning has inner pieces, then local pieces, then the rest.
  int master_shift; /*  this is the odd spinor's shift, i.e. the index of the first odd entry in the full geometry */
  int ncopies_spinor;
  int ncopies_gauge;
  int *copy_from, *copy_to, *copy_len;
  int copy_shift; /*  this is the odd spinor's shift, i.e. the index of the first odd copy in the full geometry */
  int nbuffers_spinor;
  int nbuffers_gauge;
  int *rbuf_len, *sbuf_len;
  int *rbuf_from_proc, *rbuf_start;
  int *sbuf_to_proc, *sbuf_start;
  int gsize_spinor;
  int gsize_gauge;
} geometry_descriptor;


#include "hr_omp.h"

//Loop over pieces of given type
#define _PIECE_FOR(type,ip) \
_OMP_PRAGMA ( _omp_parallel )\
  for(int ip=0;\
      ip<(type)->local_master_pieces;\
      ip++ )

//Loop over sites of piece ip of given type
#define _SITE_FOR_RED(type,ip,is,redop1,redop2) \
_OMP_PRAGMA ( _omp_for redop1 redop2  )\
  for(int is=(type)->master_start[ip]; \
      is<=(type)->master_end[ip]; \
      is++ )

#define _SITE_FOR(type,ip,is) _SITE_FOR_RED(type,ip,is,,)
#define _SITE_FOR_SUM(type,ip,is,...) _SITE_FOR_RED(type,ip,is,_omp_sum(__VA_ARGS__),)
#define _SITE_FOR_MAX(type,ip,is,...) _SITE_FOR_RED(type,ip,is,_omp_max(__VA_ARGS__),)
#define _SITE_FOR_MIN(type,ip,is,...) _SITE_FOR_RED(type,ip,is,_omp_min(__VA_ARGS__),)


//Loop over sites of all pieces of given type
//do an openmp sum over the variables in the redop list
#define _MASTER_FOR_RED(type,is,redop1,redop2) \
  _PIECE_FOR((type),_master_for_ip_##is)    \
  _SITE_FOR_RED((type),_master_for_ip_##is,is,redop1,redop2)

#define _MASTER_FOR(type,is) _MASTER_FOR_RED(type,is,,)
#define _MASTER_FOR_SUM(type,is,...) _MASTER_FOR_RED(type,is,_omp_sum(__VA_ARGS__),)
#define _MASTER_FOR_MAX(type,is,...) _MASTER_FOR_RED(type,is,_omp_max(__VA_ARGS__),)
#define _MASTER_FOR_MIN(type,is,...) _MASTER_FOR_RED(type,is,_omp_min(__VA_ARGS__),)



int setup_process(int *argc, char ***argv);
int setup_replicas();
int finalize_process(void);

void origin_coord(int *c);
void other_proc_origin_coord(int * proc_coord, int *c);
void glb_to_proc(int *g, int *p);

int geometry_init(void);
void geometry_mpi_eo(void);
void geometry_mem_alloc(void);
int proc_up(int id, int dir);
int proc_dn(int id, int dir);

void init_geometry_SAP(void);
void test_geometry_mpi(void);
void test_geometry_mpi_eo(void);
void print_wdmatrix(char *filename);


#endif
