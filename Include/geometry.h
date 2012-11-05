/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

/* this define the width of the borders for parallel dimensions
 * For different actions it must be modified
 */
#define BORDERSIZE 1


typedef struct _geometry_descriptor {
  int inner_master_pieces; /* number of inner pieces (1 o 2 for even odd or no_eo)) */
  int local_master_pieces;
  int total_spinor_master_pieces;
  int total_gauge_master_pieces;
  int *master_start, *master_end; 
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

#define _PIECE_FOR(type,i) \
  for(_PIECE_INDEX(i)=0;_PIECE_INDEX(i)<(type)->local_master_pieces;_PIECE_INDEX(i)++)

#define _SITE_FOR(type,i) \
  for(i=(type)->master_start[_PIECE_INDEX(i)]; \
      i<=(type)->master_end[_PIECE_INDEX(i)]; \
      i++ \
     )

#define _MASTER_FOR(type,i) \
  _PIECE_FOR((type),i)    \
_SITE_FOR((type),i)


int setup_process(int *argc, char ***argv);
int finalize_process(void);

void origin_coord(int *c);
void glb_to_proc(int *g, int *p);

int geometry_init(void);
void geometry_mpi_eo(void);
void geometry_mem_alloc(void);
int proc_up(int id, int dir);
int proc_dn(int id, int dir);

void test_geometry_mpi(void);
void test_geometry_mpi_eo(void);
void print_wdmatrix(char *filename);

#endif
