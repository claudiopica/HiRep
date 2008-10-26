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
  unsigned int local_master_pieces, total_master_pieces;
  unsigned int *master_start, *master_end;
  unsigned int ncopies;
  unsigned int *copy_from, *copy_to, *copy_len;
  unsigned int nbuffers;
  unsigned int *rbuf_len, *sbuf_len;
  unsigned int *rbuf_from_proc, *rbuf_start;
  unsigned int *sbuf_to_proc, *sbuf_start;
  unsigned int gsize;
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

void test_geometry_mpi(void);
void test_geometry_mpi_eo(void);
void print_wdmatrix(char *filename);

#endif
