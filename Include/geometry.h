#ifndef GEOMETRY_H
#define GEOMETRY_H

typedef struct {
  unsigned int local_pieces; /* reticolo ristretto + output buffers */
  unsigned int tot_pieces; /* above + input buffers */
  unsigned int *start; /* starting index for each piece. len=tot_pieces */
  unsigned int *len; /* length for each local piece for unique points. len=local_pieces */
  unsigned int *tlen; /* total lenght for each piece. len=tot_pieces */
  unsigned int ncopy; /* tot number of buffer copies needed to fill the holes. */
  unsigned int *copy_from, *copy_to, *copy_len; /* len=ncopy */
  int *sid; /* cartesian id of the node to which the local buffers must be sent. len=local_pieces*/
  int *rid; /* cartesian id of the node from which the input buffers are received. len=tot_pieces-local_pieces*/
  unsigned int gsize; /* global size of the associated array = \sum_i tlen_i */
} geometry_descriptor;

typedef struct {
  int size;
  int local_pieces; /* reticolo ristretto + output buffers */
  int *start; /* starting index for each piece. len=tot_pieces */
  int *end; /* ending index for each piece. len=tot_pieces */
  int *len; /* length for each local piece for unique points. len=local_pieces */
} spinor_descriptor;

#define LOCAL_SD_FOR(sd,i,j) \
	for(i=0;i<sd->local_pieces;i++) for(j=sd->start[i];j<sd->start[i]+sd->len[i];j++)

void geometry_mpi(void);
void geometry_init(void);
void geometry_set(void);
void geometry_mem_alloc(geometry_descriptor*);
void geometry_blocked(void);
void geometry_blocked_noT(void);
void geometry_lexi(void);
void geometry_eo_lexi(void);

void test_geometry(void);
void print_geometry(void);
void print_wdmatrix(char *filename);

#endif
