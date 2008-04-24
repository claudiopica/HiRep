#ifndef GEOMETRY_H
#define GEOMETRY_H

typedef struct {
  unsigned int local_master_pieces, total_master_pieces;
  unsigned int *master_start, *master_end;
  unsigned int ncopies;
  unsigned int *copy_from, *copy_to, *copy_len;
  unsigned int nbuffers;
  unsigned int *buf_len;
  unsigned int *rbuf_from_proc, *rbuf_start;
  unsigned int *sbuf_to_proc, *sbuf_start;
  unsigned int gsize;
} geometry_descriptor;


typedef struct {
  int size;
  int local_pieces; /* reticolo ristretto + output buffers */
  int *start; /* starting index for each piece. len=tot_pieces */
  int *end; /* ending index for each piece. len=tot_pieces */
  int *len; /* length for each local piece for unique points. len=local_pieces */
} spinor_descriptor;




#define LOCAL_GD_FOR(gd,i,j) \
	for(i=0;i<(gd)->local_pieces;i++) for(j=(gd)->start[i];j<(gd)->start[i]+(gd)->len[i];j++)

void geometry_mpi(void);
void geometry_mpi_eo(void);
void geometry_set(void);
void geometry_mem_alloc(void);

void test_geometry_mpi(void);
void print_wdmatrix(char *filename);










/* THESE STUFF ARE SUPPOSED OT BE REMOVED IN A NEAR FUTURE */

#define np_x 3
#define np_y 3
#define np_z 3
#define np_t 3
#define myid 1 + np_t + np_x*np_t + np_y*np_x*np_t
#define myid_sign 0

/* THESE STUFF ARE SUPPOSED OT BE REMOVED IN A NEAR FUTURE  (UP TO HERE)*/


#endif
