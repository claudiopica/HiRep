#ifndef GEOMETRY_H
#define GEOMETRY_H

typedef struct {
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



#define _MASTER_FOR(type,i) \
	for(_PIECE_INDEX(i)=0;_PIECE_INDEX(i)<(type)->local_master_pieces;_PIECE_INDEX(i)++) \
	for(i=(type)->master_start[_PIECE_INDEX(i)]; \
	    i<=(type)->master_end[_PIECE_INDEX(i)]; \
	    i++ \
	   )

void geometry_mpi(void);
void geometry_mpi_eo(void);
void geometry_set(void);
void geometry_mem_alloc(void);

void test_geometry_mpi(void);
void test_geometry_mpi_eo(void);
void print_wdmatrix(char *filename);










/* THESE STUFF ARE SUPPOSED OT BE REMOVED IN A NEAR FUTURE */

#define SAFE_MOD(x,y) ( (x)>=0 ? (x)%(y) : ((y)-(abs(x)%(y)))%(y) )


#define np_x 3
#define np_y 3
#define np_z 3
#define np_t 3
#define myid 42
#define myid_sign 1
#define proc_t(id) ((id)%np_t)
#define proc_x(id) (((id)/np_t)%np_x)
#define proc_y(id) (((id)/(np_t*np_x))%np_y)
#define proc_z(id) ((id)/(np_t*np_x*np_y))
#define global_t(id,lt) SAFE_MOD((proc_t(id)*T+(lt)),GLOBAL_T)
#define global_x(id,lx) SAFE_MOD((proc_x(id)*X+(lx)),GLOBAL_X)
#define global_y(id,ly) SAFE_MOD((proc_y(id)*Y+(ly)),GLOBAL_Y)
#define global_z(id,lz) SAFE_MOD((proc_z(id)*Z+(lz)),GLOBAL_Z)

/* THESE STUFF ARE SUPPOSED OT BE REMOVED IN A NEAR FUTURE  (UP TO HERE)*/


#endif
