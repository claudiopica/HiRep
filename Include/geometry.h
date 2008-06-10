#ifndef GEOMETRY_H
#define GEOMETRY_H

typedef struct {
  unsigned int local_master_pieces, total_master_pieces;
  unsigned int *master_start, *master_end;
  /*
  * index in (master_start[0],master_end[0]) are all the inner master indexes
  * index in (master_start[i],master_end[i]) i=1 ... local_master_pieces-1
    are all the border master indexes
  * index in (master_start[i],master_end[i]) i=local_master_pieces ... total_master_pieces-1
    are all the buffer master indexes
  */
  unsigned int ncopies;
  unsigned int *copy_from, *copy_to, *copy_len;
  /*
  the master indexes in (copy_from[i],copy_from[i]+copy_len[i]-1)
  must be copied in the copy indexes in (copy_to[i],copy_to[i]+copy_len[i]-1)
  i=0 ... ncopies-1
  */
  unsigned int nbuffers;
  unsigned int *rbuf_len, *sbuf_len;
  unsigned int *rbuf_from_proc, *rbuf_start;
  unsigned int *sbuf_to_proc, *sbuf_start;
  /*
  The i-th trasmission process (i=1 ... nbuffers-1) controls a flow of
  data from the processor A to the processor B.
  
  The processor A stores these data in the indexes
  (sbuf_start[i],sbuf_start[i]+sbuf_len[i]-1). It must send these data
  through a channel identified by i to the processor sbuf_to_proc[i]=B.
  
  The processor B receive the data through the channel identified by i,
  from the processor rbuf_from_proc[i]=A, and stores them in the
  indexes (rbuf_start[i],rbuf_start[i]+rbuf_len[i]-1).
  */
  unsigned int gsize; /*the number of all the indexes in the global geometry*/
} geometry_descriptor;

#define _PIECE_INDEX(i) i##_pindex
#define _DECLARE_INT_ITERATOR(i) int i, _PIECE_INDEX(i)
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
