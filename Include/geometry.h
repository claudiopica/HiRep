#ifndef GEOMETRY_H
#define GEOMETRY_H

typedef struct {
	int size;
   int local_pieces; /* reticolo ristretto + output buffers */
   int *start; /* starting index for each piece. len=tot_pieces */
   int *end; /* ending index for each piece. len=tot_pieces */
   int *len; /* length for each local piece for unique points. len=local_pieces */
} spinor_descriptor;

#define LOCAL_SD_FOR(sd,i,j) \
	for(i=0;i<sd->local_pieces;i++) for(j=sd->start[i];j<sd->start[i]+sd->len[i];j++)

void geometry_init();
void geometry_blocked(void);
void geometry_blocked_noT(void);
void geometry_lexi(void);
void geometry_eo_lexi(void);

void test_geometry(void);
void print_geometry(void);
void print_wdmatrix(char *filename);

#endif
