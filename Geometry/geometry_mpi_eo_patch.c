/*******************************************************************************
*
* File geometry_mpi_eo_patch.c
*
* Calculate dimension of buffers
*
*******************************************************************************/

#ifdef GEOMETRY_MPI_EO_PATCH

#include <stdlib.h>
#include "geometry.h"
#include "global.h"
#include "memory.h"
#include "logger.h"



static int** coord=NULL;

static void init_buf_dim(geometry_descriptor *gd) {
  int n,i;
  
  gd->buf_dim=malloc(sizeof(int)*gd->nbuffers);
  
  for(n=0; n<gd->nbuffers; n++) {
    i=gd->rbuf_start[n];
    gd->buf_dim[n]=4;
    if(coord[i][0]<0 || coord[i][0]>=T) gd->buf_dim[n]--;
    if(coord[i][1]<0 || coord[i][1]>=X) gd->buf_dim[n]--;
    if(coord[i][2]<0 || coord[i][2]>=Y) gd->buf_dim[n]--;
    if(coord[i][3]<0 || coord[i][3]>=Z) gd->buf_dim[n]--;
  }
}


static void print_geometry_mpi_eo_patch(geometry_descriptor *gd) {
  int n, i;
  lprintf("GEOMETRY_PATCH",0,"#######################################\n");
  for(n=0; n<gd->nbuffers; n++) {
    lprintf("GEOMETRY_PATCH",0,"n=%d\t",n);
    i=gd->rbuf_start[n];
    lprintf("GEOMETRY_PATCH",0,"rbuf_start=%d(%d,%d,%d,%d)\t",i,coord[i][0],coord[i][1],coord[i][2],coord[i][3]);
    lprintf("GEOMETRY_PATCH",0,"rbuf_len=%d\t",gd->rbuf_len[n]);
    i=gd->sbuf_start[n];
    lprintf("GEOMETRY_PATCH",0,"sbuf_start=%d\t",i);
    lprintf("GEOMETRY_PATCH",0,"sbuf_len=%d\t",gd->sbuf_len[n]);
    lprintf("GEOMETRY_PATCH",0,"dim=%d\n",gd->buf_dim[n]);
  }
  lprintf("GEOMETRY_PATCH",0,"#######################################\n");
}


void init_geometry_mpi_eo_patch() {
  if(coord!=NULL) return;
  
  int c[4], i;

	coord = (int**)malloc(sizeof(int*)*glattice.gsize);
	coord[0] = (int*)malloc(sizeof(int)*glattice.gsize*4);
	for(i=0; i<glattice.gsize; i++)
		coord[i] = coord[0] + 4*i;

	/* Set coordinates */
	for(c[0]=-T_BORDER;c[0]<T+T_BORDER;c[0]++)
	for(c[1]=-X_BORDER;c[1]<X+X_BORDER;c[1]++)
	for(c[2]=-Y_BORDER;c[2]<Y+Y_BORDER;c[2]++)
	for(c[3]=-Z_BORDER;c[3]<Z+Z_BORDER;c[3]++) {
	  i=ipt(c[0],c[1],c[2],c[3]);
	  if(i==-1) continue;
	  coord[i][0]=c[0];
	  coord[i][1]=c[1];
	  coord[i][2]=c[2];
	  coord[i][3]=c[3];
	}
  
	init_buf_dim(&glattice);
	init_buf_dim(&glat_even);
	init_buf_dim(&glat_odd);
	
	print_geometry_mpi_eo_patch(&glattice);
	print_geometry_mpi_eo_patch(&glat_even);
	print_geometry_mpi_eo_patch(&glat_odd);	
}


void free_geometry_mpi_eo_patch() {
  if(coord==NULL) return;
  
  free(coord[0]);
  free(coord);
  
  free(glattice.buf_dim);
  free(glat_even.buf_dim);
  free(glat_odd.buf_dim);
}


#endif /*GEOMETRY_MPI_EO_PATCH*/

