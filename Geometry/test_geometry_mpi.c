/*******************************************************************************
*
* File test_geometry_mpi.c
*
* Test geometry_mpi
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "geometry.h"
#include "global.h"
#include "safe_mod.h"
#include "logger.h"
#include "error.h"

#define true 1
#define false 0

static int proc_up(int id, int dir)
{
  int ix,iy,iz,it;
  
  it = id%np_t;
  ix = (id/np_t)%np_x;
  iy = (id/(np_t*np_x))%np_y;
  iz = id/(np_t*np_x*np_y);
  
  if(dir == 0) it = (it +1)%np_t;
  if(dir == 1) ix = (ix +1)%np_x;
  if(dir == 2) iy = (iy +1)%np_y;
  if(dir == 3) iz = (iz +1)%np_z;
  
  return it + ix*np_t + iy*np_x*np_t + iz*np_y*np_x*np_t;
}

static int proc_down(int id, int dir)
{
  int ix,iy,iz,it;
  
  it = id%np_t;
  ix = (id/np_t)%np_x;
  iy = (id/(np_t*np_x))%np_y;
  iz = id/(np_t*np_x*np_y);
  
  if(dir == 0) it = (it -1+np_t)%np_t;
  if(dir == 1) ix = (ix -1+np_x)%np_x;
  if(dir == 2) iy = (iy -1+np_y)%np_y;
  if(dir == 3) iz = (iz -1+np_z)%np_z;
  
  return it + ix*np_t + iy*np_x*np_t + iz*np_y*np_x*np_t;
}

/* glattice is the geometry_descriptor instance fo the global lattice */
static int local_size[4];
static int buffer_thickness[4];
static int periodic_q[4];

static int** coord;

#define NOT_ASSIGNED 0
#define ORIGINAL     1
#define DUPLICATE    2

#define INNER    0
#define LBORDER  1
#define RBORDER  2
#define LBUFFER  3
#define RBUFFER  4

typedef struct {
	unsigned int c_type;    /* NOT_ASSIGNED ; ORIGINAL; DUPLICATE */
	unsigned int b_type[4]; /* INNER ; LBORDER ; RBORDER ; LBUFFER ; RBUFFER */
	unsigned int n_inner;   /* 0 ... 4 */
	unsigned int n_border;  /* 0 ... 4 */
	unsigned int n_buffer;  /* 0 ... 4 */
	unsigned int t_flag;
} site_info;
static site_info* descr;

static int loglevel = 0;

static void initialize_test() {
	int i;
	
	local_size[0] = T; local_size[1] = X; local_size[2] = Y; local_size[3] = Z;
	buffer_thickness[0] = T_BORDER; buffer_thickness[1] = X_BORDER; buffer_thickness[2] = Y_BORDER; buffer_thickness[3] = Z_BORDER;
	periodic_q[0] = (np_t==1) ? true : false; periodic_q[1] = (np_x==1) ? true : false; periodic_q[2] = (np_y==1) ? true : false; periodic_q[3] = (np_z==1) ? true : false; 
	
	coord = (int**)malloc(sizeof(int*)*glattice.gsize);
	coord[0] = (int*)malloc(sizeof(int)*glattice.gsize*4);
	for(i=0; i<glattice.gsize; i++)
		coord[i] = coord[0] + 4*i;
	
	descr = (site_info*)malloc(sizeof(site_info)*glattice.gsize);
}

static void finalize_test() {
	free(coord[0]);
	free(coord);
	free(descr);
}

static int in_glattice_q(int c[4]) {
	int howmanyoutside = 0;
	int i;
	int flag = 0;
	for(i=0; i<4; i++) {
		if(c[i] < 0 && c[i] >= -buffer_thickness[i]) howmanyoutside++;
		else if(c[i] >= local_size[i] && c[i] < local_size[i]+buffer_thickness[i]) howmanyoutside++;
		else if(c[i] < -buffer_thickness[i] || c[i] >= local_size[i]+buffer_thickness[i]){
			flag = 1;
			break;
		}
	}
	if(howmanyoutside > 2 && flag ==0) flag = 2;
	
	if(flag == 0) return true;
	return false;
}


static int set_nb_coordinates(int x) {
   int i, nb;
   int c[4];
   int counter = 0;

   if(x<0 || x >= glattice.gsize) {
      lprintf("TEST_GEOMETRY.C",0,"set_nb_coordinates(%d)\n",x);
      error(1,1,"set_nb_coordinates","I should not be here.");
   }

   if(descr[x].c_type == NOT_ASSIGNED) return 0;
   
   for(i=0; i<4; i++) {
      nb = iup(x,i);
      memcpy(c,coord[x],sizeof(int)*4);
      c[i]++;
      if(periodic_q[i]) c[i] = safe_mod(c[i],local_size[i]);
      if(!in_glattice_q(c)) continue;
      if(nb >= 0 && nb < glattice.gsize) {
         if(descr[nb].c_type == NOT_ASSIGNED) {
            descr[nb].c_type = ORIGINAL;
            memcpy(coord[nb],c,sizeof(int)*4);
            counter++;
         }
      } else {
         lprintf("TEST_GEOMETRY.C",0,"Bad candidate index %d for coordinates (%d,%d,%d,%d), reached from (%d,%d,%d,%d) with iup[%d].\n",
                 nb,c[0],c[1],c[2],c[3],coord[x][0],coord[x][1],coord[x][2],coord[x][3],i);
         return -1;
      }
   }

   for(i=0; i<4; i++) {
      nb = idn(x,i);
      memcpy(c,coord[x],sizeof(int)*4);
      c[i]--;
      if(periodic_q[i]) c[i] = safe_mod(c[i],local_size[i]);
      if(!in_glattice_q(c)) continue;
      if(nb >= 0 && nb < glattice.gsize) {
         if(descr[nb].c_type == NOT_ASSIGNED) {
            descr[nb].c_type = ORIGINAL;
            memcpy(coord[nb],c,sizeof(int)*4);
            counter++;
         }
      } else {
         lprintf("TEST_GEOMETRY.C",0,"Bad candidate index %d for coordinates (%d,%d,%d,%d), reached from (%d,%d,%d,%d) with idn[%d].\n",
                 nb,c[0],c[1],c[2],c[3],coord[x][0],coord[x][1],coord[x][2],coord[x][3],i);
         return -1;
      }
   }
   
   return counter;
}



int even_q(int c[4]) {
   return ((c[0]+c[1]+c[2]+c[3]
           +buffer_thickness[0]+buffer_thickness[1]+buffer_thickness[2]+buffer_thickness[3]
           +myid_sign)&1);
}

int odd_q(int c[4]) {
   return !((c[0]+c[1]+c[2]+c[3]
           +buffer_thickness[0]+buffer_thickness[1]+buffer_thickness[2]+buffer_thickness[3]
           +myid_sign)&1);
}



void test_geometry_descriptor(geometry_descriptor *gd, int(*in_subset_q)(int*)) {
	int i, j;
	int x, y;
	int *cx, *cy;
	int test_q;
	int origin[4] = {0,0,0,0};

	lprintf("TEST_GEOMETRY.C",loglevel,"gsize = %d\n", gd->gsize);

  if((*in_subset_q)(origin))
    lprintf("TEST_GEOMETRY.C",loglevel,"(*in_subset_q)(origin) = true\n");
  else
    lprintf("TEST_GEOMETRY.C",loglevel,"(*in_subset_q)(origin) = false\n");


	/* TEST: gd->gsize <= glattice.gsize */

	error(gd->gsize>glattice.gsize,1,"test_geometry.c","gd->gsize <= glattice.gsize... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"gd->gsize <= glattice.gsize... OK\n");


	lprintf("TEST_GEOMETRY.C",loglevel,"local_master_pieces = %d\n", gd->local_master_pieces);
	lprintf("TEST_GEOMETRY.C",loglevel,"total_master_pieces = %d\n", gd->total_master_pieces);
	lprintf("TEST_GEOMETRY.C",loglevel,"master_start = %p ... %p\n", gd->master_start, gd->master_start+gd->local_master_pieces-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"master_end = %p ... %p\n", gd->master_end, gd->master_end+gd->local_master_pieces-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"master_start[0] = %d\n", gd->master_start[0]);
	lprintf("TEST_GEOMETRY.C",loglevel,"master_end[0] = %d\n", gd->master_end[0]);
	


	/* TEST: Members master_start[0],master_end[0] describes all and only the inner sites */

	int inner_pieces = 1;

	for(x=0; x<glattice.gsize; x++)
		descr[x].t_flag = 0;
	
	test_q = true;
	for(i=0; i<inner_pieces; i++) {
		for(x=gd->master_start[i]; x<= gd->master_end[i] ; x++) {
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d in [%d,%d] (i=%d).\n",x,gd->master_start[i],gd->master_end[i],i);
				test_q = false;
			} else {
				if(descr[x].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included twice in the inner pieces.\n",x);
					test_q = false;
				}
				descr[x].t_flag = 1;
			}
		}
	}
	
	for(i=0; i<inner_pieces; i++)
		for(x=glattice.master_start[i]; x<= glattice.master_end[i] ; x++)
		  descr[x].t_flag += 2;
	
	for(x=0; x<glattice.gsize; x++) {
		cx = coord[x];
		if(descr[x].t_flag == 1) {
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included in gd inner pieces, but not in glattice inner pieces.\n",x);
			test_q = false;
		} else if(descr[x].t_flag == 2 && (*in_subset_q)(cx)) {
		  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not included in the inner pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		} else if(descr[x].t_flag == 3 && !(*in_subset_q)(cx)) {
		  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is included in the inner pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		}
	}
	
	error(!test_q,1,"test_geometry.c","Members master_start[0],master_end[0] describes all and only the inner sites... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members master_start[0],master_end[0] describes all and only the inner sites... OK\n");


	/* TEST: Members local_master_pieces,master_start,master_end describes all and only the local sites */

	for(x=0; x<glattice.gsize; x++)
		descr[x].t_flag = 0;
	
	test_q = true;
	for(i=0; i<gd->local_master_pieces; i++) {
/* 	  if(gd->master_end[i]<3000) */
/* 			lprintf("TEST_GEOMETRY.C",0,"[%d,%d] (i=%d) [%d,%d].\n",gd->master_start[i],gd->master_end[i],i,glattice.master_start[i],glattice.master_end[i]); */

		for(x=gd->master_start[i]; x<= gd->master_end[i] ; x++) {
		
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d in [%d,%d] (i=%d).\n",x,gd->master_start[i],gd->master_end[i],i);
				test_q = false;
			} else {
				if(descr[x].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included twice in the local pieces.\n",x);
					test_q = false;
				}
				descr[x].t_flag = 1;
			}
		}
	}
	
	for(i=0; i<glattice.local_master_pieces; i++)
		for(x=glattice.master_start[i]; x<= glattice.master_end[i] ; x++)
		  descr[x].t_flag += 2;
	
	for(x=0; x<glattice.gsize; x++) {
		cx = coord[x];
		if(descr[x].t_flag == 1) {
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included in gd local pieces, but not in glattice local pieces.\n",x);
			test_q = false;
		} else if(descr[x].t_flag == 2 && (*in_subset_q)(cx)) {
		  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not included in the local pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		} else if(descr[x].t_flag == 3 && !(*in_subset_q)(cx)) {
		  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is included in the local pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		}
	}
	
	error(!test_q,1,"test_geometry.c","Members local_master_pieces,master_start,master_end describes all and only the local sites... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members local_master_pieces,master_start,master_end describes all and only the local sites... OK\n");


	/* TEST: Members total_master_pieces,master_start,master_end describes all the sites */

	for(x=0; x<glattice.gsize; x++)
		descr[x].t_flag = 0;
	
	test_q = true;
	for(i=0; i<gd->total_master_pieces; i++) {
		for(x=gd->master_start[i]; x<= gd->master_end[i] ; x++) {
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d in [%d,%d] (i=%d).\n",x,gd->master_start[i],gd->master_end[i],i);
				test_q = false;
			} else {
				if(descr[x].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included twice in the total pieces.\n",x);
					test_q = false;
				}
				descr[x].t_flag = 1;
			}
		}
	}
	
	for(i=0; i<glattice.total_master_pieces; i++)
		for(x=glattice.master_start[i]; x<= glattice.master_end[i] ; x++)
		  descr[x].t_flag += 2;
	
	for(x=0; x<glattice.gsize; x++) {
		cx = coord[x];
		if(descr[x].t_flag == 1) {
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included in gd total pieces, but not in glattice total pieces.\n",x);
			test_q = false;
		} else if(descr[x].t_flag == 2 && (*in_subset_q)(cx)) {
		  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not included in the total pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		} else if(descr[x].t_flag == 3 && !(*in_subset_q)(cx)) {
		  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is included in the total pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		}
	}
	
	error(!test_q,1,"test_geometry.c","Members total_master_pieces,master_start,master_end describes all the sites... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members total_master_pieces,master_start,master_end describes all the sites... OK\n");


  /* TEST: Members ncopies,copy_from,copy_to,copy_len describes all the copy indexes */
	lprintf("TEST_GEOMETRY.C",loglevel,"ncopies = %d\n", gd->ncopies);
	lprintf("TEST_GEOMETRY.C",loglevel,"copy_len = %p ... %p\n", gd->copy_len, gd->copy_len+gd->ncopies-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"copy_from = %p ... %p\n", gd->copy_from, gd->copy_from+gd->ncopies-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"copy_to = %p ... %p\n", gd->copy_to, gd->copy_to+gd->ncopies-1);
  
	for(x=0; x<glattice.gsize; x++)
		descr[x].t_flag = 0;
	
	test_q = true;
	for(i=0; i<gd->ncopies; i++) {
		for(j=0; j<gd->copy_len[i]; j++) {
		  x=gd->copy_from[i]+j;
		  y=gd->copy_to[i]+j;
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d (i=%d, j=%d).\n",x,i,j);
				test_q = false;
			} else if(y<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"y = %d (i=%d, j=%d).\n",y,i,j);
				test_q = false;
			} else {
				if(descr[y].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is written twice in the copy process.\n",y);
					test_q = false;
				}
				descr[y].t_flag = 1;
				cx=coord[x];
				cy=coord[y];
				if(cx[0]!=cy[0] || cx[1]!=cy[1] || cx[2]!=cy[2] || cx[3]!=cy[3]) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d (%d,%d,%d,%d) is copied in index %d (%d,%d,%d,%d).\n",
					        x,cx[0],cx[1],cx[2],cx[3],
					        y,cy[0],cy[1],cy[2],cy[3]);
					test_q = false;
				}
			}
		}
	}
	
	for(i=0; i<glattice.ncopies; i++)
		for(j=0; j<glattice.copy_len[i]; j++) {
		  y=glattice.copy_to[i]+j;
		  descr[y].t_flag += 2;
		}
	
	for(x=0; x<glattice.gsize; x++) {
		cx = coord[x];
		if(descr[x].t_flag == 1) {
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is a gd copy index, but not a glattice copy index.\n",x);
			test_q = false;
		} else if(descr[x].t_flag == 2 && (*in_subset_q)(cx)) {
		  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not a gd copy index.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		} else if(descr[x].t_flag == 3 && !(*in_subset_q)(cx)) {
		  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is a gd copy index.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		}
	}
	
	error(!test_q,1,"test_geometry.c","Members ncopies,copy_from,copy_to,copy_len describes all the copy indexes... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members ncopies,copy_from,copy_to,copy_len describes all the copy indexes... OK\n");

 	/* TEST: Members nbuffers,rbuf*,sbuf* describes the sending/receiving buffers */
	lprintf("TEST_GEOMETRY.C",loglevel,"nbuffers = %d\n", gd->nbuffers);
	lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_len = %p ... %p\n", gd->rbuf_len, gd->rbuf_len+gd->nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_len = %p ... %p\n", gd->sbuf_len, gd->sbuf_len+gd->nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_from_proc = %p ... %p\n", gd->rbuf_from_proc, gd->rbuf_from_proc+gd->nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_start = %p ... %p\n", gd->rbuf_start, gd->rbuf_start+gd->nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_to_proc = %p ... %p\n", gd->sbuf_to_proc, gd->sbuf_to_proc+gd->nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_start = %p ... %p\n", gd->sbuf_start, gd->sbuf_start+gd->nbuffers-1);
  
 	test_q = true;
  for(i=0; i<gd->nbuffers; i++) {
    if(gd->rbuf_from_proc[i] != glattice.rbuf_from_proc[i]) {
				lprintf("TEST_GEOMETRY.C",0,"gd->rbuf_from_proc[%d] = %d  but glattice.rbuf_from_proc[%d] = %d\n",
				        i,gd->rbuf_from_proc[i],i,glattice.rbuf_from_proc[i]);
				test_q = false;
    } else if(gd->sbuf_to_proc[i] != glattice.sbuf_to_proc[i]) {
				lprintf("TEST_GEOMETRY.C",0,"gd->sbuf_to_proc[%d] = %d  but glattice.sbuf_to_proc[%d] = %d\n",
				        i,gd->sbuf_to_proc[i],i,glattice.sbuf_to_proc[i]);
				test_q = false;
    }
  }

  for(i=0; i<gd->nbuffers; i++) {
  
   	for(x=0; x<glattice.gsize; x++)
    	descr[x].t_flag = 0;

		for(j=0; j<gd->sbuf_len[i]; j++) {
			x=gd->sbuf_start[i]+j;
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d sbuf*(i=%d, j=%d).\n",x,i,j);
				test_q = false;
			} else {
				if(descr[x].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is twice in the gd->sbuf*[%d].\n",x,i);
					test_q = false;
				}
				descr[x].t_flag = 1;
			}
		}

		for(j=0; j<glattice.sbuf_len[i]; j++) {
		  x=glattice.sbuf_start[i]+j;
		  descr[x].t_flag += 2;
		}

		for(j=0; j<glattice.sbuf_len[i]; j++) {
		  x=glattice.sbuf_start[i]+j;
		  cx = coord[x];
		  if(descr[x].t_flag == 1) {
			  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is in gd->sbuf*[%d], but not in glattice.sbuf*[%i].\n",x,i,i);
			  test_q = false;
		  } else if(descr[x].t_flag == 2 && (*in_subset_q)(cx)) {
		    lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not in gd->sbuf*[%d].\n",x,cx[0],cx[1],cx[2],cx[3],i);
			  test_q = false;
		  } else if(descr[x].t_flag == 3 && !(*in_subset_q)(cx)) {
		    lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is in gd->sbuf*[%d].\n",x,cx[0],cx[1],cx[2],cx[3],i);
			  test_q = false;
		  }
	  }

	}
	
  for(i=0; i<gd->nbuffers; i++) {
  
   	for(x=0; x<glattice.gsize; x++)
    	descr[x].t_flag = 0;

		for(j=0; j<gd->rbuf_len[i]; j++) {
			x=gd->rbuf_start[i]+j;
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d rbuf*(i=%d, j=%d).\n",x,i,j);
				test_q = false;
			} else {
				if(descr[x].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is twice in the gd->rbuf*[%d].\n",x,i);
					test_q = false;
				}
				descr[x].t_flag = 1;
			}
		}

		for(j=0; j<glattice.rbuf_len[i]; j++) {
		  x=glattice.rbuf_start[i]+j;
		  descr[x].t_flag += 2;
		}

		for(j=0; j<glattice.rbuf_len[i]; j++) {
		  x=glattice.rbuf_start[i]+j;
		  cx = coord[x];
		  if(descr[x].t_flag == 1) {
			  lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is in gd->rbuf*[%d], but not in glattice.rbuf*[%i].\n",x,i,i);
			  test_q = false;
		  } else if(descr[x].t_flag == 2 && (*in_subset_q)(cx)) {
		    lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not in gd->rbuf*[%d].\n",x,cx[0],cx[1],cx[2],cx[3],i);
			  test_q = false;
		  } else if(descr[x].t_flag == 3 && !(*in_subset_q)(cx)) {
		    lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is in gd->rbuf*[%d].\n",x,cx[0],cx[1],cx[2],cx[3],i);
			  test_q = false;
		  }
	  }

	}
	
	error(!test_q,1,"test_geometry.c","Members nbuffers,rbuf*,sbuf* describes the sending/receiving buffers... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members nbuffers,rbuf*,sbuf* describes the sending/receiving buffers... OK\n");

}


void test_glattice() {
	int i, j, k;
	int x, y;
	int *cx, *cy;
	int test_q;
	
	lprintf("TEST_GEOMETRY.C",loglevel,"Memory allocation... OK\n");
	lprintf("TEST_GEOMETRY.C",loglevel,"gsize = %d\n", glattice.gsize);
	lprintf("TEST_GEOMETRY.C",loglevel,"myid_sign = %d\n", myid_sign);
	lprintf("TEST_GEOMETRY.C",loglevel,"local_size = {%d,%d,%d,%d}\n", local_size[0], local_size[1], local_size[2], local_size[3]);
	lprintf("TEST_GEOMETRY.C",loglevel,"buffer_thickness = {%d,%d,%d,%d}\n", buffer_thickness[0], buffer_thickness[1], buffer_thickness[2], buffer_thickness[3]);
	lprintf("TEST_GEOMETRY.C",loglevel,"periodic_q = {%d,%d,%d,%d}\n", periodic_q[0], periodic_q[1], periodic_q[2], periodic_q[3]);
	lprintf("TEST_GEOMETRY.C",loglevel,"coord = %p ... %p\n", coord, coord+glattice.gsize-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"coord[0] = %p ... %p\n", coord[0], coord[0]+4*glattice.gsize-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"descr = %p ... %p\n", descr, descr+glattice.gsize-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"iup = %p ... %p\n", iup, iup+4*glattice.gsize-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"idn = %p ... %p\n", idn, idn+4*glattice.gsize-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"ipt = %p ... %p\n", ipt, ipt+VOLUME-1);
	
	memset(descr,'\0',sizeof(site_info)*glattice.gsize);
	
	/* Set the origin */
	int origin = ipt(0,0,0,0);
	coord[origin][0] = coord[origin][1] = coord[origin][2] = coord[origin][3] = 0;
	descr[origin].c_type = ORIGINAL;
	lprintf("TEST_GEOMETRY.C",loglevel,"Origin found at index %d.\n",origin);
	
	/* Set the coordinates for the points reached by iup, idn */
   int howmany;
   test_q=true;
   do {
      howmany=0;
      for(x=0; x<glattice.gsize; x++) {
         int ret = set_nb_coordinates(x);
         if(ret==-1) test_q = false;
         else howmany += ret;
      }
   } while(howmany != 0 && test_q == true);
	
	error(!test_q,1,"test_geometry.c","Set the coordinates for the points reached by iup, idn... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Set the coordinates for the points reached by iup, idn... OK\n");
	
	/* TEST: No duplicate coordinates for the original indexes */
	test_q = true;
	for(x=0; x<glattice.gsize; x++) {
		if(descr[x].c_type == NOT_ASSIGNED) continue;
		cx = coord[x];
		for(y=x+1; y<glattice.gsize; y++) {
			if(descr[y].c_type == NOT_ASSIGNED) continue;
			cy = coord[y];
			if(cx[0]==cy[0] && cx[1]==cy[1] && cx[2]==cy[2] && cx[3]==cy[3]) {
				lprintf("TEST_GEOMETRY.C",0,"Duplicate coordinates for the indexes %d and %d\n",x,y);
				lprintf("TEST_GEOMETRY.C",0,"Duplicate coordinates: ( %d , %d , %d , %d )\n",cx[0],cx[1],cx[2],cx[3]);
				test_q = false;
			}
		}
	}
	error(!test_q,1,"test_geometry.c","No duplicate coordinates for the original indexes... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"No duplicate coordinates for the original indexes... OK\n");
	
	/* Set the coordinates for the duplicate points */
	/* TEST: Each duplicate point has just one original point */
	lprintf("TEST_GEOMETRY.C",loglevel,"ncopies = %d\n", glattice.ncopies);
	lprintf("TEST_GEOMETRY.C",loglevel,"copy_len = %p ... %p\n", glattice.copy_len, glattice.copy_len+glattice.ncopies-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"copy_from = %p ... %p\n", glattice.copy_from, glattice.copy_from+glattice.ncopies-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"copy_to = %p ... %p\n", glattice.copy_to, glattice.copy_to+glattice.ncopies-1);

	test_q = true;
	for(i=0; i<glattice.ncopies; i++) {
		for(j=0; j<glattice.copy_len[i]; j++) {
			x = glattice.copy_from[i]+j;
			y = glattice.copy_to[i]+j;
			if(x<0 || x>=glattice.gsize || y<0 || y>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"I'm trying to copy from %d to %d.\n",x,y);
				test_q = false;
			} else {
				cx = coord[x];
				cy = coord[y];
				if(descr[y].c_type != NOT_ASSIGNED) {
					lprintf("TEST_GEOMETRY.C",0,"Coordinates for the index %d already assigned (i=%d,j=%d)\n",y,i,j);
					lprintf("TEST_GEOMETRY.C",0,"Coordinates of %d: ( %d , %d , %d , %d )\n",x,cx[0],cx[1],cx[2],cx[3]);
					lprintf("TEST_GEOMETRY.C",0,"Coordinates of %d: ( %d , %d , %d , %d )\n",y,cy[0],cy[1],cy[2],cy[3]);
					test_q = false;
				}
				descr[y].c_type = DUPLICATE;
				memcpy(cy,cx,4*sizeof(int));
			}

		}
	}
	error(!test_q,1,"test_geometry.c","Each duplicate point has just one original point... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Each duplicate point has just one original point... OK\n");

	/* TEST: Coordinates are assigned to each point */
	test_q = true;
	for(x=0; x<glattice.gsize; x++) {
		if(descr[x].c_type == NOT_ASSIGNED) {
			lprintf("TEST_GEOMETRY.C",0,"Coordinates are not assigned to the index %d\n",x);
			test_q = false;
		}
	}
	error(!test_q,1,"test_geometry.c","Coordinates are assigned to each point... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Coordinates are assigned to each point... OK\n");
	
	/* Select borders and buffers */
	/* TEST: Assigned coordinates are in the right ranges */
	test_q = true;
	for(x=0; x<glattice.gsize; x++) {
		cx = coord[x];
		descr[x].n_inner = 0;
		descr[x].n_border = 0;
		descr[x].n_buffer = 0;
		for(i=0; i<4; i++) {
			if(cx[i] >= buffer_thickness[i] && cx[i] < local_size[i]-buffer_thickness[i]){
				descr[x].b_type[i] = INNER;
				descr[x].n_inner++;
			} else if(cx[i] >= 0 && cx[i] < buffer_thickness[i]) {
				descr[x].b_type[i] = LBORDER;
				descr[x].n_border++;
			} else if(cx[i] >= local_size[i]-buffer_thickness[i] && cx[i] < local_size[i]) {
				descr[x].b_type[i] = RBORDER;
				descr[x].n_border++;
			} else if(cx[i] >= -buffer_thickness[i] && cx[i] < 0) {
				descr[x].b_type[i] = LBUFFER;
				descr[x].n_buffer++;
			} else if(cx[i] >= local_size[i] && cx[i] < local_size[i]+buffer_thickness[i]) {
				descr[x].b_type[i] = RBUFFER;
				descr[x].n_buffer++;
			} else {
				lprintf("TEST_GEOMETRY.C",0,"Coordinates of %d: ( %d , %d , %d , %d )\n",x,cx[0],cx[1],cx[2],cx[3]);
				test_q = false;
			}
		}
	}
	error(!test_q,1,"test_geometry.c","Assigned coordinates are in the right ranges... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Assigned coordinates are in the right ranges... OK\n");
	
	/* TEST: ipt gives the right map between coordinates and local indexes */
	test_q = true;
	
	for(x=0; x<glattice.gsize; x++) {
		if(descr[x].n_buffer != 0 || descr[x].c_type != ORIGINAL) continue;
		cx = coord[x];
		if(ipt(cx[0],cx[1],cx[2],cx[3]) != x) {
			lprintf("TEST_GEOMETRY.C",loglevel,"ipt(%d,%d,%d,%d)=%d ; but x=%d .\n",cx[0],cx[1],cx[2],cx[3],ipt(cx[0],cx[1],cx[2],cx[3]),x);
			test_q = false;
		}
	}
	
	int c0, c1, c2, c3;
	for(c0=0; c0<local_size[0]; c0++)
	for(c1=0; c1<local_size[1]; c1++)
	for(c2=0; c2<local_size[2]; c2++)
	for(c3=0; c3<local_size[3]; c3++) {
		x = ipt(c0,c1,c2,c3);
		if(x<0) {
			lprintf("TEST_GEOMETRY.C",loglevel,"ipt(%d,%d,%d,%d)=%d < 0 .\n",c0,c1,c2,c3,x);
			test_q = false;
		} else if(x>=glattice.gsize) {
			lprintf("TEST_GEOMETRY.C",loglevel,"ipt(%d,%d,%d,%d)=%d >= gsize .\n",c0,c1,c2,c3,x);
			test_q = false;
		} else if(descr[x].n_buffer != 0) {
			lprintf("TEST_GEOMETRY.C",loglevel,"ipt(%d,%d,%d,%d)=%d which is non-local.\n",c0,c1,c2,c3,x);
			test_q = false;
		} else if(descr[x].c_type != ORIGINAL) {
			lprintf("TEST_GEOMETRY.C",loglevel,"ipt(%d,%d,%d,%d)=%d which is a duplicate.\n",c0,c1,c2,c3,x);
			test_q = false;
		} else {
			cx = coord[x];
			if(c0 != cx[0] || c1 != cx[1] || c2 != cx[2] || c3 != cx[3]) {
				lprintf("TEST_GEOMETRY.C",loglevel,"ipt(%d,%d,%d,%d)=%d but coord[%x]=(%d,%d,%d,%d).\n",c0,c1,c2,c3,x,x,cx[0],cx[1],cx[2],cx[3]);
				test_q = false;
			}
		}
	}
	
	error(!test_q,1,"test_geometry.c","ipt gives the right map between coordinates and local indexes... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"ipt gives the right map between coordinates and local indexes... OK\n");
	
	/* TEST: iup, idn give the right maps between neighbours */
	int id[4][4];
	for(i=0; i<4; i++) {
		for(j=0; j<4; j++) id[i][j] = 0;
		id[i][i] = 1;
	}

	test_q = true;
	for(x=0; x<glattice.gsize; x++) {
		if(descr[x].c_type != ORIGINAL) continue;
		cx = coord[x];
		for(i=0; i<4; i++) {
			y = iup(x,i);
			int c[4];
			memcpy(c,coord[x],sizeof(int)*4);
			c[i]++;
			if(periodic_q[i]) c[i] = safe_mod(c[i],local_size[i]);
			if(!in_glattice_q(c)) continue;
			cy = coord[y];
			if(cy[0] != cx[0]+id[i][0] &&
			   cy[1] != cx[1]+id[i][1] &&
			   cy[2] != cx[2]+id[i][2] &&
			   cy[3] != cx[3]+id[i][3]) {
					lprintf("TEST_GEOMETRY.C",loglevel,"iup[%d](%d,%d,%d,%d)=(%d,%d,%d,%d).\n",i,cx[0],cx[1],cx[2],cx[3],cy[0],cy[1],cy[2],cy[3]);
					test_q = false;
			}
		}
		for(i=0; i<4; i++) {
			y = idn(x,i);
			int c[4];
			memcpy(c,coord[x],sizeof(int)*4);
			c[i]--;
			if(periodic_q[i]) c[i] = safe_mod(c[i],local_size[i]);
			if(!in_glattice_q(c)) continue;
			cy = coord[y];
			if(cy[0] != cx[0]-id[i][0] &&
			   cy[1] != cx[1]-id[i][1] &&
			   cy[2] != cx[2]-id[i][2] &&
			   cy[3] != cx[3]-id[i][3]) {
					lprintf("TEST_GEOMETRY.C",loglevel,"idn[%d](%d,%d,%d,%d)=(%d,%d,%d,%d).\n",i,cx[0],cx[1],cx[2],cx[3],cy[0],cy[1],cy[2],cy[3]);
					test_q = false;
			}
		}
	}
	error(!test_q,1,"test_geometry.c","iup, idn give the right maps between neighbours... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"iup, idn give the right maps between neighbours... OK\n");

	lprintf("TEST_GEOMETRY.C",loglevel,"local_master_pieces = %d\n", glattice.local_master_pieces);
	lprintf("TEST_GEOMETRY.C",loglevel,"total_master_pieces = %d\n", glattice.total_master_pieces);
	lprintf("TEST_GEOMETRY.C",loglevel,"master_start = %p ... %p\n", glattice.master_start, glattice.master_start+glattice.local_master_pieces-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"master_end = %p ... %p\n", glattice.master_end, glattice.master_end+glattice.local_master_pieces-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"master_start[0] = %d\n", glattice.master_start[0]);
	lprintf("TEST_GEOMETRY.C",loglevel,"master_end[0] = %d\n", glattice.master_end[0]);

	/* TEST: Members master_start[0],master_end[0] describes all and only the inner sites */

	int inner_pieces = 1;
	
	test_q = true;
	int inner_volume = (local_size[0]-2*buffer_thickness[0])*(local_size[1]-2*buffer_thickness[1])*
	                   (local_size[2]-2*buffer_thickness[2])*(local_size[3]-2*buffer_thickness[3]);
	if(glattice.master_end[0]-glattice.master_start[0]+1 != inner_volume) {
	   test_q = false;
	   lprintf("TEST_GEOMETRY.C",loglevel,"master_end[0]-master_start[0]+1 = %d but it should be %d\n",
	           glattice.master_end[0]-glattice.master_start[0]+1,inner_volume);
	}
	error(!test_q,1,"test_geometry.c","Members master_start[0],master_end[0] describes all and only the inner sites... FAILED");

	for(x=0; x<glattice.gsize; x++)
		descr[x].t_flag = 0;
	
	test_q = true;
	for(i=0; i<inner_pieces; i++) {
		for(x=glattice.master_start[i]; x<= glattice.master_end[i] ; x++) {
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d in [%d,%d] (i=%d).\n",x,glattice.master_start[i],glattice.master_end[i],i);
				test_q = false;
			} else {
				if(descr[x].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included twice in the inner pieces.\n",x);
					test_q = false;
				}
				descr[x].t_flag = 1;
			}
		}
	}
	
	for(x=0; x<glattice.gsize; x++) {
		if(descr[x].n_inner == 4 && descr[x].t_flag == 0 && descr[x].c_type == ORIGINAL) {
			cx = coord[x];
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not included in the inner pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		} else if((descr[x].n_inner != 4 || descr[x].c_type != ORIGINAL) && descr[x].t_flag == 1) {
			cx = coord[x];
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is included in the inner pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		}
	}
	error(!test_q,1,"test_geometry.c","Members master_start[0],master_end[0] describes all and only the inner sites... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members master_start[0],master_end[0] describes all and only the inner sites... OK\n");


	/* TEST: Members local_master_pieces,master_start,master_end describes all and only the local sites */
	
	test_q = true;
	int local_volume = local_size[0]*local_size[1]*local_size[2]*local_size[3];
	int sum_local_len = 0;
	for(i=0; i<glattice.local_master_pieces; i++)
		sum_local_len += glattice.master_end[i]-glattice.master_start[i]+1;
	
	if(sum_local_len != local_volume) {
	   test_q = false;
	   lprintf("TEST_GEOMETRY.C",loglevel,"The sum of the lengths of all the local pieces is %d but should be %d\n",sum_local_len,local_volume);
	}
	error(!test_q,1,"test_geometry.c","Members local_master_pieces,master_start,master_end describes all and only the local site... FAILED");

	for(x=0; x<glattice.gsize; x++)
		descr[x].t_flag = 0;
	
	test_q = true;
	for(i=0; i<glattice.local_master_pieces; i++) {
		for(x=glattice.master_start[i]; x<=glattice.master_end[i]; x++) {
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d in [%d,%d] (i=%d).\n",x,glattice.master_start[i],glattice.master_end[i],i);
				test_q = false;
			} else {
				if(descr[x].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included twice in the local pieces.\n",x);
					test_q = false;
				}
				descr[x].t_flag = 1;
			}
		}
	}
	
	for(x=0; x<glattice.gsize; x++) {
		if(descr[x].n_buffer == 0 && descr[x].t_flag == 0 && descr[x].c_type == ORIGINAL) {
			cx = coord[x];
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not included in the local pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		} else if((descr[x].n_buffer != 0 || descr[x].c_type != ORIGINAL) && descr[x].t_flag == 1) {
			cx = coord[x];
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is included in the local pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		}
	}
	error(!test_q,1,"test_geometry.c","Members local_master_pieces,master_start,master_end describes all and only the local site... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members local_master_pieces,master_start,master_end describes all and only the local site... OK\n");


	/* TEST: Members total_master_pieces,master_start,master_end describes all the sites */

	test_q = true;
	int total_volume = local_size[0]*local_size[1]*local_size[2]*local_size[3]+
	                   2*buffer_thickness[0]*local_size[1]*local_size[2]*local_size[3]+
	                   2*local_size[0]*buffer_thickness[1]*local_size[2]*local_size[3]+
	                   2*local_size[0]*local_size[1]*buffer_thickness[2]*local_size[3]+
	                   2*local_size[0]*local_size[1]*local_size[2]*buffer_thickness[3]+
	                   4*buffer_thickness[0]*buffer_thickness[1]*local_size[2]*local_size[3]+
	                   4*buffer_thickness[0]*local_size[1]*buffer_thickness[2]*local_size[3]+
	                   4*buffer_thickness[0]*local_size[1]*local_size[2]*buffer_thickness[3]+
	                   4*local_size[0]*buffer_thickness[1]*buffer_thickness[2]*local_size[3]+
	                   4*local_size[0]*buffer_thickness[1]*local_size[2]*buffer_thickness[3]+
	                   4*local_size[0]*local_size[1]*buffer_thickness[2]*buffer_thickness[3];
	int sum_total_len = 0;
	for(i=0; i<glattice.total_master_pieces; i++)
		sum_total_len += glattice.master_end[i]-glattice.master_start[i]+1;
	
	if(sum_total_len != total_volume) {
	   test_q = false;
	   lprintf("TEST_GEOMETRY.C",loglevel,"The sum of the lengths of all the total pieces is %d but should be %d\n",sum_total_len,total_volume);
	}
	error(!test_q,1,"test_geometry.c","Members total_master_pieces,master_start,master_end describes all the sites... FAILED");

	for(x=0; x<glattice.gsize; x++)
		descr[x].t_flag = 0;
	
	test_q = true;
	for(i=0; i<glattice.total_master_pieces; i++) {
		for(x=glattice.master_start[i]; x<=glattice.master_end[i]; x++) {
			if(x<0 || x>=glattice.gsize) {
				lprintf("TEST_GEOMETRY.C",0,"x = %d in [%d,%d] (i=%d).\n",x,glattice.master_start[i],glattice.master_end[i],i);
				test_q = false;
			} else {
				if(descr[x].t_flag == 1) {
					lprintf("TEST_GEOMETRY.C",loglevel,"Index %d is included twice in the global pieces.\n",x);
					test_q = false;
				}
				descr[x].t_flag = 1;
			}
		}
	}
	
	for(x=0; x<glattice.gsize; x++) {
		if(descr[x].t_flag == 0 && descr[x].c_type == ORIGINAL) {
			cx = coord[x];
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is not included in the global pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		} else if(descr[x].c_type != ORIGINAL && descr[x].t_flag == 1) {
			cx = coord[x];
			lprintf("TEST_GEOMETRY.C",loglevel,"Index %d of coordinates (%d,%d,%d,%d) is included in the global pieces.\n",x,cx[0],cx[1],cx[2],cx[3]);
			test_q = false;
		}
	}
	error(!test_q,1,"test_geometry.c","Members total_master_pieces,master_start,master_end describes all the sites... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members total_master_pieces,master_start,master_end describes all the sites... OK\n");

	/* TEST: Members nbuffers,rbuf*,sbuf* describes the sending/receiving buffers */
	lprintf("TEST_GEOMETRY.C",loglevel,"nbuffers = %d\n", glattice.nbuffers);
	lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_len = %p ... %p\n", glattice.rbuf_len, glattice.rbuf_len+glattice.nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_len = %p ... %p\n", glattice.sbuf_len, glattice.sbuf_len+glattice.nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_from_proc = %p ... %p\n", glattice.rbuf_from_proc, glattice.rbuf_from_proc+glattice.nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_start = %p ... %p\n", glattice.rbuf_start, glattice.rbuf_start+glattice.nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_to_proc = %p ... %p\n", glattice.sbuf_to_proc, glattice.sbuf_to_proc+glattice.nbuffers-1);
	lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_start = %p ... %p\n", glattice.sbuf_start, glattice.sbuf_start+glattice.nbuffers-1);


	test_q = true;
	for(i=0; i<glattice.nbuffers; i++) {

#define LOCAL 100
		int sbuf_mask[4], rbuf_mask[4];
		int dimension = 0;

		/* glattice.rbuf_len[i] must be equal to glattice.sbuf_len[i] */

    if(glattice.rbuf_len[i] != glattice.sbuf_len[i]) {
			lprintf("TEST_GEOMETRY.C",loglevel,"glattice.rbuf_len[%d]=%d ; glattice.sbuf_len[%d]=%d .\n",i,glattice.rbuf_len[i],i,glattice.sbuf_len[i]);
			test_q = false;
    }
		
		/* Determine which border is described by sbuf*[i] - result in sbuf_mask, dimension */
		
		x = glattice.sbuf_start[i];
		for(k=0;k<4;k++)
			sbuf_mask[k] = descr[x].b_type[k];

		for(j=0; j<glattice.sbuf_len[i]; j++) {
			x = glattice.sbuf_start[i]+j;
			for(k=0;k<4;k++) {
				if(descr[x].b_type[k] == LBUFFER || descr[x].b_type[k] == RBUFFER) {
					lprintf("TEST_GEOMETRY.C",loglevel,"glattice.sbuf_start[%d]+%d=%x but it is in a buffer.\n",i,j,x);
					test_q = false;
				} else if(descr[x].b_type[k] != sbuf_mask[k])
					sbuf_mask[k] = LOCAL;
			}
		}
		
		/* Accordingly with sid*[i], determine which buffer is expected to be described by rbuf*[i] - result in rbuf_mask */
		for(k=0;k<4;k++) {
			if(sbuf_mask[k] == LBORDER) rbuf_mask[k] = RBUFFER;
			else if(sbuf_mask[k] == RBORDER) rbuf_mask[k] = LBUFFER;
			else {
				rbuf_mask[k] = LOCAL;
				dimension++;
			}
		}
		
		/* TEST: dimension can be only 2 or 3 */
		if(dimension != 3 && dimension != 2) {
			lprintf("TEST_GEOMETRY.C",loglevel,"sbuf*[%d] has dimension=%d and sbuf_mask={%d,%d,%d,%d}.\n",i,dimension,
			        sbuf_mask[0],sbuf_mask[1],sbuf_mask[2],sbuf_mask[3]);
			test_q = false;
		}

		/* TEST: rbuf*[i] describes the right buffer, accordingly to rbuf_mask */
		for(j=0; j<glattice.rbuf_len[i]; j++) {
			x = glattice.rbuf_start[i]+j;
			cx = coord[x];
			for(k=0;k<4;k++) {
				if(rbuf_mask[k] == RBUFFER && descr[x].b_type[k] != RBUFFER) {
					lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_mask[%d]=RBUFFER but coord[%d=glattice.rbuf_start[%d]+%d]=(%d,%d,%d,%d).\n",
					        k,x,i,j,cx[0],cx[1],cx[2],cx[3]);
					test_q = false;
				} else if(rbuf_mask[k] == LBUFFER && descr[x].b_type[k] != LBUFFER) {
					lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_mask[%d]=LBUFFER but coord[%d=glattice.rbuf_start[%d]+%d]=(%d,%d,%d,%d).\n",
					        k,x,i,j,cx[0],cx[1],cx[2],cx[3]);
					test_q = false;
				} else if(rbuf_mask[k] == LOCAL && (descr[x].b_type[k] == LBUFFER || descr[x].b_type[k] == RBUFFER)) {
					lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_mask[%d]=LOCAL but coord[%d=glattice.rbuf_start[%d]+%d]=(%d,%d,%d,%d).\n",
					        k,x,i,j,cx[0],cx[1],cx[2],cx[3]);
					test_q = false;
				}
			}
		}

		/* TEST: All the points of the determined border are included in sbuf*[i], at least once */		
		int sid_crange[2][4];
		for(k=0;k<4;k++) {
			if(sbuf_mask[k] == LBORDER) {
				sid_crange[0][k] = 0;
				sid_crange[1][k] = buffer_thickness[k]-1;
			} else if(sbuf_mask[k] == RBORDER) {
				sid_crange[0][k] = local_size[k]-buffer_thickness[k];
				sid_crange[1][k] = local_size[k]-1;
			} else {
				sid_crange[0][k] = buffer_thickness[k];
				sid_crange[1][k] = local_size[k]-buffer_thickness[k]-1;
			}
		}
		
		for(c0=sid_crange[0][0]; c0<=sid_crange[1][0]; c0++)
		for(c1=sid_crange[0][1]; c1<=sid_crange[1][1]; c1++)
		for(c2=sid_crange[0][2]; c2<=sid_crange[1][2]; c2++)
		for(c3=sid_crange[0][3]; c3<=sid_crange[1][3]; c3++) {
			int isthere_q = false;
			for(j=0; j<glattice.sbuf_len[i] && isthere_q == false; j++) {
				x = glattice.sbuf_start[i]+j;
				cx = coord[x];
				if(c0 == cx[0] && c1 == cx[1] && c2 == cx[2] && c3 == cx[3])
					isthere_q = true;
			}
			if(!isthere_q) {
				lprintf("TEST_GEOMETRY.C",loglevel,"The point (%d,%d,%d,%d) has the right sbuf_mask={%d,%d,%d,%d} but it was not found in sid*[%d].\n",
				        c0,c1,c2,c3,sbuf_mask[0],sbuf_mask[1],sbuf_mask[2],sbuf_mask[3],i);
				test_q = false;
			}
		}
		
		/* TEST: rbuf*[i] is obtained from sbuf*[i] by a shift */
		int shift[4];
		for(k=0;k<4;k++) {
			if(sbuf_mask[k] == LBORDER) shift[k] = local_size[k];
			else if(sbuf_mask[k] == RBORDER) shift[k] = -local_size[k];
			else shift[k] = 0;
		}

		for(j=0; j<glattice.rbuf_len[i]; j++) {
			x = glattice.sbuf_start[i]+j;
			y = glattice.rbuf_start[i]+j;
			cx = coord[x];
			cy = coord[y];
			if(cy[0] != cx[0]+shift[0] ||
			   cy[1] != cx[1]+shift[1] ||
			   cy[2] != cx[2]+shift[2] ||
			   cy[3] != cx[3]+shift[3]) {
				lprintf("TEST_GEOMETRY.C",loglevel,"The point (i=%d,j=%d) (%d,%d,%d,%d) is sent to (%d,%d,%d,%d) != (%d+%d,%d+%d,%d+%d,%d+%d).\n",
				        i,j,cx[0],cx[1],cx[2],cx[3],cy[0],cy[1],cy[2],cy[3],cx[0],shift[0],cx[1],shift[1],cx[2],shift[2],cx[3],shift[3]);
				test_q = false;
			}
		}
		
		/* TEST: rid*[i] is received from the right processor */
		int id;
		
		if(glattice.rbuf_from_proc[i] == myid && !periodic_q[i]) {
			lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_from_proc[%d]=myid=%d (rbuf_mask={%d,%d,%d,%d}) but not periodic BC.\n",i,myid,rbuf_mask[0],rbuf_mask[1],rbuf_mask[2],rbuf_mask[3]);
			test_q = false;
		}
		if(glattice.rbuf_from_proc[i] < 0 || glattice.rbuf_from_proc[i] >= np_x*np_y*np_z*np_t) {
			lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_from_proc[%d]=%d (rbuf_mask={%d,%d,%d,%d}) out of range.\n", i,glattice.rbuf_from_proc[i],rbuf_mask[0],rbuf_mask[1],rbuf_mask[2],rbuf_mask[3]);
			test_q = false;
		}
		id = myid;
		for(k=0; k<4; k++) {
			if(shift[k] > 0) id = proc_up(id,k);
			else if(shift[k] < 0) id = proc_down(id,k);
		}
		if(id != glattice.rbuf_from_proc[i]) {
			lprintf("TEST_GEOMETRY.C",loglevel,"rbuf_from_proc[%d]=%d (rbuf_mask={%d,%d,%d,%d}) , expected %d.\n", i,glattice.rbuf_from_proc[i],rbuf_mask[0],rbuf_mask[1],rbuf_mask[2],rbuf_mask[3],id);
			test_q = false;
		}

		/* TEST: sid*[i] is sent to the right processor */
		if(glattice.sbuf_to_proc[i] == myid && !periodic_q[i]) {
			lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_to_proc[%d]=myid=%d (sbuf_mask={%d,%d,%d,%d}) but not periodic BC.\n",i,myid,sbuf_mask[0],sbuf_mask[1],sbuf_mask[2],sbuf_mask[3]);
			test_q = false;
		}
		if(glattice.sbuf_to_proc[i] < 0 || glattice.sbuf_to_proc[i] >= np_x*np_y*np_z*np_t) {
			lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_to_proc[%d]=%d (sbuf_mask={%d,%d,%d,%d}) out of range.\n", i,glattice.sbuf_to_proc[i],sbuf_mask[0],sbuf_mask[1],sbuf_mask[2],sbuf_mask[3]);
			test_q = false;
		}
		id = myid;
		for(k=0; k<4; k++) {
			if(shift[k] > 0) id = proc_down(id,k);
			else if(shift[k] < 0) id = proc_up(id,k);
		}
		if(id != glattice.sbuf_to_proc[i]) {
			lprintf("TEST_GEOMETRY.C",loglevel,"sbuf_to_proc[%d]=%d (sbuf_mask={%d,%d,%d,%d}) , expected %d.\n", i,glattice.sbuf_to_proc[i],sbuf_mask[0],sbuf_mask[1],sbuf_mask[2],sbuf_mask[3],id);
			test_q = false;
		}


#undef LOCAL		
		
	}
	error(!test_q,1,"test_geometry.c","Members nbuffers,rbuf*,sbuf* describes the sending/receiving buffers... FAILED");
	lprintf("TEST_GEOMETRY.C",loglevel,"Members nbuffers,rbuf*,sbuf* describes the sending/receiving buffers... OK\n");

}


void test_geometry_mpi() {

	initialize_test();


	/* TEST: ipt, iup, idown, glattice */
	
	lprintf("TEST_GEOMETRY.C",loglevel,"\nTESTING ipt, iup, idown, glattice\n");
	test_glattice();
	
	
	finalize_test();
}


void test_geometry_mpi_eo() {

	initialize_test();


	/* TEST: ipt, iup, idown, glattice */
	
	lprintf("TEST_GEOMETRY.C",loglevel,"\nTESTING ipt, iup, idown, glattice\n");
	test_glattice();


	/* TEST: glat_even */
	
	lprintf("TEST_GEOMETRY.C",loglevel,"\nTESTING glat_even\n");
	test_geometry_descriptor(&glat_even, &even_q);
	
	
	/* TEST: glat_odd */
	
	lprintf("TEST_GEOMETRY.C",loglevel,"\nTESTING glat_odd\n");
	test_geometry_descriptor(&glat_odd, &odd_q);

	
	finalize_test();
}
