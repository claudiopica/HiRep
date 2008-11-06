/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File geometry_init.c
 *
 * Inizialization of geometry structures
 *
 *******************************************************************************/

#include "geometry.h"
#include "global.h"
#include "error.h"
#include "logger.h"
#include <stdlib.h>

static int *alloc_mem=NULL;

static void free_memory() {
  if(alloc_mem!=NULL) {
    free(alloc_mem);
    alloc_mem=NULL;
    iup=idn=NULL;
    ipt=NULL;
  }
}

/* compute the size of the local lattice in the given direction given 
 * the cartesian coordinate of the node, global lattice size and the 
 * number of processes on that coordinate */
/* can return 0 if the npx is too big for the requested lattice */
static int compute_dim(int cx, int glb_x, int npx) {
  int n=glb_x/npx;
  int r=glb_x-n*npx; /* gbl_x%npx */

  /* do some sanity checks... */
  /* require a minimum size of 2*BORDERSIZE sites on all nodes */
  if (n<=(2*BORDERSIZE) && npx>1) { 
    error(1,1,"compute_dim " __FILE__,"The local lattice size is too small!!!");
  }


  /* return dimension */
  if (cx<r) ++n;
  return n;

}

/* this function return the coordinates of the local site (0,0,0,0)
*/
void origin_coord(int *c) {
  int d,r;

  d=GLB_T/NP_T; r=GLB_T-d*NP_T;
  c[0]=(COORD[0]<r)?(d+1)*COORD[0]:d*COORD[0]+r;
  d=GLB_X/NP_X; r=GLB_X-d*NP_X;
  c[1]=(COORD[1]<r)?(d+1)*COORD[1]:d*COORD[1]+r;
  d=GLB_Y/NP_Y; r=GLB_Y-d*NP_Y;
  c[2]=(COORD[2]<r)?(d+1)*COORD[2]:d*COORD[2]+r;
  d=GLB_Z/NP_Z; r=GLB_Z-d*NP_Z;
  c[3]=(COORD[3]<r)?(d+1)*COORD[3]:d*COORD[3]+r;
}

/* this function compute the parity of this process.
 * To do this we need to compute the global position of local site (0,0,0,0).
 * It is assumed that local dimensions are assigned by the above function compute_dim.
 *
 * INPUTS:
 * COORD[4], GLB_T, ..., GLB_Z, NP_T, ..., NP_Z
 * OUTPUTS:
 * PSIGN
 */
static void compute_psign() {
  int c[4];

  origin_coord(c);
  PSIGN=(c[0]+c[1]+c[2]+c[3])&1;
}


/* given the global coordinate of a site, this function returns the cartesian coordinates 
 * of the process on which the site is
 */
void glb_to_proc(int *g, int *p) {
  int d,r,c;

  /* sanity checks */
  if (!((g[0]<GLB_T)&&(g[1]<GLB_X)&&(g[2]<GLB_Y)&&(g[3]<GLB_Z))  
      || ((g[0]<0)||(g[1]<0)||(g[2]<0)||(g[3]<0))) { 
    error(1,1,"glb_to_proc " __FILE__,"The global coordinates are invalid!!!");
  }

  d=GLB_T/NP_T; r=GLB_T-d*NP_T;
  p[0]=0; c=-1; while(1) {c+=d; if (p[0]<r) c++; if(g[0]>c) p[0]++; else break; }
  d=GLB_X/NP_X; r=GLB_X-d*NP_X;
  p[1]=0; c=-1; while(1) {c+=d; if (p[1]<r) c++; if(g[1]>c) p[1]++; else break; }
  d=GLB_Y/NP_Y; r=GLB_Y-d*NP_Y;
  p[2]=0; c=-1; while(1) {c+=d; if (p[2]<r) c++; if(g[2]>c) p[2]++; else break; }
  d=GLB_Z/NP_Z; r=GLB_Z-d*NP_Z;
  p[3]=0; c=-1; while(1) {c+=d; if (p[3]<r) c++; if(g[3]>c) p[3]++; else break; }
}


/* this function set up the cartesian communicator */
/* return codes:
 * 0 => success
 * 1 => there is nothing to be done on this process node. 
 *      The caller should wait for other processes to finish.
 * 2 => error during inizialization
 *
 * INPUTS:
 * NP_T, NP_X, NP_Y, NP_Z
 * GLB_T, GLB_X, GLB_Y, GLB_Z
 *
 * OUTPUTS:
 * T, X, Y, Z, VOL3, VOLUME
 * MPI: cart_comm, CID, COORDS, CART_SIZE, PSIGN
 */
int geometry_init() {

#ifdef WITH_MPI
  /* MPI variables */
  int dims[4];
  int periods[4]={1,1,1,1}; /* all directions are periodic */
  int reorder=1; /* reassign ids */
  int mpiret;
#endif

#ifdef WITH_MPI
  MPI_Initialized(&mpiret);
  if(!mpiret) {
    lprintf("MPI",0,"ERROR: MPI has not been initialized!!!\n");
    error(1,1,"setup_process " __FILE__,"Cannot create cartesian communicator");
  }

  /* create the cartesian communicator */
  if(NP_T<2 && NP_X<2 && NP_Y<2 && NP_Z<2) {
    lprintf("GEOMETRY",0,"WARNING: NO PARALLEL DIMENSIONS SPECIFIED!!!\n");
    lprintf("GEOMETRY",0,"WARNING: THE MPI CODE SHOULD NOT BE USED IN THIS CASE!!!\n");
  }
  dims[0]=NP_T;
  dims[1]=NP_X;
  dims[2]=NP_Y;
  dims[3]=NP_Z;
  mpiret=MPI_Cart_create(MPI_COMM_WORLD, 4, dims, periods, reorder, &cart_comm);
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"setup_process " __FILE__,"Cannot create cartesian communicator");
  }
  if (cart_comm == MPI_COMM_NULL) {
    lprintf("GEOMETRY",0,"WARNING: Nothing to be done for PID %d!\n",PID);
    return 1;
  }
  MPI_Comm_size(cart_comm, &CART_SIZE);
  MPI_Comm_rank(cart_comm, &CID);
  MPI_Cart_coords(cart_comm, CID, 4, COORD);
  lprintf("GEOMETRY",0,"CId = %d (%d,%d,%d,%d)\n",CID,COORD[0],COORD[1],COORD[2],COORD[3]);
  if (CART_SIZE!=WORLD_SIZE)
    lprintf("GEOMETRY",0,"WARNING: Cartesian size %d != world size %d\n",CART_SIZE,WORLD_SIZE);
  /* this is for the parity of this process */
  compute_psign();
#else
  CID=0;
  COORD[0]=COORD[1]=COORD[2]=COORD[3]=0;
  NP_T=NP_X=NP_Y=NP_Z=1;
  PSIGN=0;
#endif


  /* from here on there are no more free parameters to set */

  /* compute local lattice size */
  T=compute_dim(COORD[0],GLB_T,NP_T);
  X=compute_dim(COORD[1],GLB_X,NP_X);
  Y=compute_dim(COORD[2],GLB_Y,NP_Y);
  Z=compute_dim(COORD[3],GLB_Z,NP_Z);

  VOL3=X*Y*Z;
  VOLUME=VOL3*T;

  return 0;
}

/* setup_process
 * Assign a unique PID to each process and setup 
 * communications if necessary
 *
 * return codes:
 * 0 => success
 *
 * OUTPUS:
 * MPI: PID, WORLD_SIZE
 */
int setup_process(int *argc, char ***argv) {

#ifdef WITH_MPI
  /* MPI variables */
  int mpiret;
#endif

#ifdef WITH_MPI
  /* init MPI */
  mpiret=MPI_Init(argc,argv);
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"setup_process " __FILE__,"MPI inizialization failed");
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &PID);
  MPI_Comm_size(MPI_COMM_WORLD, &WORLD_SIZE);
#else
  PID=0;
  WORLD_SIZE=1;
#endif

  return 0;

}


/* this function is intended to clean up before process ending
 *
 * return codes:
 * 0 => success
 */
int finalize_process() { 
#ifdef WITH_MPI
  /* MPI variables */
  int init;
#endif

#ifdef WITH_MPI
  MPI_Initialized(&init);
  if(init) MPI_Finalize();
#endif

  return 0;

}

void geometry_mem_alloc() {
  static int init=1;

  if (init) {
    int *cur;
    size_t req_mem=0;
    unsigned int VOL_SIZE=glattice.gsize;

    req_mem+=2*4*VOL_SIZE; /* for iup and idn */
    req_mem+=(X+2*X_BORDER)*(Y+2*Y_BORDER)*(Z+2*Z_BORDER)*(T+2*T_BORDER);     /* for ipt */

    alloc_mem=malloc(req_mem*sizeof(int));
    error((alloc_mem==NULL),1,"geometry_init [geometry_init.c]",
	"Cannot allocate memory");

    cur=alloc_mem;
#define ALLOC(ptr,size) ptr=cur; cur+=(size) 

    /* iup and idn */
    ALLOC(iup,4*VOL_SIZE);
    ALLOC(idn,4*VOL_SIZE);
    /* ipt */

    ALLOC(ipt,(X+2*X_BORDER)*(Y+2*Y_BORDER)*(Z+2*Z_BORDER)*(T+2*T_BORDER));
    /* ipt_4d */


    atexit(&free_memory);

    init=0;
#undef ALLOC

  }
}


int proc_up(int id, int dir) 
{
#ifdef WITH_MPI
  int coords[4];
	int outid;
	
	MPI_Cart_coords(cart_comm, id, 4, coords);
	++coords[dir];
	MPI_Cart_rank(cart_comm, coords, &outid);
  
	return outid;
#else
	return 0;
#endif
}

int proc_dn(int id, int dir)
{
#ifdef WITH_MPI
  int coords[4];
	int outid;
	
	MPI_Cart_coords(cart_comm, id, 4, coords);
	--coords[dir];
	MPI_Cart_rank(cart_comm, coords, &outid);
  
	return outid;
#else
	return 0;
#endif
}
