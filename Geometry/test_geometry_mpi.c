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
#include "memory.h"

#define true 1
#define false 0


static int loglevel = 0;


/* glattice is the geometry_descriptor instance fo the global lattice */
static int local_size[4];
static int buffer_thickness[4];
static int periodic_q[4];
static int nindices;
static int ncopies;
static int* copy_from;
static int* copy_to;

/*
static int** global_coord;
static int** rec_global_coord;
*/

#define NOT_ASSIGNED 0

#define ORIGINAL     1
#define DUPLICATE    2

#define INNER    1
#define LBORDER  2
#define RBORDER  3
#define LBUFFER  4
#define RBUFFER  5
#define LOCAL    6

#define EVEN 1
#define ODD  2

struct _site_info;
typedef struct _site_info {
  int index;
  int coord[4];                    /* local coordinates */
  int glb_coord[4];                /* global coordinates */
  unsigned int parity;             /* NOT_ASSIGNED ; EVEN ; ODD */
  unsigned int c_type;             /* NOT_ASSIGNED ; ORIGINAL; DUPLICATE */
  unsigned int local;              /* false ; true */
  struct _site_info* original;     /* if it is a DUPLICATE: pointer to its original */
  unsigned int ncopies;            /* if it is an ORIGINAL: number of copies */
  struct _site_info** copies;      /* if it is an ORIGINAL: list of copies */
  unsigned int b_type[4];          /* NOT_ASSIGNED ; INNER ; LBORDER ; RBORDER ; LBUFFER ; RBUFFER */
  int test;
} site_info;
static site_info* sites=NULL;
static site_info* origin=NULL;


struct _block_info;
typedef struct _block_info {
  int index;
  unsigned int mask[4];        /* NOT_ASSIGNED ; INNER ; LBORDER ; RBORDER ; LBUFFER ; RBUFFER ; LOCAL */  
  unsigned int parity;         /* NOT_ASSIGNED ; EVEN ; ODD */
  unsigned int length;
  site_info* start;
  unsigned int ninners;        /* how many directions are labeled as INNER */
  unsigned int nlocals;        /* how many directions are labeled as LOCAL */
  unsigned int nborders;       /* how many directions are labeled as L/RBORDER */
  unsigned int nbuffers;       /* how many directions are labeled as L/RBUFFER */
  unsigned int nsenders;
  struct _block_info** senders;    /* communication possible senders */
  struct _block_info* receiver;    /* communication receiver */
  unsigned int test;
} block_info;
static unsigned int nblocks;
static block_info* blocks=NULL;


typedef struct {
  block_info* from;
  block_info* to;
} communication_info;

/*
static const char* site_status(int x) {
  static const char label[4][256]={"NOT_ASSIGNED","ORIGINAL","DUPLICATE","NULL"};
  if(descr==NULL) return label[3];
  return label[descr[x].c_type];
}
*/

static void print_coordinates(site_info* s) {
  lprintf("TEST_GEOMETRY",loglevel," [S%d] Coordinates: ",s->index);
  if(s->c_type == NOT_ASSIGNED)
    lprintf("TEST_GEOMETRY",loglevel,"NOT_ASSIGNED");
  else if(s->c_type == ORIGINAL)
    lprintf("TEST_GEOMETRY",loglevel,"ORIGINAL");
  else if(s->c_type == DUPLICATE)
    lprintf("TEST_GEOMETRY",loglevel,"DUPLICATE");
  else
    lprintf("TEST_GEOMETRY",loglevel,"INVALID");
  lprintf("TEST_GEOMETRY",loglevel," ( %d , %d , %d , %d )\n",
    s->coord[0],s->coord[1],s->coord[2],s->coord[3]);
}

static void print_global_coordinates(site_info* s) {
  lprintf("TEST_GEOMETRY",loglevel," [S%d] Global coordinates: ",s->index);
  lprintf("TEST_GEOMETRY",loglevel,"( %d , %d , %d , %d )\n",
    s->glb_coord[0],s->glb_coord[1],s->glb_coord[2],s->glb_coord[3]);
}

static void print_parity(site_info* s) {
  lprintf("TEST_GEOMETRY",loglevel," [S%d] Parity ",s->index);
  if(s->parity == NOT_ASSIGNED)
    lprintf("TEST_GEOMETRY",loglevel,"NOT_ASSIGNED\n");
  else if(s->parity == EVEN)
    lprintf("TEST_GEOMETRY",loglevel,"EVEN\n");
  else if(s->parity == ODD)
    lprintf("TEST_GEOMETRY",loglevel,"ODD\n");
  else
    lprintf("TEST_GEOMETRY",loglevel,"INVALID\n");
}

static void print_copies_info(site_info* s) {
  int k;
  
  lprintf("TEST_GEOMETRY",loglevel," [S%d] Its original: %d\n",s->index,s->original->index);
  if(s->original->index != s->index) print_coordinates(s->original);
  
  lprintf("TEST_GEOMETRY",loglevel," [S%d] Number of copies: %d\n",s->index,s->ncopies);
  
  if(s->ncopies>0) {
    lprintf("TEST_GEOMETRY",loglevel," [S%d] Copies:",s->index);
    for(k=0; k<s->ncopies; k++)
      lprintf("TEST_GEOMETRY",loglevel," %d",s->copies[k]->index);
    lprintf("TEST_GEOMETRY",loglevel,"\n");
    
    for(k=0; k<s->ncopies; k++)
      print_coordinates(s->copies[k]);
  }
}


static void print_geometric_info(site_info* s) {
  int i;
  
  lprintf("TEST_GEOMETRY",loglevel," [S%d] Location:",s->index);
  for(i=0; i<4; i++) {
    if(s->b_type[i] == LBUFFER)
      lprintf("TEST_GEOMETRY",loglevel," LBUFFER");
    else if(s->b_type[i] == LBORDER)
      lprintf("TEST_GEOMETRY",loglevel," LBORDER");
    else if(s->b_type[i] == INNER)
      lprintf("TEST_GEOMETRY",loglevel," INNER");
    else if(s->b_type[i] == RBORDER)
      lprintf("TEST_GEOMETRY",loglevel," RBORDER");
    else if(s->b_type[i] == RBUFFER)
      lprintf("TEST_GEOMETRY",loglevel," RBUFFER");
    else if(s->b_type[i] == NOT_ASSIGNED)
      lprintf("TEST_GEOMETRY",loglevel," NOT_ASSIGNED");
    else
      lprintf("TEST_GEOMETRY",loglevel," INVALID");
  }
  lprintf("TEST_GEOMETRY",loglevel,"\n");
  if(s->local) lprintf("TEST_GEOMETRY",loglevel," [S%d] Local",s->index);
  else lprintf("TEST_GEOMETRY",loglevel," [S%d] Nonlocal",s->index);
  lprintf("TEST_GEOMETRY",loglevel,"\n");
}


static void print_full_site_info(site_info* s) {
//  lprintf("TEST_GEOMETRY",loglevel,"Full info for site %d\n",s->index);
  
  print_coordinates(s);
  
  print_global_coordinates(s);
  
  print_parity(s);
  
  print_copies_info(s);
  
  print_geometric_info(s);

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

static void set_originals_coords(site_info* s, int *test_q) {
  int i, nb;
  int c[4];

  if(s->c_type == NOT_ASSIGNED) {
    lprintf("TEST_GEOMETRY",loglevel,"ERROR set_originals_coords: coordinates not assiged yet\n");
    print_full_site_info(s);
    *test_q = false;
    exit(-1);
  }
  
  for(i=0; i<4; i++) {
    nb = iup(s->index,i);
    memcpy(c,s->coord,sizeof(int)*4);
    c[i]++;
    if(periodic_q[i]) c[i] = safe_mod(c[i],local_size[i]);
    if(!in_glattice_q(c)) continue;
    if(nb >= 0 && nb < nindices) {
      if(sites[nb].c_type == NOT_ASSIGNED) {
        sites[nb].c_type = ORIGINAL;
        memcpy(sites[nb].coord,c,sizeof(int)*4);
        set_originals_coords(&sites[nb], test_q);
      }
    } else {
      lprintf("TEST_GEOMETRY",0,"Bad candidate index %d for coordinates (%d,%d,%d,%d), reached from (%d,%d,%d,%d) with iup[%d].\n",
              nb,c[0],c[1],c[2],c[3],s->coord[0],s->coord[1],s->coord[2],s->coord[3],i);
      *test_q = false;
    }
  }
  
  for(i=0; i<4; i++) {
    nb = idn(s->index,i);
    memcpy(c,s->coord,sizeof(int)*4);
    c[i]--;
    if(periodic_q[i]) c[i] = safe_mod(c[i],local_size[i]);
    if(!in_glattice_q(c)) continue;
    if(nb >= 0 && nb < nindices) {
      if(sites[nb].c_type == NOT_ASSIGNED) {
        sites[nb].c_type = ORIGINAL;
        memcpy(sites[nb].coord,c,sizeof(int)*4);
        set_originals_coords(&sites[nb], test_q);
      }
    } else {
      lprintf("TEST_GEOMETRY",0,"Bad candidate index %d for coordinates (%d,%d,%d,%d), reached from (%d,%d,%d,%d) with idn[%d].\n",
              nb,c[0],c[1],c[2],c[3],s->coord[0],s->coord[1],s->coord[2],s->coord[3],i);
      *test_q = false;
    }
  }

  s->parity = ((s->coord[0]+s->coord[1]+s->coord[2]+s->coord[3]+PSIGN+1)&1) ? EVEN : ODD;
  
  s->local = true;
  for(i=0; i<4; i++) {
    if(s->coord[i]<0) {
      s->b_type[i] = LBUFFER;
      s->local = false;
    } else if(s->coord[i]<buffer_thickness[i]) {
      s->b_type[i] = LBORDER;
    } else if(s->coord[i]<local_size[i]-buffer_thickness[i]) {
      s->b_type[i] = INNER;
    } else if(s->coord[i]<local_size[i]) {
      s->b_type[i] = RBORDER;
    } else {
      s->b_type[i] = RBUFFER;
      s->local = false;
    }
  }

  s->glb_coord[0] = (zerocoord[0]+s->coord[0] + GLB_T)%GLB_T;
  s->glb_coord[1] = (zerocoord[1]+s->coord[1] + GLB_X)%GLB_X;
  s->glb_coord[2] = (zerocoord[2]+s->coord[2] + GLB_Y)%GLB_Y;
  s->glb_coord[3] = (zerocoord[3]+s->coord[3] + GLB_Z)%GLB_Z;
}




static void set_duplicates_coords(int* test_q) {
  int i, k;
  site_info *from;
  site_info *to;
  site_info** tmp;
  
  for(i=0; i<ncopies; i++) {
    from = &sites[copy_from[i]];
    to = &sites[copy_to[i]];
    
    if(to->c_type != NOT_ASSIGNED) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR set_duplicates_coords: attempting copy from %d to %d\n",copy_from[i],copy_to[i]);
      lprintf("TEST_GEOMETRY",loglevel,"ERROR set_duplicates_coords: but destination has been already initialized\n");
      print_full_site_info(from);
      print_full_site_info(to);
      *test_q = false;
      exit(-1);
    }
    if(from->c_type != ORIGINAL) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR set_duplicates_coords: attempting copy from %d to %d\n",copy_from[i],copy_to[i]);
      lprintf("TEST_GEOMETRY",loglevel,"ERROR set_duplicates_coords: but source is not an ORIGINAL\n");
      print_full_site_info(from);
      print_full_site_info(to);
      *test_q = false;
      exit(-1);
    }
    
    for(k=0; k<4; k++) {
      to->coord[k] = from->coord[k];
      to->b_type[k] = from->b_type[k];
      to->glb_coord[k] = from->glb_coord[k];
    }
    to->parity = from->parity;
    to->original = from;
    to->c_type = DUPLICATE;
    
    tmp = from->copies;
    from->ncopies++;
    from->copies = malloc(sizeof(site_info*)*from->ncopies);
    if(tmp!=NULL) {
      memcpy(from->copies,tmp,sizeof(site_info*)*(from->ncopies-1));
      free(tmp);
    }
    from->copies[from->ncopies-1] = to;    
  }
}




static void initialize_sites() {
  int i, j, k, test_q;
  site_info *s;
  
  local_size[0] = T; local_size[1] = X; local_size[2] = Y; local_size[3] = Z;
  buffer_thickness[0] = T_BORDER; buffer_thickness[1] = X_BORDER; buffer_thickness[2] = Y_BORDER; buffer_thickness[3] = Z_BORDER;
  periodic_q[0] = (NP_T==1) ? true : false; periodic_q[1] = (NP_X==1) ? true : false; periodic_q[2] = (NP_Y==1) ? true : false; periodic_q[3] = (NP_Z==1) ? true : false; 
  
  nindices = glattice.gsize_gauge;
  
  ncopies = 0;
  for(i=0 ; i<glattice.ncopies_gauge ; i++)
    ncopies += glattice.copy_len[i];
  
  copy_from = copy_to = NULL;
  if(ncopies != 0) {
    copy_from = malloc(sizeof(int)*ncopies);
    copy_to = malloc(sizeof(int)*ncopies);
    j=0;
    for(i=0 ; i<glattice.ncopies_gauge ; i++) {
      for(k=0; k<glattice.copy_len[i]; k++) {
        copy_from[j] = glattice.copy_from[i]+k;
        copy_to[j] = glattice.copy_to[i]+k;
        j++;
      }
    }
  }


  sites = (site_info*)malloc(sizeof(site_info)*nindices);
  for(i=0; i<nindices; i++) {
    sites[i].index = i;
    sites[i].parity = NOT_ASSIGNED;
    sites[i].c_type = NOT_ASSIGNED;
    sites[i].original = &sites[i];
    sites[i].ncopies = 0;
    sites[i].copies = NULL;
    sites[i].b_type[0] = NOT_ASSIGNED;
    sites[i].b_type[1] = NOT_ASSIGNED;
    sites[i].b_type[2] = NOT_ASSIGNED;
    sites[i].b_type[3] = NOT_ASSIGNED;
    sites[i].local = false;
  }

  lprintf("TEST_GEOMETRY",loglevel,"Memory allocation... OK\n");
  lprintf("TEST_GEOMETRY",loglevel,"nindices = %d\n", nindices);
  lprintf("TEST_GEOMETRY",loglevel,"ncopies = %d\n", ncopies);
  lprintf("TEST_GEOMETRY",loglevel,"PSIGN = %d\n", PSIGN);
  lprintf("TEST_GEOMETRY",loglevel,"local_size = {%d,%d,%d,%d}\n", local_size[0], local_size[1], local_size[2], local_size[3]);
  lprintf("TEST_GEOMETRY",loglevel,"buffer_thickness = {%d,%d,%d,%d}\n", buffer_thickness[0], buffer_thickness[1], buffer_thickness[2], buffer_thickness[3]);
  lprintf("TEST_GEOMETRY",loglevel,"periodic_q = {%d,%d,%d,%d}\n", periodic_q[0], periodic_q[1], periodic_q[2], periodic_q[3]);
  lprintf("TEST_GEOMETRY",loglevel,"sites = %p ... %p\n", sites, sites+nindices-1);
  lprintf("TEST_GEOMETRY",loglevel,"iup = %p ... %p\n", iup, iup+4*nindices-1);
  lprintf("TEST_GEOMETRY",loglevel,"idn = %p ... %p\n", idn, idn+4*nindices-1);
  lprintf("TEST_GEOMETRY",loglevel,"ipt = %p ... %p\n", ipt, ipt+VOLUME-1);


  /* Set the coordinates for the points reached by iup, idn - iteratively from the origin */
  origin = &sites[ipt(0,0,0,0)];
  origin->coord[0] = origin->coord[1] = origin->coord[2] = origin->coord[3] = 0;
  origin->c_type = ORIGINAL;
  test_q = true;
  set_originals_coords(origin, &test_q);

  error(!test_q,1,"test_geometry.c","Set the coordinates for points reached by iup, idn... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Set the coordinates for points reached by iup, idn... OK\n");
  
  
  /* Set the coordinates for the duplicate points */
  test_q = true;
  set_duplicates_coords(&test_q);
  error(!test_q,1,"test_geometry.c","Set the coordinates for duplicate points... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Set the coordinates for duplicate points... OK\n");  


  lprintf("TEST_GEOMETRY",loglevel,"Origin\n");
  print_full_site_info(origin);

  /* TEST: No duplicate coordinates for the original indices */
  test_q = true;
  for(s=sites; s<sites+nindices; s++) {
    if(s->c_type != ORIGINAL) continue;
    site_info *s2;
    for(s2=s+1; s2<sites+nindices; s2++) {
      if(s2->c_type != ORIGINAL) continue;
      if(s->coord[0]==s2->coord[0] && s->coord[1]==s2->coord[1] && s->coord[2]==s2->coord[2] && s->coord[3]==s2->coord[3]) {
        lprintf("TEST_GEOMETRY",0,"Original sites %d and %d have same coordinates\n",s->index,s2->index);
        print_full_site_info(s);
        print_full_site_info(s2);
        test_q = false;
      }
    }
  }
  error(!test_q,1,"test_geometry.c","No duplicate coordinates for the original indexes... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"No duplicate coordinates for the original indexes... OK\n");

  /* TEST: Coordinates are assigned to each point */
  test_q = true;
  for(s=sites; s<sites+nindices; s++) {
    if(s->c_type == NOT_ASSIGNED) {
      lprintf("TEST_GEOMETRY",0,"Coordinates are not assigned to the index %d\n",s->index);
      test_q = false;
    }
  }
  error(!test_q,1,"test_geometry.c","Coordinates are assigned to each point... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Coordinates are assigned to each point... OK\n");

  /* TEST: Coordinates are in the right range */
  test_q = true;
  for(s=sites; s<sites+nindices; s++) {
    if(s->coord[0] < -buffer_thickness[0] || s->coord[0] >= local_size[0]+buffer_thickness[0]
      || s->coord[1] < -buffer_thickness[1] || s->coord[1] >= local_size[1]+buffer_thickness[1]
      || s->coord[2] < -buffer_thickness[2] || s->coord[2] >= local_size[2]+buffer_thickness[2]
      || s->coord[3] < -buffer_thickness[3] || s->coord[3] >= local_size[3]+buffer_thickness[3]) {
      lprintf("TEST_GEOMETRY",0,"Coordinates of %d out of range\n",s->index);
      print_full_site_info(s);
      test_q = false;
    }
  }
  error(!test_q,1,"test_geometry.c","Coordinates are in the right range... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Coordinates are in the right range... OK\n");

  
  /* TEST: ipt gives the right map between coordinates and local indexes */
  test_q = true;
  
  for(s=sites; s<sites+nindices; s++) {
    if(!s->local || s->c_type != ORIGINAL) continue;
    int *cx = s->coord;
    if(ipt(cx[0],cx[1],cx[2],cx[3]) != s->index) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR ipt(%d,%d,%d,%d)=%d ; but s->index=%d .\n",cx[0],cx[1],cx[2],cx[3],ipt(cx[0],cx[1],cx[2],cx[3]),s->index);
      print_full_site_info(s);
      test_q = false;
    }
  }
  
  int c0, c1, c2, c3;
  for(c0=0; c0<local_size[0]; c0++)
  for(c1=0; c1<local_size[1]; c1++)
  for(c2=0; c2<local_size[2]; c2++)
  for(c3=0; c3<local_size[3]; c3++) {
    int x = ipt(c0,c1,c2,c3);
    if(x<0) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR ipt(%d,%d,%d,%d)=%d < 0\n",c0,c1,c2,c3,x);
      test_q = false;
    } else if(x>=nindices) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR ipt(%d,%d,%d,%d)=%d >= nindices\n",c0,c1,c2,c3,x);
      test_q = false;
    } else if(!sites[x].local) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR ipt(%d,%d,%d,%d)=%d which is nonlocal\n",c0,c1,c2,c3,x);
      print_full_site_info(sites+x);
      test_q = false;
    } else if(sites[x].c_type != ORIGINAL) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR ipt(%d,%d,%d,%d)=%d which is not ORIGINAL\n",c0,c1,c2,c3,x);
      print_full_site_info(sites+x);
      test_q = false;
    } else {
      int *cx = sites[x].coord;
      if(c0 != cx[0] || c1 != cx[1] || c2 != cx[2] || c3 != cx[3]) {
        lprintf("TEST_GEOMETRY",loglevel,"ipt(%d,%d,%d,%d)=%d but coord[%x]=(%d,%d,%d,%d).\n",c0,c1,c2,c3,x,x,cx[0],cx[1],cx[2],cx[3]);
        print_full_site_info(sites+x);
        test_q = false;
      }
    }
  }
  
  error(!test_q,1,"test_geometry.c","ipt gives the right map between coordinates and local indexes... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"ipt gives the right map between coordinates and local indexes... OK\n");
  
  /* TEST: iup, idn give the right maps between neighbours */
  int id[4][4];
  for(i=0; i<4; i++) {
    for(j=0; j<4; j++) id[i][j] = 0;
    id[i][i] = 1;
  }

  test_q = true;
  for(s=sites; s<sites+nindices; s++) {
    if(s->c_type != ORIGINAL) continue;
    int *cx = s->coord;
    for(i=0; i<4; i++) {
      int y = iup(s->index,i);
      int c[4];
      memcpy(c,s->coord,sizeof(int)*4);
      c[i]++;
      if(periodic_q[i]) c[i] = safe_mod(c[i],local_size[i]);
      if(!in_glattice_q(c)) continue;
      int *cy = sites[y].coord;
      if(cy[0] != cx[0]+id[i][0] &&
         cy[1] != cx[1]+id[i][1] &&
         cy[2] != cx[2]+id[i][2] &&
         cy[3] != cx[3]+id[i][3]) {
          lprintf("TEST_GEOMETRY",loglevel,"ERROR iup[%d](%d,%d,%d,%d)=(%d,%d,%d,%d)\n",i,cx[0],cx[1],cx[2],cx[3],cy[0],cy[1],cy[2],cy[3]);
          test_q = false;
      }
    }
    for(i=0; i<4; i++) {
      int y = idn(s->index,i);
      int c[4];
      memcpy(c,s->coord,sizeof(int)*4);
      c[i]--;
      if(periodic_q[i]) c[i] = safe_mod(c[i],local_size[i]);
      if(!in_glattice_q(c)) continue;
      int *cy = sites[y].coord;
      if(cy[0] != cx[0]-id[i][0] &&
         cy[1] != cx[1]-id[i][1] &&
         cy[2] != cx[2]-id[i][2] &&
         cy[3] != cx[3]-id[i][3]) {
          lprintf("TEST_GEOMETRY",loglevel,"ERROR idn[%d](%d,%d,%d,%d)=(%d,%d,%d,%d)\n",i,cx[0],cx[1],cx[2],cx[3],cy[0],cy[1],cy[2],cy[3]);
          test_q = false;
      }
    }
  }
  error(!test_q,1,"test_geometry.c","iup, idn give the right maps between neighbours... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"iup, idn give the right maps between neighbours... OK\n");
}



static void finalize_test() {
  int i;
  for(i=0; i<nindices; i++) {
    if(sites[i].copies != NULL) free(sites[i].copies);
  }
  free(sites);
}



static int compare_to_b_mask(const site_info *s, const unsigned int b_mask[4]) {
  if(b_mask[0]!=NOT_ASSIGNED && b_mask[0]!=LOCAL && b_mask[0]!=s->b_type[0]) return false;
  if(b_mask[1]!=NOT_ASSIGNED && b_mask[1]!=LOCAL && b_mask[1]!=s->b_type[1]) return false;
  if(b_mask[2]!=NOT_ASSIGNED && b_mask[2]!=LOCAL && b_mask[2]!=s->b_type[2]) return false;
  if(b_mask[3]!=NOT_ASSIGNED && b_mask[3]!=LOCAL && b_mask[3]!=s->b_type[3]) return false;
  if(b_mask[0]==LOCAL && s->b_type[0]!=INNER && s->b_type[0]!=LBORDER && s->b_type[0]!=RBORDER) return false;
  if(b_mask[1]==LOCAL && s->b_type[1]!=INNER && s->b_type[1]!=LBORDER && s->b_type[1]!=RBORDER) return false;
  if(b_mask[2]==LOCAL && s->b_type[2]!=INNER && s->b_type[2]!=LBORDER && s->b_type[2]!=RBORDER) return false;
  if(b_mask[3]==LOCAL && s->b_type[3]!=INNER && s->b_type[3]!=LBORDER && s->b_type[3]!=RBORDER) return false;
  return true;
}

static int compare_to_parity(const site_info *s, const int parity) {
  if(parity!=NOT_ASSIGNED && parity!=s->parity) return false;
  return true;
}


static void search_first_block(block_info *block) {
  site_info *s, *s2, *s3, *starting_site;
  int flag, n, k;
  unsigned int mask2[4]={block->mask[0],block->mask[1],block->mask[2],block->mask[3]},parity2;
  int loc[4]={T,X,Y,Z};
  
  block->ninners = 0;
  block->nlocals = 0;
  block->nbuffers = 0;
  block->nborders = 0;
  for(k=0;k<4;k++) {
    if(block->mask[k]==LOCAL) block->nlocals++;
    if(block->mask[k]==INNER) block->ninners++;
    if(block->mask[k]==LBORDER || block->mask[k]==RBORDER) block->nborders++;
    if(block->mask[k]==LBUFFER || block->mask[k]==RBUFFER) block->nbuffers++;
  }
  
  /* Calculate block size */
  parity2=block->parity;
  for(k=0;k<4;k++) {
    if(mask2[k]==LBUFFER){
      mask2[k]=RBORDER;
      if(loc[k]%2==1) parity2=EVEN+ODD-parity2;
    }
    if(mask2[k]==RBUFFER){
      mask2[k]=LBORDER;
      if(loc[k]%2==1) parity2=EVEN+ODD-parity2;
    }
  }

  block->length=0;
  for(s=sites; s<sites+nindices; s++) {
    if(compare_to_b_mask(s,mask2) && compare_to_parity(s,parity2) && s->c_type==ORIGINAL)
      block->length++;
  }
    
  if(block->length==0) {
    block->start=NULL;
    return;
  }
  
  /* Identify consecutive sites with desired parity and b_type */
  starting_site=block->start;
  if(starting_site==NULL) starting_site=sites;
  
  n=0;
  for(s=sites+nindices-1; s>=starting_site; s--) {
    if(compare_to_b_mask(s,block->mask) && compare_to_parity(s,block->parity))
      n++;
    else
      n=0;
    s->test=n;
  }
  
  for(s=starting_site; s<sites+nindices; s++) {
    if(s->test < block->length) continue;
    if(!compare_to_b_mask(s,block->mask) || !compare_to_parity(s,block->parity)) {
      lprintf("TEST_GEOMETRY",loglevel,"Searching block parity=%d mask=(%d,%d,%d,%d) length=%d\n",block->parity,block->mask[0],block->mask[1],block->mask[2],block->mask[3],block->length);
      lprintf("TEST_GEOMETRY",loglevel,"ERROR search_block: I shouldn't be here (1)");
      exit(-1);
    }
    flag=true;
    for(s2=s+1; s2<s+block->length; s2++) {
      if(!compare_to_b_mask(s2,block->mask) || !compare_to_parity(s2,block->parity)) {
        lprintf("TEST_GEOMETRY",loglevel,"Searching block parity=%d mask=(%d,%d,%d,%d) length=%d\n",block->parity,block->mask[0],block->mask[1],block->mask[2],block->mask[3],block->length);
        lprintf("TEST_GEOMETRY",loglevel,"ERROR search_block: I shouldn't be here (2)");
        exit(-1);
      }
    }
    /* Blocks must contain each point only once */
    for(s2=s; s2<s+block->length; s2++)
    for(s3=s2+1; s3<s+block->length; s3++) {
      int *c2 = s2->coord;
      int *c3 = s3->coord;
      if(c3[0]==c2[0] && c3[1]==c2[1] && c3[2]==c2[2] && c3[3]==c2[3]) {
        flag=false;
        break;
      }
    }
    if(flag) {
      block->start=s;
      return;
    }
  }

  block->start=NULL;
}


static void append_block(block_info* nb) {
  block_info *tmp, *b;
  
  tmp = malloc(sizeof(block_info)*(nblocks+1));
  if(nblocks>0) {
    memcpy(tmp,blocks,sizeof(block_info)*nblocks);
    free(blocks);
  }
  blocks = tmp;
  b = blocks+nblocks;

  b->index = nblocks;
  b->ninners = nb->ninners;
  b->nlocals = nb->nlocals;
  b->nbuffers = nb->nbuffers;
  b->nborders = nb->nborders;
  b->mask[0] = nb->mask[0];
  b->mask[1] = nb->mask[1];
  b->mask[2] = nb->mask[2];
  b->mask[3] = nb->mask[3];
  b->parity = nb->parity;
  b->length = nb->length;
  b->start = nb->start;
  b->nsenders = nb->nsenders;
  b->senders = nb->senders;
  b->receiver = nb->receiver;

  nblocks++;
}


static void print_short_block_info(block_info *b) {
  int i;
  
  lprintf("TEST_GEOMETRY",loglevel," [B%d]  \tStart %d  \tLength %d", b->index, b->start->index, b->length);
  
  lprintf("TEST_GEOMETRY",loglevel,"  \tFROM",b->index);
  if(b->nsenders>0) {
    for(i=0; i<b->nsenders; i++)
      lprintf("TEST_GEOMETRY",loglevel," %d",b->senders[i]->index);
  } else
    lprintf("TEST_GEOMETRY",loglevel," *  ");
  
  lprintf("TEST_GEOMETRY",loglevel,"  \tTO",b->index);
  if(b->receiver!=NULL) {
    lprintf("TEST_GEOMETRY",loglevel," %d",b->receiver->index);
  } else
    lprintf("TEST_GEOMETRY",loglevel," *  ");

  lprintf("TEST_GEOMETRY",loglevel,"     \tMask (");
  for(i=0; i<4; i++) {
    if(b->mask[i] == LBUFFER)
      lprintf("TEST_GEOMETRY",loglevel," LBUFFER");
    else if(b->mask[i] == LBORDER)
      lprintf("TEST_GEOMETRY",loglevel," LBORDER");
    else if(b->mask[i] == INNER)
      lprintf("TEST_GEOMETRY",loglevel," INNER");
    else if(b->mask[i] == RBORDER)
      lprintf("TEST_GEOMETRY",loglevel," RBORDER");
    else if(b->mask[i] == RBUFFER)
      lprintf("TEST_GEOMETRY",loglevel," RBUFFER");
    else if(b->mask[i] == LOCAL)
      lprintf("TEST_GEOMETRY",loglevel," LOCAL");
    else if(b->mask[i] == NOT_ASSIGNED)
      lprintf("TEST_GEOMETRY",loglevel," NOT_ASSIGNED");
    else
      lprintf("TEST_GEOMETRY",loglevel," INVALID");
  }
  
  lprintf("TEST_GEOMETRY",loglevel," )  \tParity ");
  if(b->parity == NOT_ASSIGNED)
    lprintf("TEST_GEOMETRY",loglevel,"NOT_ASSIGNED");
  else if(b->parity == EVEN)
    lprintf("TEST_GEOMETRY",loglevel,"EVEN");
  else if(b->parity == ODD)
    lprintf("TEST_GEOMETRY",loglevel,"ODD");
  else
    lprintf("TEST_GEOMETRY",loglevel,"INVALID");
  
  lprintf("TEST_GEOMETRY",loglevel,"\n");
}


static void print_block_info(block_info *b) {
  int i;
  
  lprintf("TEST_GEOMETRY",loglevel," [B%d] block info\n",b->index);

  lprintf("TEST_GEOMETRY",loglevel," [B%d] mask:",b->index,b->index);
  for(i=0; i<4; i++) {
    if(b->mask[i] == LBUFFER)
      lprintf("TEST_GEOMETRY",loglevel," LBUFFER");
    else if(b->mask[i] == LBORDER)
      lprintf("TEST_GEOMETRY",loglevel," LBORDER");
    else if(b->mask[i] == INNER)
      lprintf("TEST_GEOMETRY",loglevel," INNER");
    else if(b->mask[i] == RBORDER)
      lprintf("TEST_GEOMETRY",loglevel," RBORDER");
    else if(b->mask[i] == RBUFFER)
      lprintf("TEST_GEOMETRY",loglevel," RBUFFER");
    else if(b->mask[i] == LOCAL)
      lprintf("TEST_GEOMETRY",loglevel," LOCAL");
    else if(b->mask[i] == NOT_ASSIGNED)
      lprintf("TEST_GEOMETRY",loglevel," NOT_ASSIGNED");
    else
      lprintf("TEST_GEOMETRY",loglevel," INVALID");
  }
  lprintf("TEST_GEOMETRY",loglevel,"\n");
  
  lprintf("TEST_GEOMETRY",loglevel," [B%d] ninners = %d\n",b->index,b->ninners);
  lprintf("TEST_GEOMETRY",loglevel," [B%d] nlocals = %d\n",b->index,b->nlocals);
  lprintf("TEST_GEOMETRY",loglevel," [B%d] nborders = %d\n",b->index,b->nborders);
  lprintf("TEST_GEOMETRY",loglevel," [B%d] nbuffers = %d\n",b->index,b->nbuffers);

  lprintf("TEST_GEOMETRY",loglevel," [B%d] Parity: ",b->index);
  if(b->parity == NOT_ASSIGNED)
    lprintf("TEST_GEOMETRY",loglevel,"NOT_ASSIGNED\n");
  else if(b->parity == EVEN)
    lprintf("TEST_GEOMETRY",loglevel,"EVEN\n");
  else if(b->parity == ODD)
    lprintf("TEST_GEOMETRY",loglevel,"ODD\n");
  else
    lprintf("TEST_GEOMETRY",loglevel,"INVALID\n");

  lprintf("TEST_GEOMETRY",loglevel," [B%d] Length: %d\n",b->index,b->length);
  if(b->start==NULL) {
    lprintf("TEST_GEOMETRY",loglevel," [B%d] First index: NULL\n",b->index);
  } else {
    lprintf("TEST_GEOMETRY",loglevel," [B%d] First index: %d\n",b->index,b->start->index);
    print_full_site_info(b->start);
  }
  
  
  if(b->nsenders>0) {
    lprintf("TEST_GEOMETRY",loglevel," [B%d] Senders:",b->index);
    for(i=0; i<b->nsenders; i++)
      lprintf("TEST_GEOMETRY",loglevel," %d",b->senders[i]->index);
    lprintf("TEST_GEOMETRY",loglevel,"\n");
  }
  
  if(b->receiver!=NULL) {
    lprintf("TEST_GEOMETRY",loglevel," [B%d] Receiver: %d\n",b->index,b->receiver->index);
  }
}


static void search_all_blocks(unsigned int parity, unsigned int mask[4]) {
  block_info tmp;
  int n;
  
  tmp.index = 0;
  tmp.ninners = 0;
  tmp.nlocals = 0;
  tmp.nbuffers = 0;
  tmp.nborders = 0;
  tmp.mask[0] = mask[0];
  tmp.mask[1] = mask[1];
  tmp.mask[2] = mask[2];
  tmp.mask[3] = mask[3];
  tmp.parity = parity;
  tmp.length = 0;
  tmp.start = NULL;
  tmp.nsenders = 0;
  tmp.senders = NULL;
  tmp.receiver = NULL;
  
  n=0;
  while(1) {
    search_first_block(&tmp);
    if(tmp.start!=NULL) {
      append_block(&tmp);
      tmp.start++;
      n++;
    } else
      break;
  }
  
  /*
  if(n==0) {
    lprintf("TEST_GEOMETRY",loglevel,"ERROR search_all_blocks: Block not found");
    print_block_info(&tmp);
    error(1,1,"test_geometry.c","Block initialization... FAILED");
  }
  */
}



static void initialize_blocks() {
  int n, j, k, bj, bk, test_q;
  unsigned int mask[4], parity;
  site_info *s;
  block_info *b;
  
  nblocks=0;
  
  mask[0]=mask[1]=mask[2]=mask[3]=INNER;
  parity=EVEN;
  search_all_blocks(parity, mask);

  mask[0]=mask[1]=mask[2]=mask[3]=INNER;
  parity=ODD;
  search_all_blocks(parity, mask);

/*
#define INNER    1
#define LBORDER  2
#define RBORDER  3
#define LBUFFER  4
#define RBUFFER  5

#define EVEN 1
#define ODD  2
*/
  
  for(j=0; j<4; j++) {
    if(periodic_q[j]) continue;
    for(bj=LBORDER; bj<=RBUFFER; bj++)
    for(parity=EVEN; parity<=ODD; parity++) {
      mask[0]=mask[1]=mask[2]=mask[3]=LOCAL;
      mask[j]=bj;
      search_all_blocks(parity, mask);
    }
  }
  
  for(j=0; j<4; j++)
  for(k=j+1; k<4; k++) {
    if(periodic_q[j] || periodic_q[k]) continue;
    for(bj=LBORDER; bj<=RBUFFER; bj++)
    for(bk=LBORDER; bk<=RBUFFER; bk++)
    for(parity=EVEN; parity<=ODD; parity++) {
      if((bj==LBORDER || bj==RBORDER) && (bk==LBUFFER || bk==RBUFFER)) continue;
      if((bk==LBORDER || bk==RBORDER) && (bj==LBUFFER || bj==RBUFFER)) continue;
      mask[0]=mask[1]=mask[2]=mask[3]=LOCAL;
      mask[j]=bj;
      mask[k]=bk;
      search_all_blocks(parity, mask);
    }
  }
  
  /* Find communication pairing */
  for(b=blocks;b<blocks+nblocks;b++) {
    if(b->nbuffers==0) continue;
    int twinmask[4]={b->mask[0],b->mask[1],b->mask[2],b->mask[3]};
    int twinparity=b->parity;
    for(k=0;k<4;k++) {
      if(twinmask[k]==LBUFFER){
        twinmask[k]=RBORDER;
        if(local_size[k]%2==1) twinparity=EVEN+ODD-twinparity;
      }
      if(twinmask[k]==RBUFFER){
        twinmask[k]=LBORDER;
        if(local_size[k]%2==1) twinparity=EVEN+ODD-twinparity;
      }
    }
    block_info* b2;
    for(b2=blocks;b2<blocks+nblocks;b2++) {
      /* Find sender with the right parity, mask and length*/
      if(b2->mask[0]!=twinmask[0] || b2->mask[1]!=twinmask[1] || b2->mask[2]!=twinmask[2] || b2->mask[3]!=twinmask[3] || b2->parity!=twinparity || b2->length!=b->length)
        continue;
      
      /* Sender must be obtained by shifting the original block */
      int shift[4];
      shift[0]=b2->start[0].coord[0]-b->start[0].coord[0];
      shift[1]=b2->start[0].coord[1]-b->start[0].coord[1];
      shift[2]=b2->start[0].coord[2]-b->start[0].coord[2];
      shift[3]=b2->start[0].coord[3]-b->start[0].coord[3];
      
      test_q = true;
      for(j=0;j<b->length;j++) {
        if(shift[0]!=b2->start[j].coord[0]-b->start[j].coord[0] ||
            shift[1]!=b2->start[j].coord[1]-b->start[j].coord[1] ||
            shift[2]!=b2->start[j].coord[2]-b->start[j].coord[2] ||
            shift[3]!=b2->start[j].coord[3]-b->start[j].coord[3])
          test_q = false;
      }
      if(!test_q) continue;
    
      block_info **tmp=malloc(sizeof(block_info*)*(b->nsenders+1));
      if(b->nsenders!=0) {
        memcpy(tmp,b->senders,sizeof(block_info*)*b->nsenders);
        free(b->senders);
      }
      b->senders=tmp;
      b->senders[b->nsenders]=b2;
      b->nsenders++;
      
      /* Sender sends only to ONE receiver */
      if(b2->receiver!=NULL) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR initialize_blocks: Block %d has already got a receiver (%d)\n",b2->index,b2->receiver->index);
        print_block_info(b2);
        print_block_info(b);
	      test_q=false;
      }
      b2->receiver=b;
    }
  }
  lprintf("TEST_GEOMETRY",loglevel,"Find communication senders... OK\n");

  
  /* TEST: All buffers have a sender */
  test_q=true;
  for(b=blocks;b<blocks+nblocks;b++) {
    if(b->nbuffers>0 && b->nsenders==0) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR initialize_blocks: Block %d is a buffer but does not have a sender\n",b-blocks);
      print_block_info(b);
      test_q=false;
    }
  }
  error(!test_q,1,"test_geometry.c","All buffers have a sender... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"All buffers have a sender... OK\n");
  
  
  /* TEST: Every site belongs to a block */
  for(s=sites; s<sites+nindices; s++)
    s->test=false;

  for(b=blocks;b<blocks+nblocks;b++) {
    for(s=b->start; s<b->start+b->length; s++)
      s->test=true;
  }

  test_q=true;
  for(s=sites; s<sites+nindices; s++) {
    if(!s->test) {
      lprintf("TEST_GEOMETRY",loglevel,"This site does not belong to any valid block:\n");
      print_full_site_info(s);
    }
  }
  error(!test_q,1,"test_geometry.c","Every site belongs to a valid block... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Every site belongs to a valid block... OK\n");
  
   
  lprintf("TEST_GEOMETRY",loglevel,"List of blocks (nblocks=%d)\n",nblocks);
  for(k=0;k<nindices;k++)
  for(n=0;n<nblocks;n++) {
    if(blocks[n].start->index==k)
      print_short_block_info(&blocks[n]);
  }
}


/*
typedef struct _geometry_descriptor {
  int inner_master_pieces;
  int local_master_pieces;
  int total_spinor_master_pieces;
  int total_gauge_master_pieces;
  int *master_start, *master_end;
  int master_shift;
  int ncopies_spinor;
  int ncopies_gauge;
  int *copy_from, *copy_to, *copy_len;
  int copy_odd_shift;
  int nbuffers_spinor;
  int nbuffers_gauge;
  int *rbuf_len, *sbuf_len;
  int *rbuf_from_proc, *rbuf_start;
  int *sbuf_to_proc, *sbuf_start;
  int gsize_spinor;
  int gsize_gauge;
} geometry_descriptor;
*/


static void test_goemetry_descriptor(geometry_descriptor *gd, int parity) {
  int first, last, howmany;
  int n, i;
    int test_q;
  site_info *s;
  block_info *b;
  
  
  /* TEST: gsize_gauge */

  lprintf("TEST_GEOMETRY",loglevel,"gsize_gauge = %d\n",gd->gsize_gauge);

  if(gd->gsize_gauge!=-1) {
    first=nindices;
    last=-1;  
    for(b=blocks; b<blocks+nblocks; b++) {
      if(parity!=NOT_ASSIGNED && parity!=b->parity) continue;
      if(b->ninners+b->nlocals>=2) {
        if(b->start[0].index<first) first=b->start[0].index;
        if(b->start[b->length-1].index>last) last=b->start[b->length-1].index;
      }
    }
    howmany=last-first+1;
  
    if(gd->gsize_gauge != howmany) {
      lprintf("TEST_GEOMETRY",loglevel,"gsize_gauge should be %d\n",howmany);
      error(1,1,"test_geometry.c","Check gsize_gauge... FAILED");
    }
    lprintf("TEST_GEOMETRY",loglevel,"Check gsize_gauge... OK\n");
  }

  
  /* TEST: gsize_spinor */

  lprintf("TEST_GEOMETRY",loglevel,"gsize_spinor = %d\n",gd->gsize_spinor);

  first=nindices;
  last=-1;  
  for(b=blocks; b<blocks+nblocks; b++) {
    if(parity!=NOT_ASSIGNED && parity!=b->parity) continue;
    if(b->ninners+b->nlocals>=3) {
      if(b->start[0].index<first) first=b->start[0].index;
      if(b->start[b->length-1].index>last) last=b->start[b->length-1].index;
    }
  }
  howmany=last-first+1;

  if(gd->gsize_spinor != howmany) {
    lprintf("TEST_GEOMETRY",loglevel,"gsize_spinor should be %d\n",howmany);
    error(1,1,"test_geometry.c","Check gsize_spinor... FAILED");
  }
  lprintf("TEST_GEOMETRY",loglevel,"Check gsize_spinor... OK\n");
  
  
  /* TEST: inner_master_pieces spans all the original inner points, each once*/
    
  lprintf("TEST_GEOMETRY",loglevel,"inner_master_pieces = %d\n",gd->inner_master_pieces);
  lprintf("TEST_GEOMETRY",loglevel,"First master index = %d\n",gd->master_start[0]);
  lprintf("TEST_GEOMETRY",loglevel,"master_shift = %d\n",gd->master_shift);

  for(s=sites; s<sites+nindices; s++) s->test=0;
  
  test_q=true;
  for(n=0;n<gd->inner_master_pieces;n++)
  for(i=gd->master_start[n]; i<=gd->master_end[n]; i++) {
    sites[i].test++;
    
    if(sites[i].c_type!=ORIGINAL) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be ORIGINAL\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
    if(parity!=NOT_ASSIGNED && sites[i].parity!=parity) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has wrong parity\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
    if(sites[i].b_type[0]!=INNER || sites[i].b_type[1]!=INNER || sites[i].b_type[2]!=INNER || sites[i].b_type[3]!=INNER) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be INNER (found in master_piece=%d)\n",i,n);
      print_full_site_info(&sites[i]);
      test_q=false;
    }

    if(sites[i].test>1) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has already appeared in master pieces\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
  }
  
  for(s=sites; s<sites+nindices; s++) {
    if(s->c_type!=ORIGINAL) continue;
    if(parity!=NOT_ASSIGNED && s->parity!=parity) continue;
    if(s->b_type[0]!=INNER || s->b_type[1]!=INNER || s->b_type[2]!=INNER || s->b_type[3]!=INNER) continue;
    if(s->test==0) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d missing in master pieces\n",s->index);
      print_full_site_info(s);
      test_q=false;
    }
  }
  
  error(!test_q,1,"test_geometry.c","Check master inner pieces... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Check master inner pieces... OK\n");
  
  
  /* TEST: local_master_pieces spans all the original local points, each once*/
  
  lprintf("TEST_GEOMETRY",loglevel,"local_master_pieces = %d\n",gd->local_master_pieces);

  for(s=sites; s<sites+nindices; s++) s->test=0;
  
  test_q=true;
  for(n=0;n<gd->local_master_pieces;n++)
  for(i=gd->master_start[n]; i<=gd->master_end[n]; i++) {
    int nlocals;
    
    sites[i].test++;
    
    if(sites[i].c_type!=ORIGINAL) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be ORIGINAL\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
    
    if(parity!=NOT_ASSIGNED && sites[i].parity!=parity) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has wrong parity\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
    
    nlocals=0;
    if(sites[i].b_type[0]!=LBUFFER && sites[i].b_type[0]!=RBUFFER) nlocals++;
    if(sites[i].b_type[1]!=LBUFFER && sites[i].b_type[1]!=RBUFFER) nlocals++;
    if(sites[i].b_type[2]!=LBUFFER && sites[i].b_type[2]!=RBUFFER) nlocals++;
    if(sites[i].b_type[3]!=LBUFFER && sites[i].b_type[3]!=RBUFFER) nlocals++;
    if(nlocals!=4) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be LOCAL\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
    
    if(sites[i].test>1) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has already appeared in master pieces\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
  }
  
  for(s=sites; s<sites+nindices; s++) {
    int nlocals;
    
    if(s->c_type!=ORIGINAL) continue;
    
    if(parity!=NOT_ASSIGNED && s->parity!=parity) continue;
    
    nlocals=0;
    if(s->b_type[0]!=LBUFFER && s->b_type[0]!=RBUFFER) nlocals++;
    if(s->b_type[1]!=LBUFFER && s->b_type[1]!=RBUFFER) nlocals++;
    if(s->b_type[2]!=LBUFFER && s->b_type[2]!=RBUFFER) nlocals++;
    if(s->b_type[3]!=LBUFFER && s->b_type[3]!=RBUFFER) nlocals++;
    if(nlocals!=4) continue;
    if(s->test==0) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d missing in master pieces\n",s->index);
      print_full_site_info(s);
      test_q=false;
    }
  }
  
  error(!test_q,1,"test_geometry.c","Check master local pieces... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Check master local pieces... OK\n");
  
  
  /* TEST: total_spinor_master_pieces spans all the original local points and the 3D buffers, each once*/
  
  lprintf("TEST_GEOMETRY",loglevel,"total_spinor_master_pieces = %d\n",gd->total_spinor_master_pieces);
  
  for(s=sites; s<sites+nindices; s++) s->test=0;
  
  test_q=true;
  for(n=0;n<gd->total_spinor_master_pieces;n++)
  for(i=gd->master_start[n]; i<=gd->master_end[n]; i++) {
    int nlocals;
    
    sites[i].test++;
    
    if(sites[i].c_type!=ORIGINAL) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be ORIGINAL\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
    
    if(parity!=NOT_ASSIGNED && sites[i].parity!=parity) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has wrong parity\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
    
    nlocals=0;
    if(sites[i].b_type[0]!=LBUFFER && sites[i].b_type[0]!=RBUFFER) nlocals++;
    if(sites[i].b_type[1]!=LBUFFER && sites[i].b_type[1]!=RBUFFER) nlocals++;
    if(sites[i].b_type[2]!=LBUFFER && sites[i].b_type[2]!=RBUFFER) nlocals++;
    if(sites[i].b_type[3]!=LBUFFER && sites[i].b_type[3]!=RBUFFER) nlocals++;
    if(nlocals<3) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be LOCAL or in a 3D-buffer\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
    
    if(sites[i].test>1) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has already appeared in master pieces\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
  }
  
  for(s=sites; s<sites+nindices; s++) {
    int nlocals;
    
    if(s->c_type!=ORIGINAL) continue;
    
    if(parity!=NOT_ASSIGNED && s->parity!=parity) continue;
    
    nlocals=0;
    if(s->b_type[0]!=LBUFFER && s->b_type[0]!=RBUFFER) nlocals++;
    if(s->b_type[1]!=LBUFFER && s->b_type[1]!=RBUFFER) nlocals++;
    if(s->b_type[2]!=LBUFFER && s->b_type[2]!=RBUFFER) nlocals++;
    if(s->b_type[3]!=LBUFFER && s->b_type[3]!=RBUFFER) nlocals++;
    if(nlocals<3) continue;
    if(s->test==0) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d missing in master pieces\n",s->index);
      print_full_site_info(s);
      test_q=false;
    }
  }
  
  error(!test_q,1,"test_geometry.c","Check master local pieces + 3D buffers... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Check master local pieces + 3D buffers... OK\n");
  
  
  /* TEST: total_gauge_master_pieces spans all the original local points and the 2/3D buffers, each once*/
  
  lprintf("TEST_GEOMETRY",loglevel,"total_gauge_master_pieces = %d\n",gd->total_gauge_master_pieces);
  
  if(gd->total_gauge_master_pieces!=-1) {
    for(s=sites; s<sites+nindices; s++) s->test=0;
    
    test_q=true;
    for(n=0;n<gd->total_gauge_master_pieces;n++)
    for(i=gd->master_start[n]; i<=gd->master_end[n]; i++) {
      int nlocals;
      
      sites[i].test++;
      
      if(sites[i].c_type!=ORIGINAL) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be ORIGINAL\n",i);
        print_full_site_info(&sites[i]);
        test_q=false;
      }
      
      if(parity!=NOT_ASSIGNED && sites[i].parity!=parity) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has wrong parity\n",i);
        print_full_site_info(&sites[i]);
        test_q=false;
      }
      
      nlocals=0;
      if(sites[i].b_type[0]!=LBUFFER && sites[i].b_type[0]!=RBUFFER) nlocals++;
      if(sites[i].b_type[1]!=LBUFFER && sites[i].b_type[1]!=RBUFFER) nlocals++;
      if(sites[i].b_type[2]!=LBUFFER && sites[i].b_type[2]!=RBUFFER) nlocals++;
      if(sites[i].b_type[3]!=LBUFFER && sites[i].b_type[3]!=RBUFFER) nlocals++;
      if(nlocals<2) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be LOCAL or in a 3D-buffer\n",i);
        print_full_site_info(&sites[i]);
        test_q=false;
      }
      
      if(sites[i].test>1) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has already appeared in master pieces\n",i);
        print_full_site_info(&sites[i]);
        test_q=false;
      }
    }
    
    for(s=sites; s<sites+nindices; s++) {
      int nlocals;
      
      if(s->c_type!=ORIGINAL) continue;
      
      if(parity!=NOT_ASSIGNED && s->parity!=parity) continue;
      
      nlocals=0;
      if(s->b_type[0]!=LBUFFER && s->b_type[0]!=RBUFFER) nlocals++;
      if(s->b_type[1]!=LBUFFER && s->b_type[1]!=RBUFFER) nlocals++;
      if(s->b_type[2]!=LBUFFER && s->b_type[2]!=RBUFFER) nlocals++;
      if(s->b_type[3]!=LBUFFER && s->b_type[3]!=RBUFFER) nlocals++;
      if(nlocals<2) continue;
      if(s->test==0) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d missing in master pieces\n",s->index);
        print_full_site_info(s);
        test_q=false;
      }
    }
    
    error(!test_q,1,"test_geometry.c","Check master local pieces + 2/3D buffers... FAILED");
    lprintf("TEST_GEOMETRY",loglevel,"Check master local pieces + 2/3D buffers... OK\n");
  }


  /* TEST: all points in 3D borders are sent where they should */

  lprintf("TEST_GEOMETRY",loglevel,"nbuffers_spinor = %d\n",gd->nbuffers_spinor);
  
  for(b=blocks; b<blocks+nblocks; b++) b->test=0;

  test_q=true;
  for(n=0;n<gd->nbuffers_spinor;n++) {
    int shift[4], proc_to, mu;
    
    /* Search corresponding block */
    block_info *this=NULL;
    for(b=blocks; b<blocks+nblocks; b++) {
      if(b->start->index==gd->sbuf_start[n] && b->length==gd->sbuf_len[n] && b->receiver!=NULL) {
        this=b;
        if(this->receiver->test==1) {
          lprintf("TEST_GEOMETRY",loglevel,"sbuf[%d] start=%d len=%d\n",n,gd->sbuf_start[n],gd->sbuf_len[n]);
          lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Trying to send block %d to %d, but block %d has already received data!",this->index,this->receiver->index,this->receiver->index);
          print_block_info(this);
          print_block_info(this->receiver);
          test_q=false;
        }
        this->receiver->test=1;
        break;
      }
    }
    
    if(this==NULL) {
      lprintf("TEST_GEOMETRY",loglevel,"sbuf[%d] start=%d len=%d\n",n,gd->sbuf_start[n],gd->sbuf_len[n]);
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: No block found!");
      test_q=false;
    } else {
      /* Find the paired processor */
      shift[0]=sites[this->receiver->start->index].coord[0]-sites[this->start->index].coord[0];
      shift[1]=sites[this->receiver->start->index].coord[1]-sites[this->start->index].coord[1];
      shift[2]=sites[this->receiver->start->index].coord[2]-sites[this->start->index].coord[2];
      shift[3]=sites[this->receiver->start->index].coord[3]-sites[this->start->index].coord[3];
      
      proc_to=CID;
      mu=0;
      if(shift[mu]>0) {
        proc_to=proc_dn(proc_to,mu);
      } else if(shift[mu]<0) {
        proc_to=proc_up(proc_to,mu);
      }
      mu=1;
      if(shift[mu]>0) {
        proc_to=proc_dn(proc_to,mu);
      } else if(shift[mu]<0) {
        proc_to=proc_up(proc_to,mu);
      }
      mu=2;
      if(shift[mu]>0) {
        proc_to=proc_dn(proc_to,mu);
      } else if(shift[mu]<0) {
        proc_to=proc_up(proc_to,mu);
      }
      mu=3;
      if(shift[mu]>0) {
        proc_to=proc_dn(proc_to,mu);
      } else if(shift[mu]<0) {
        proc_to=proc_up(proc_to,mu);
      }
      
      if(gd->sbuf_to_proc[n]!=proc_to) {
        lprintf("TEST_GEOMETRY",loglevel,"sbuf[%d] start=%d len=%d\n",n,gd->sbuf_start[n],gd->sbuf_len[n]);
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Send data to processor %d, should be %d\n",gd->sbuf_to_proc[n],proc_to);
        test_q=false;
      }
    }
  }
    
  for(b=blocks; b<blocks+nblocks; b++) {
    if(b->nbuffers!=1) continue;
    if(parity!=NOT_ASSIGNED && b->senders[0]->parity!=parity) continue;
    if(b->test==0) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Block %d has not been received\n",b->index);
      print_block_info(b);
      test_q=false;
    }
  }
  
  error(!test_q,1,"test_geometry.c","Check 3D send communications... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Check 3D send communications... OK\n");




  /* TEST: all points in 3D buffers are received where they should from */

  for(b=blocks; b<blocks+nblocks; b++) b->test=0;

  test_q=true;
  for(n=0;n<gd->nbuffers_spinor;n++) {
    int shift[4], proc_from, mu;
lprintf("TEST_GEOMETRY",loglevel,"rbuf[%d] start=%d len=%d\n",n,gd->rbuf_start[n],gd->rbuf_len[n]);
    
    /* Search corresponding block */
    block_info *this=NULL;
    for(b=blocks; b<blocks+nblocks; b++) {
      if(b->start->index==gd->rbuf_start[n] && b->length==gd->rbuf_len[n] && b->nsenders!=0) {
        this=b;
        if(this->test==1) {
          lprintf("TEST_GEOMETRY",loglevel,"rbuf[%d] start=%d len=%d\n",n,gd->rbuf_start[n],gd->rbuf_len[n]);
          lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Trying to receive block %d, but it has already received data!\n",this->index);
          print_block_info(this);
          test_q=false;
        }
        this->test=1;
        break;
      }
    }
    
    if(this==NULL) {
      lprintf("TEST_GEOMETRY",loglevel,"rbuf[%d] start=%d len=%d\n",n,gd->rbuf_start[n],gd->rbuf_len[n]);
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: No block found!\n");
      print_full_site_info(sites+gd->rbuf_start[n]);
      test_q=false;
    } else {
      /* Find the paired processor */
      shift[0]=sites[this->start->index].coord[0]-sites[this->senders[0]->start->index].coord[0];
      shift[1]=sites[this->start->index].coord[1]-sites[this->senders[0]->start->index].coord[1];
      shift[2]=sites[this->start->index].coord[2]-sites[this->senders[0]->start->index].coord[2];
      shift[3]=sites[this->start->index].coord[3]-sites[this->senders[0]->start->index].coord[3];
      
      proc_from=CID;
      mu=0;
      if(shift[mu]>0) {
        proc_from=proc_up(proc_from,mu);
      } else if(shift[mu]<0) {
        proc_from=proc_dn(proc_from,mu);
      }
      mu=1;
      if(shift[mu]>0) {
        proc_from=proc_up(proc_from,mu);
      } else if(shift[mu]<0) {
        proc_from=proc_dn(proc_from,mu);
      }
      mu=2;
      if(shift[mu]>0) {
        proc_from=proc_up(proc_from,mu);
      } else if(shift[mu]<0) {
        proc_from=proc_dn(proc_from,mu);
      }
      mu=3;
      if(shift[mu]>0) {
        proc_from=proc_up(proc_from,mu);
      } else if(shift[mu]<0) {
        proc_from=proc_dn(proc_from,mu);
      }
      
      if(gd->rbuf_from_proc[n]!=proc_from) {
        lprintf("TEST_GEOMETRY",loglevel,"rbuf[%d] start=%d len=%d\n",n,gd->rbuf_start[n],gd->rbuf_len[n]);
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Receive data from processor %d, should be %d\n",gd->rbuf_from_proc[n],proc_from);
        test_q=false;
      }
    }
  }
    
  for(b=blocks; b<blocks+nblocks; b++) {
    if(b->nbuffers!=1) continue;
    if(parity!=NOT_ASSIGNED && b->parity!=parity) continue;
    if(b->test==0) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Block %d has not been received\n",b->index);
      print_block_info(b);
      test_q=false;
    }
  }
  
  error(!test_q,1,"test_geometry.c","Check 3D receive communications... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Check 3D receive communications... OK\n");
  


  /* TEST: all points in 2/3D borders are sent where they should */

  lprintf("TEST_GEOMETRY",loglevel,"nbuffers_gauge = %d\n",gd->nbuffers_gauge);
  
  if(gd->nbuffers_gauge>0){
  
    for(b=blocks; b<blocks+nblocks; b++) b->test=0;
  
    test_q=true;
    for(n=0;n<gd->nbuffers_gauge;n++) {
      int shift[4], proc_to, mu;
      
      /* Search corresponding block */
      block_info *this=NULL;
      for(b=blocks; b<blocks+nblocks; b++) {
        if(b->start->index==gd->sbuf_start[n] && b->length==gd->sbuf_len[n] && b->receiver!=NULL) {
          this=b;
          if(this->receiver->test==1) {
            lprintf("TEST_GEOMETRY",loglevel,"sbuf[%d] start=%d len=%d\n",n,gd->sbuf_start[n],gd->sbuf_len[n]);
            lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Trying to send block %d to %d, but block %d has already received data!",this->index,this->receiver->index,this->receiver->index);
            print_block_info(this);
            print_block_info(this->receiver);
            test_q=false;
          }
          this->receiver->test=1;
          break;
        }
      }
      
      if(this==NULL) {
        lprintf("TEST_GEOMETRY",loglevel,"sbuf[%d] start=%d len=%d\n",n,gd->sbuf_start[n],gd->sbuf_len[n]);
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: No block found!\n");
        test_q=false;
      } else {
        /* Find the paired processor */
        shift[0]=sites[this->receiver->start->index].coord[0]-sites[this->start->index].coord[0];
        shift[1]=sites[this->receiver->start->index].coord[1]-sites[this->start->index].coord[1];
        shift[2]=sites[this->receiver->start->index].coord[2]-sites[this->start->index].coord[2];
        shift[3]=sites[this->receiver->start->index].coord[3]-sites[this->start->index].coord[3];
        
        proc_to=CID;
        mu=0;
        if(shift[mu]>0) {
          proc_to=proc_dn(proc_to,mu);
        } else if(shift[mu]<0) {
          proc_to=proc_up(proc_to,mu);
        }
        mu=1;
        if(shift[mu]>0) {
          proc_to=proc_dn(proc_to,mu);
        } else if(shift[mu]<0) {
          proc_to=proc_up(proc_to,mu);
        }
        mu=2;
        if(shift[mu]>0) {
          proc_to=proc_dn(proc_to,mu);
        } else if(shift[mu]<0) {
          proc_to=proc_up(proc_to,mu);
        }
        mu=3;
        if(shift[mu]>0) {
          proc_to=proc_dn(proc_to,mu);
        } else if(shift[mu]<0) {
          proc_to=proc_up(proc_to,mu);
        }
        
        if(gd->sbuf_to_proc[n]!=proc_to) {
          lprintf("TEST_GEOMETRY",loglevel,"sbuf[%d] start=%d len=%d\n",n,gd->sbuf_start[n],gd->sbuf_len[n]);
          lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Send data to processor %d, should be %d\n",gd->sbuf_to_proc[n],proc_to);
          test_q=false;
        }
      }
    }
      
    for(b=blocks; b<blocks+nblocks; b++) {
      if(b->nbuffers!=1 && b->nbuffers!=2) continue;
      if(parity!=NOT_ASSIGNED && b->senders[0]->parity!=parity) continue;
      if(b->test==0) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Block %d has not been received\n",b->index);
        print_block_info(b);
        test_q=false;
      }
    }
    
    error(!test_q,1,"test_geometry.c","Check 2/3D send communications... FAILED");
    lprintf("TEST_GEOMETRY",loglevel,"Check 2/3D send communications... OK\n");

  }



  /* TEST: all points in 2/3D buffers are received where they should from */

  if(gd->nbuffers_gauge>0){
    for(b=blocks; b<blocks+nblocks; b++) b->test=0;
  
    test_q=true;
    for(n=0;n<gd->nbuffers_gauge;n++) {
      int shift[4], proc_from, mu;
      
      /* Search corresponding block */
      block_info *this=NULL;
      for(b=blocks; b<blocks+nblocks; b++) {
        if(b->start->index==gd->rbuf_start[n] && b->length==gd->rbuf_len[n] && b->nsenders!=0) {
          this=b;
          if(this->test==1) {
            lprintf("TEST_GEOMETRY",loglevel,"rbuf[%d] start=%d len=%d\n",n,gd->rbuf_start[n],gd->rbuf_len[n]);
            lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Trying to receive block %d, but it has already received data!\n",this->index);
            print_block_info(this);
            test_q=false;
          }
          this->test=1;
          break;
        }
      }
      
      if(this==NULL) {
        lprintf("TEST_GEOMETRY",loglevel,"rbuf[%d] start=%d len=%d\n",n,gd->rbuf_start[n],gd->rbuf_len[n]);
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: No block found!\n");
	print_full_site_info(sites+gd->rbuf_start[n]);
        test_q=false;
      } else {
        /* Find the paired processor */
        shift[0]=sites[this->start->index].coord[0]-sites[this->senders[0]->start->index].coord[0];
        shift[1]=sites[this->start->index].coord[1]-sites[this->senders[0]->start->index].coord[1];
        shift[2]=sites[this->start->index].coord[2]-sites[this->senders[0]->start->index].coord[2];
        shift[3]=sites[this->start->index].coord[3]-sites[this->senders[0]->start->index].coord[3];
        
        proc_from=CID;
        mu=0;
        if(shift[mu]>0) {
          proc_from=proc_up(proc_from,mu);
        } else if(shift[mu]<0) {
          proc_from=proc_dn(proc_from,mu);
        }
        mu=1;
        if(shift[mu]>0) {
          proc_from=proc_up(proc_from,mu);
        } else if(shift[mu]<0) {
          proc_from=proc_dn(proc_from,mu);
        }
        mu=2;
        if(shift[mu]>0) {
          proc_from=proc_up(proc_from,mu);
        } else if(shift[mu]<0) {
          proc_from=proc_dn(proc_from,mu);
        }
        mu=3;
        if(shift[mu]>0) {
          proc_from=proc_up(proc_from,mu);
        } else if(shift[mu]<0) {
          proc_from=proc_dn(proc_from,mu);
        }
        
        if(gd->rbuf_from_proc[n]!=proc_from) {
          lprintf("TEST_GEOMETRY",loglevel,"rbuf[%d] start=%d len=%d\n",n,gd->rbuf_start[n],gd->rbuf_len[n]);
          lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Receive data from processor %d, should be %d\n",gd->rbuf_from_proc[n],proc_from);
          test_q=false;
        }
      }
    }
      
    for(b=blocks; b<blocks+nblocks; b++) {
      if(b->nbuffers!=1 && b->nbuffers!=2) continue;
      if(parity!=NOT_ASSIGNED && b->parity!=parity) continue;
      if(b->test==0) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Block %d has not been received\n",b->index);
        print_block_info(b);
        test_q=false;
      }
    }
    
    error(!test_q,1,"test_geometry.c","Check 2/3D receive communications... FAILED");
    lprintf("TEST_GEOMETRY",loglevel,"Check 2/3D receive communications... OK\n");
  }
  
  
  
  
  /* TEST: ncopies_spinor spans all the duplicate points in 3D borders, the pairing with the original is correct, each copy is done only once */
  
  lprintf("TEST_GEOMETRY",loglevel,"ncopies_spinor = %d\n",gd->ncopies_spinor);
  lprintf("TEST_GEOMETRY",loglevel,"nbuffers_spinor = %d\n",gd->nbuffers_spinor);
  
  for(s=sites; s<sites+nindices; s++) s->test=0;
  
  test_q=true;
  for(n=0;n<gd->ncopies_spinor;n++)
  for(i=0; i<gd->copy_len[n]; i++) {
    int from, to;
    int nlocals, ninners;
    
    from=gd->copy_from[n]+i;
    to=gd->copy_to[n]+i;

    sites[to].test++;
    
    if(sites[from].c_type!=ORIGINAL || sites[to].c_type!=DUPLICATE || sites[to].original!=&sites[from]) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be ORIGINAL, and site %d should be one of its DUPLICATEs\n",from,to);
      print_full_site_info(&sites[from]);
      print_full_site_info(&sites[to]);
      test_q=false;
    }
    
    if(parity!=NOT_ASSIGNED && sites[from].parity!=parity) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has wrong parity\n",from);
      print_full_site_info(&sites[from]);
      test_q=false;
    }
    
    nlocals=0;
    if(sites[from].b_type[0]!=LBUFFER && sites[from].b_type[0]!=RBUFFER) nlocals++;
    if(sites[from].b_type[1]!=LBUFFER && sites[from].b_type[1]!=RBUFFER) nlocals++;
    if(sites[from].b_type[2]!=LBUFFER && sites[from].b_type[2]!=RBUFFER) nlocals++;
    if(sites[from].b_type[3]!=LBUFFER && sites[from].b_type[3]!=RBUFFER) nlocals++;
    ninners=0;
    if(sites[from].b_type[0]==INNER) ninners++;
    if(sites[from].b_type[1]==INNER) ninners++;
    if(sites[from].b_type[2]==INNER) ninners++;
    if(sites[from].b_type[3]==INNER) ninners++;
    if(nlocals!=4 || ninners>3) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be in a 3D-border\n",from);
      print_full_site_info(&sites[from]);
      test_q=false;
    }
    
    if(sites[to].test>1) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has already been copied\n",to);
      print_full_site_info(&sites[to]);
      test_q=false;
    }
  }
  
  for(n=0;n<gd->nbuffers_spinor;n++)
  for(i=gd->sbuf_start[n]; i<gd->sbuf_start[n]+gd->sbuf_len[n]; i++) {
    if(sites[i].c_type==DUPLICATE && sites[i].test==0) {
      lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d is DUPLICATE, but has not been copied\n",i);
      print_full_site_info(&sites[i]);
      test_q=false;
    }
  }
  
  error(!test_q,1,"test_geometry.c","Check copies in 3D borders... FAILED");
  lprintf("TEST_GEOMETRY",loglevel,"Check copies in 3D borders... OK\n");
  
  
  /* TEST: ncopies_gauge spans all the duplicate points in 2/3D borders, the pairing with the original is correct, each copy is done only once */
  
  lprintf("TEST_GEOMETRY",loglevel,"ncopies_gauge = %d\n",gd->ncopies_gauge);
  lprintf("TEST_GEOMETRY",loglevel,"nbuffers_gauge = %d\n",gd->nbuffers_gauge);
  
  if(gd->ncopies_gauge!=-1) {
    for(s=sites; s<sites+nindices; s++) s->test=0;
    
    test_q=true;
    for(n=0;n<gd->ncopies_gauge;n++)
    for(i=0; i<gd->copy_len[n]; i++) {
      int from, to;
      int nlocals, ninners;
      
      from=gd->copy_from[n]+i;
      to=gd->copy_to[n]+i;
  
      sites[to].test++;
      
      if(sites[from].c_type!=ORIGINAL || sites[to].c_type!=DUPLICATE || sites[to].original!=&sites[from]) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be ORIGINAL, and site %d should be one of its DUPLICATEs\n",from,to);
        print_full_site_info(&sites[from]);
        print_full_site_info(&sites[to]);
        test_q=false;
      }
      
      if(parity!=NOT_ASSIGNED && sites[from].parity!=parity) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has wrong parity\n",from);
        print_full_site_info(&sites[from]);
        test_q=false;
      }

      nlocals=0;
      if(sites[from].b_type[0]!=LBUFFER && sites[from].b_type[0]!=RBUFFER) nlocals++;
      if(sites[from].b_type[1]!=LBUFFER && sites[from].b_type[1]!=RBUFFER) nlocals++;
      if(sites[from].b_type[2]!=LBUFFER && sites[from].b_type[2]!=RBUFFER) nlocals++;
      if(sites[from].b_type[3]!=LBUFFER && sites[from].b_type[3]!=RBUFFER) nlocals++;
      ninners=0;
      if(sites[from].b_type[0]==INNER) ninners++;
      if(sites[from].b_type[1]==INNER) ninners++;
      if(sites[from].b_type[2]==INNER) ninners++;
      if(sites[from].b_type[3]==INNER) ninners++;
      if(nlocals!=4 || ninners>3) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d should be in a 2/3D-border\n",from);
        print_full_site_info(&sites[from]);
        test_q=false;
      }
      
      if(sites[to].test>1) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d has already been copied\n",to);
        print_full_site_info(&sites[to]);
        test_q=false;
      }
    }
    
    for(n=0;n<gd->nbuffers_gauge;n++)
    for(i=gd->sbuf_start[n]; i<gd->sbuf_start[n]+gd->sbuf_len[n]; i++) {
      if(sites[i].c_type==DUPLICATE && sites[i].test==0) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Site %d is DUPLICATE, but has not been copied\n",i);
        print_full_site_info(&sites[i]);
        test_q=false;
      }
    }

    
    error(!test_q,1,"test_geometry.c","Check copies in 2/3D borders... FAILED");
    lprintf("TEST_GEOMETRY",loglevel,"Check copies in 2/3D borders... OK\n");
  }
                                           

  /* TEST: simulate copy and communications */

#ifdef WITH_MPI
   if(gd->gsize_gauge!=-1) {
    int  mpiret;
    MPI_Request *comm_req;
    int nreq=4*gd->nbuffers_gauge;
    int* glb_coord;
    int *indices;
    
    indices=amalloc(gd->gsize_gauge*sizeof(int),ALIGN);
    glb_coord=amalloc(gd->gsize_gauge*sizeof(int)*4,ALIGN);

    for(n=0;n<gd->gsize_gauge;n++) indices[n]=sites[n].index;

    for(n=0;n<gd->gsize_gauge;n++)
      if(sites[n+gd->master_shift].c_type==ORIGINAL && sites[n+gd->master_shift].local)
        memcpy(glb_coord+4*n,sites[n+gd->master_shift].glb_coord,sizeof(int)*4);

    for(n=0;n<gd->ncopies_gauge;n++)
      memcpy(glb_coord+4*(gd->copy_to[n]-gd->master_shift), glb_coord+4*(gd->copy_from[n]-gd->master_shift), gd->copy_len[n]*sizeof(int)*4);
    
    if (gd->nbuffers_gauge>0) {
      comm_req=amalloc(4*gd->nbuffers_gauge*sizeof(MPI_Request),ALIGN);
      for (i=0; i<4*gd->nbuffers_gauge; ++i)
        comm_req[i]=MPI_REQUEST_NULL;
    } else {
      comm_req=NULL;
    }

    for (i=0; i<gd->nbuffers_gauge; ++i) {

      lprintf("TEST_GEOMETRY",loglevel,"Send/receive n=%d ... ",i);
    
      /* send ith buffer */
      mpiret=MPI_Isend(glb_coord+4*(gd->sbuf_start[i]-gd->master_shift), /* buffer */
                       (gd->sbuf_len[i])*4, /* lenght in units of doubles */
                       MPI_INT, /* basic datatype */
                       gd->sbuf_to_proc[i], /* cid of destination */
                       i, /* tag of communication */
                       cart_comm, /* use the cartesian communicator */
                       &(comm_req[4*i]) /* handle to communication request */
      );
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot start send buffer");
      }
  
      /* receive ith buffer */
      mpiret=MPI_Irecv(glb_coord+4*(gd->rbuf_start[i]-gd->master_shift), /* buffer */
                       (gd->rbuf_len[i])*4, /* lenght in units of doubles */
                       MPI_INT, /* basic datatype */
                       gd->rbuf_from_proc[i], /* cid of origin */
                       i, /* tag of communication */
                       cart_comm, /* use the cartesian communicator */
                       &(comm_req[4*i+1]) /* handle to communication request */
      );
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot start receive buffer");
      }
        
        /* send ith buffer */
        mpiret=MPI_Isend(indices+(gd->sbuf_start[i]-gd->master_shift), /* buffer */
                         gd->sbuf_len[i], /* lenght in units of doubles */
                         MPI_INT, /* basic datatype */
                         gd->sbuf_to_proc[i], /* cid of destination */
                         gd->nbuffers_gauge+i, /* tag of communication */
                         cart_comm, /* use the cartesian communicator */
                         &(comm_req[4*i+2]) /* handle to communication request */
                         );
        if (mpiret != MPI_SUCCESS) {
            char mesg[MPI_MAX_ERROR_STRING];
            int mesglen;
            MPI_Error_string(mpiret,mesg,&mesglen);
            lprintf("MPI",0,"ERROR: %s\n",mesg);
            error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot start send buffer");
        }
        
        /* receive ith buffer */
        mpiret=MPI_Irecv(indices+(gd->rbuf_start[i]-gd->master_shift), /* buffer */
                         gd->rbuf_len[i], /* lenght in units of doubles */
                         MPI_INT, /* basic datatype */
                         gd->rbuf_from_proc[i], /* cid of origin */
                         gd->nbuffers_gauge+i, /* tag of communication */
                         cart_comm, /* use the cartesian communicator */
                         &(comm_req[4*i+3]) /* handle to communication request */
                         );
        if (mpiret != MPI_SUCCESS) {
            char mesg[MPI_MAX_ERROR_STRING];
            int mesglen;
            MPI_Error_string(mpiret,mesg,&mesglen);
            lprintf("MPI",0,"ERROR: %s\n",mesg);
            error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot start receive buffer");
        }

      lprintf("TEST_GEOMETRY",loglevel,"DONE\n",i);
    
    }
  
    lprintf("TEST_GEOMETRY",loglevel,"arrivato\n");  
    if(nreq>0) {
      MPI_Status status[nreq];
      mpiret=MPI_Waitall(nreq, comm_req, status);
  
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen, k;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        for (k=0; k<nreq; ++k) {
          if (status[k].MPI_ERROR != MPI_SUCCESS) {
            MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
            lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
                k,
                status[k].MPI_SOURCE,
                status[k].MPI_TAG,
                mesg);
          }
        }
        error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot complete communications");
      }
    }
    
    test_q=true;


    for(n=0;n<gd->gsize_gauge;n++)
      if(memcmp(glb_coord+4*n,sites[n+gd->master_shift].glb_coord,sizeof(int)*4)!=0) {
        lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Global coordinates (%d,%d,%d,%d), site %d on neighbour processor\n",glb_coord[4*n], glb_coord[4*n+1], glb_coord[4*n+2], glb_coord[4*n+3],indices[n]);
        print_full_site_info(&sites[n+gd->master_shift]);
        print_full_site_info(&sites[indices[n]]);
        test_q=false;
      }
        

    error(!test_q,1,"test_geometry.c","Global coordinates sent and received correctly (gauge)... FAILED");
    lprintf("TEST_GEOMETRY",loglevel,"Global coordinates sent and received correctly (gauge)... OK\n");
  
    afree(comm_req);
    afree(glb_coord);
  }

  
    if(gd->gsize_spinor!=-1) {
        int  mpiret;
        MPI_Request *comm_req;
        int nreq=4*gd->nbuffers_spinor;
        int* glb_coord;
        int *indices;
        
        indices=amalloc(gd->gsize_spinor*sizeof(int),ALIGN);
        glb_coord=amalloc(gd->gsize_spinor*sizeof(int)*4,ALIGN);
        
        for(n=0;n<gd->gsize_spinor;n++) indices[n]=sites[n].index;
        
        for(n=0;n<gd->gsize_spinor;n++)
            if(sites[n+gd->master_shift].c_type==ORIGINAL && sites[n+gd->master_shift].local)
                memcpy(glb_coord+4*n,sites[n+gd->master_shift].glb_coord,sizeof(int)*4);
        
        for(n=0;n<gd->ncopies_spinor;n++)
            memcpy(glb_coord+4*(gd->copy_to[n]-gd->master_shift), glb_coord+4*(gd->copy_from[n]-gd->master_shift), gd->copy_len[n]*sizeof(int)*4);
        
        if (gd->nbuffers_spinor>0) {
            comm_req=amalloc(4*gd->nbuffers_spinor*sizeof(MPI_Request),ALIGN);
            for (i=0; i<4*gd->nbuffers_spinor; ++i)
                comm_req[i]=MPI_REQUEST_NULL;
        } else {
            comm_req=NULL;
        }
        
        for (i=0; i<gd->nbuffers_spinor; ++i) {
            
            lprintf("TEST_GEOMETRY",loglevel,"Send/receive n=%d ... ",i);
            
            /* send ith buffer */
            mpiret=MPI_Isend(glb_coord+4*(gd->sbuf_start[i]-gd->master_shift), /* buffer */
                             (gd->sbuf_len[i])*4, /* lenght in units of doubles */
                             MPI_INT, /* basic datatype */
                             gd->sbuf_to_proc[i], /* cid of destination */
                             i, /* tag of communication */
                             cart_comm, /* use the cartesian communicator */
                             &(comm_req[4*i]) /* handle to communication request */
                             );
            if (mpiret != MPI_SUCCESS) {
                char mesg[MPI_MAX_ERROR_STRING];
                int mesglen;
                MPI_Error_string(mpiret,mesg,&mesglen);
                lprintf("MPI",0,"ERROR: %s\n",mesg);
                error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot start send buffer");
            }
            
            /* receive ith buffer */
            mpiret=MPI_Irecv(glb_coord+4*(gd->rbuf_start[i]-gd->master_shift), /* buffer */
                             (gd->rbuf_len[i])*4, /* lenght in units of doubles */
                             MPI_INT, /* basic datatype */
                             gd->rbuf_from_proc[i], /* cid of origin */
                             i, /* tag of communication */
                             cart_comm, /* use the cartesian communicator */
                             &(comm_req[4*i+1]) /* handle to communication request */
                             );
            if (mpiret != MPI_SUCCESS) {
                char mesg[MPI_MAX_ERROR_STRING];
                int mesglen;
                MPI_Error_string(mpiret,mesg,&mesglen);
                lprintf("MPI",0,"ERROR: %s\n",mesg);
                error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot start receive buffer");
            }
            
            /* send ith buffer */
            mpiret=MPI_Isend(indices+(gd->sbuf_start[i]-gd->master_shift), /* buffer */
                             gd->sbuf_len[i], /* lenght in units of doubles */
                             MPI_INT, /* basic datatype */
                             gd->sbuf_to_proc[i], /* cid of destination */
                             gd->nbuffers_spinor+i, /* tag of communication */
                             cart_comm, /* use the cartesian communicator */
                             &(comm_req[4*i+2]) /* handle to communication request */
                             );
            if (mpiret != MPI_SUCCESS) {
                char mesg[MPI_MAX_ERROR_STRING];
                int mesglen;
                MPI_Error_string(mpiret,mesg,&mesglen);
                lprintf("MPI",0,"ERROR: %s\n",mesg);
                error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot start send buffer");
            }
            
            /* receive ith buffer */
            mpiret=MPI_Irecv(indices+(gd->rbuf_start[i]-gd->master_shift), /* buffer */
                             gd->rbuf_len[i], /* lenght in units of doubles */
                             MPI_INT, /* basic datatype */
                             gd->rbuf_from_proc[i], /* cid of origin */
                             gd->nbuffers_spinor+i, /* tag of communication */
                             cart_comm, /* use the cartesian communicator */
                             &(comm_req[4*i+3]) /* handle to communication request */
                             );
            if (mpiret != MPI_SUCCESS) {
                char mesg[MPI_MAX_ERROR_STRING];
                int mesglen;
                MPI_Error_string(mpiret,mesg,&mesglen);
                lprintf("MPI",0,"ERROR: %s\n",mesg);
                error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot start receive buffer");
            }
            
            lprintf("TEST_GEOMETRY",loglevel,"DONE\n",i);
            
        }
        
        
        if(nreq>0) {
            MPI_Status status[nreq];
            mpiret=MPI_Waitall(nreq, comm_req, status);
            
            if (mpiret != MPI_SUCCESS) {
                char mesg[MPI_MAX_ERROR_STRING];
                int mesglen, k;
                MPI_Error_string(mpiret,mesg,&mesglen);
                lprintf("MPI",0,"ERROR: %s\n",mesg);
                for (k=0; k<nreq; ++k) {
                    if (status[k].MPI_ERROR != MPI_SUCCESS) {
                        MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
                        lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
                                k,
                                status[k].MPI_SOURCE,
                                status[k].MPI_TAG,
                                mesg);
                    }
                }
                error(1,1,"test_goemetry_descriptor [test_geometry_mpi.c]","Cannot complete communications");
            }
        }
        
        test_q=true;
        
        for(n=0;n<gd->gsize_spinor;n++)
            if(memcmp(glb_coord+4*n,sites[n+gd->master_shift].glb_coord,sizeof(int)*4)!=0) {
                lprintf("TEST_GEOMETRY",loglevel,"ERROR test_goemetry_descriptor: Global coordinates (%d,%d,%d,%d), site %d on neighbour processor\n",glb_coord[4*n], glb_coord[4*n+1], glb_coord[4*n+2], glb_coord[4*n+3],indices[n]);
                print_full_site_info(&sites[n+gd->master_shift]);
                print_full_site_info(&sites[indices[n]]);
                test_q=false;
            }
        
        
        error(!test_q,1,"test_geometry.c","Global coordinates sent and received correctly (spinor)... FAILED");
        lprintf("TEST_GEOMETRY",loglevel,"Global coordinates sent and received correctly (spinor)... OK\n");
        
        afree(comm_req);
        afree(glb_coord);
    }

  
#endif
  
}


void test_geometry_mpi_eo() {

  initialize_sites();

  initialize_blocks();
  
  lprintf("TEST_GEOMETRY",loglevel,"********************* glattice\n");
  test_goemetry_descriptor(&glattice,NOT_ASSIGNED);
  
  lprintf("TEST_GEOMETRY",loglevel,"********************* glat_even\n");
  test_goemetry_descriptor(&glat_even,EVEN);
  
  lprintf("TEST_GEOMETRY",loglevel,"********************* glat_odd\n");
  test_goemetry_descriptor(&glat_odd,ODD);
  
  finalize_test();
}

