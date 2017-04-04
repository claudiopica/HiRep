#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "update.h"
#include "utils.h"
#include "communications.h"
#include "observables.h"
#include <math.h>
#include <string.h>

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional
#endif

#ifdef BC_XYZ_TWISTED
#error This code does not work with the twisted BCs
#endif

#ifdef BC_T_OPEN
#warning This code does not work with the open BCs
#endif

//#error "Wilson loop must be fixed, The zerocoord global location must be added"

#define _WL_4VOL_INDEX(t,x,y,z) ((t)+(x)*T+(y)*T*X+(z)*T*X*Y)

static suNg_field* ws_gtf[2]; /* gtf = gauge transformation field */
static suNg* buf_gtf[3];
static suNg* Polyakov;
static suNg_field* HYP;
static int rot12=(1==0);
static int rot13=(1==0);
static int rot23=(1==0);

struct {
  int c[3];
  int* path;
  int length;
  int** perm;
  int nperms;
  int nsteps;
} WL_path[256];
int WL_npaths=0;
int WL_max_nsteps=0;



static int WL_init=(1==0);
void WL_initialize() {
  if(WL_init) return;
  
  error(NP_X!=1 || NP_Y!=1 || NP_Z !=1,1,"WL_initialize [wilsonloops.c]","The Wilson loops code does not work with spatial paralelization!");

  HYP=alloc_gfield(&glattice);
  ws_gtf[0]=alloc_gtransf(&glattice);
  ws_gtf[1]=alloc_gtransf(&glattice);
  buf_gtf[0]=amalloc(sizeof(suNg)*T*X*Y*Z,ALIGN);
  buf_gtf[1]=amalloc(sizeof(suNg)*T*X*Y*Z,ALIGN);
  buf_gtf[2]=amalloc(sizeof(suNg)*T*X*Y*Z,ALIGN);
  Polyakov=amalloc(sizeof(suNg)*X*Y*Z,ALIGN);

  #if (defined(BC_X_ANTIPERIODIC) && defined(BC_Y_ANTIPERIODIC)) || (defined(BC_X_PERIODIC) && defined(BC_Y_PERIODIC))
    if(GLB_X==GLB_Y) rot12=(1==1);
  #endif

  #if (defined(BC_X_ANTIPERIODIC) && defined(BC_Z_ANTIPERIODIC)) || (defined(BC_X_PERIODIC) && defined(BC_Z_PERIODIC))
    if(GLB_X==GLB_Z) rot13=(1==1);
  #endif

  #if (defined(BC_Z_ANTIPERIODIC) && defined(BC_Y_ANTIPERIODIC)) || (defined(BC_Z_PERIODIC) && defined(BC_Y_PERIODIC))
    if(GLB_Z==GLB_Y) rot23=(1==1);
  #endif
  
  WL_init=(1==1);
}




void WL_free() {
  if(!WL_init) return;
  
  free_gfield(HYP);
  free_gtransf(ws_gtf[0]);
  free_gtransf(ws_gtf[1]);
  afree(buf_gtf[0]);
  afree(buf_gtf[1]);
  afree(buf_gtf[2]);
  afree(Polyakov);
  
  for(int i=0; i<WL_npaths; i++) {
    afree(WL_path[i].path);
    afree(WL_path[i].perm);
  }
}



/* Variation on hep-lat/0005018 */
static void WL_3Dpath(int c[3], int* path) {
  int p[3], q[3];
  double dist[3];
  
  c[0]=(c[0]>=0)?c[0]:-c[0];
  c[1]=(c[1]>=0)?c[1]:-c[1];
  c[2]=(c[2]>=0)?c[2]:-c[2];

  p[0]=p[1]=p[2]=0;
  for(int i=0; i<c[0]+c[1]+c[2]; i++) {
    for(int k=0; k<3; k++) {
      q[0]=p[0];q[1]=p[1];q[2]=p[2];
      q[k]++;
      #define _PROD(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
      #define _SQUARE(a) ((a)*(a))
      dist[k]=_PROD(q,q)*_PROD(c,c) - _SQUARE(_PROD(q,c));
      #undef _PROD
      #undef _SQUARE
    }
    
    if(dist[0]<=dist[1] && dist[0]<=dist[2]) path[i]=0;
    else if(dist[1]<dist[0] && dist[1]<=dist[2]) path[i]=1;
    else path[i]=2;
    
    p[path[i]]++;
  }

   error(p[0]!=c[0] || p[1]!=c[1] || p[2]!=c[2],1,"WL_make_path [wilsonloops.c]","Wrong end point");
 
}



static int abc_perm[6][3]={{0,1,2},{1,2,0},{2,0,1},{2,1,0},{1,0,2},{0,2,1}};
static int aab_perm[3][3]={{0,1,2},{0,2,1},{2,1,0}};
static int abb_perm[3][3]={{0,1,2},{1,0,2},{2,1,0}};
static int aaa_perm[1][3]={{0,1,2}};
static int aaX_perm[1][3]={{0,1,2}};
static int abX_perm[2][3]={{0,1,2},{1,0,2}};
static int aXa_perm[1][3]={{0,1,2}};
static int aXb_perm[2][3]={{0,1,2},{2,1,0}};
static int Xaa_perm[1][3]={{0,1,2}};
static int Xab_perm[2][3]={{0,1,2},{0,2,1}};

void WL_load_path(int c[3], int nsteps) {
  c[0]=(c[0]>=0)?c[0]:-c[0];
  c[1]=(c[1]>=0)?c[1]:-c[1];
  c[2]=(c[2]>=0)?c[2]:-c[2];
  
  #define _SWAP(a,b) {int tmp=a;a=b;b=tmp;}
  if(rot12) {
    if(c[0]<c[1]) _SWAP(c[0],c[1]);
  }
  if(rot23) {
    if(c[1]<c[2]) _SWAP(c[1],c[2]);
  }
  if(rot12) {
    if(c[0]<c[1]) _SWAP(c[0],c[1]);
  }
  if(rot13) {
    if(c[0]<c[2]) _SWAP(c[0],c[2]);
  }
  #undef _SWAP
  
  for(int i=0; i<WL_npaths; i++) {
    if(c[0]==WL_path[i].c[0] && c[1]==WL_path[i].c[1] && c[2]==WL_path[i].c[2]) return;
  }
  
  WL_path[WL_npaths].c[0]=c[0];
  WL_path[WL_npaths].c[1]=c[1];
  WL_path[WL_npaths].c[2]=c[2];
  WL_path[WL_npaths].length=c[0]+c[1]+c[2];
  WL_path[WL_npaths].path=amalloc(sizeof(int)*(c[0]+c[1]+c[2]),ALIGN);
  WL_3Dpath(c,WL_path[WL_npaths].path);
  
  int max_nsteps=nsteps;
  if(c[0]!=0) max_nsteps = (max_nsteps<GLB_Y/c[0])?max_nsteps:GLB_X/c[0];
  if(c[1]!=0) max_nsteps = (max_nsteps<GLB_Y/c[1])?max_nsteps:GLB_Y/c[1];
  if(c[2]!=0) max_nsteps = (max_nsteps<GLB_Z/c[2])?max_nsteps:GLB_Z/c[2];
  
  if(nsteps>max_nsteps) {
    lprintf("WILSON LOOPS",0,"WARNING!!! nsteps reduced from %d to %d for c=(%d,%d,%d)\n",nsteps,max_nsteps,c[0],c[1],c[2]);
    nsteps=max_nsteps;
  }
  error(nsteps==0,1,"WL_load_path [wilsonloops.c]","nsteps == 0!");
  
  WL_path[WL_npaths].nsteps=nsteps;
  if(nsteps>WL_max_nsteps) WL_max_nsteps=nsteps;
  
  if(rot12 && rot23) {
    if(c[0]==c[1] && c[1]==c[2]) {
      WL_path[WL_npaths].nperms=1;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=aaa_perm[w];
    } else if(c[0]==c[1]) {
      WL_path[WL_npaths].nperms=3;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=aab_perm[w];
    } else if(c[1]==c[2]) {
      WL_path[WL_npaths].nperms=3;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=abb_perm[w];
    } else {
      WL_path[WL_npaths].nperms=6;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=abc_perm[w];
    }
  } else if(rot12) {
    if(c[0]==c[1]) {
      WL_path[WL_npaths].nperms=1;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=aaX_perm[w];
    } else {
      WL_path[WL_npaths].nperms=2;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=abX_perm[w];
    }
  } else if(rot23) {
    if(c[1]==c[2]) {
      WL_path[WL_npaths].nperms=1;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=Xaa_perm[w];
    } else {
      WL_path[WL_npaths].nperms=2;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=Xab_perm[w];
    }
  } else if(rot13) {
    if(c[0]==c[2]) {
      WL_path[WL_npaths].nperms=1;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=aXa_perm[w];
    } else {
      WL_path[WL_npaths].nperms=2;
      WL_path[WL_npaths].perm=amalloc(sizeof(int*)*WL_path[WL_npaths].nperms,ALIGN);
      for(int w=0;w<WL_path[WL_npaths].nperms;w++) WL_path[WL_npaths].perm[w]=aXb_perm[w];
    }
  }
  
  lprintf("ARA WILSON LOOPS",0,"c=( %d , %d , %d ) path added (nsteps = %d ; length = %d ; nperms = %d )\n",c[0],c[1],c[2],nsteps,WL_path[WL_npaths].length,WL_path[WL_npaths].nperms);   
  
  WL_npaths++;
}



void WL_Hamiltonian_gauge(suNg_field* out, suNg_field* in) {
  int i,j;
  int x,y,z;
#ifdef WITH_MPI  
  int mpiret;
#endif /* WITH_MPI */

  /* LOC(t)=zerocoord[0]+t */
  /* ws_gtf[0](t) = U_0(LOC(0)) ... U_0(LOC(t)) */
  for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
    i=ipt(0,x,y,z);
    *_FIELD_AT(ws_gtf[0],i)=*_4FIELD_AT(in,i,0);
    for(int t=1;t<T;t++) {
      i=ipt(t,x,y,z);
      j=ipt(t-1,x,y,z);
      _suNg_times_suNg(*_FIELD_AT(ws_gtf[0],i),*_FIELD_AT(ws_gtf[0],j),*_4FIELD_AT(in,i,0));
    }
  }
  

  /* buf_gtf[1] = U_0(0) ... U_0(LOC(T-1)-T) */
#ifdef WITH_MPI  
  if(COORD[0]!=0) {
    MPI_Status status;
    mpiret=MPI_Recv(buf_gtf[1], /* buffer */
        (X*Y*Z)*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        proc_dn(CID,0), /* cid of origin */
        COORD[0]-1, /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &status);
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      if (status.MPI_ERROR != MPI_SUCCESS) {
        MPI_Error_string(status.MPI_ERROR,mesg,&mesglen);
        lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
            0, 
            status.MPI_SOURCE, 
            status.MPI_TAG, 
            mesg);
      }
      error(1,1,"WL_Hamiltonian_gauge [wilsonloops.c]","Cannot receive buf_gtf[1]");
    }
#endif /* NDEBUG */
  }
#endif /* WITH_MPI */


  /* buf_gtf[0] = U_0(0) ... U_0(LOC(T-1)) */
  if(COORD[0]==0) {
    for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
      i=ipt(T-1,x,y,z);
      buf_gtf[0][_WL_3VOL_INDEX(x,y,z)]=*_FIELD_AT(ws_gtf[0],i);
    }
  } else {
    for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
      i=ipt(T-1,x,y,z);
      _suNg_times_suNg(buf_gtf[0][_WL_3VOL_INDEX(x,y,z)],buf_gtf[1][_WL_3VOL_INDEX(x,y,z)],*_FIELD_AT(ws_gtf[0],i));
    }    
  }
  
  
#ifdef WITH_MPI  
  if(COORD[0]!=NP_T-1) {
    mpiret=MPI_Send(buf_gtf[0], /* buffer */
        (X*Y*Z)*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        proc_up(CID,0), /* cid of destination */
        COORD[0], /* tag of communication */
        cart_comm /* use the cartesian communicator */
        );
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"WL_Hamiltonian_gauge [wilsonloops.c]","Cannot send buf_gtf[0]");
    }
#endif /* NDEBUG */
  }
#endif /* WITH_MPI */


  /* ws_gtf[1](t) = U_0(0) ... U_0(LOC(t-1)) */
  if(COORD[0]==0) {
    for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
      i=ipt(0,x,y,z);
      _suNg_unit(*_FIELD_AT(ws_gtf[1],i));
      for(int t=1;t<T;t++) {
        i=ipt(t,x,y,z);
        j=ipt(t-1,x,y,z);
        *_FIELD_AT(ws_gtf[1],i)=*_FIELD_AT(ws_gtf[0],j);
      }
    }
  } else {
    for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
      i=ipt(0,x,y,z);
      *_FIELD_AT(ws_gtf[1],i)=buf_gtf[1][_WL_3VOL_INDEX(x,y,z)];
      for(int t=1;t<T;t++) {
        i=ipt(t,x,y,z);
        j=ipt(t-1,x,y,z);
        _suNg_times_suNg(*_FIELD_AT(ws_gtf[1],i),buf_gtf[1][_WL_3VOL_INDEX(x,y,z)],*_FIELD_AT(ws_gtf[0],j));
      }
    }
  }
  start_gt_sendrecv(ws_gtf[1]);
  
  _PIECE_FOR(&glattice,ixp) {
    suNg tmp;
    if(ixp==glattice.inner_master_pieces) {
      _OMP_PRAGMA( master )
      complete_gt_sendrecv(ws_gtf[1]);
      _OMP_PRAGMA( barrier )
    }
    _SITE_FOR(&glattice,ixp,ix) {
      for(int mu=0;mu<4;mu++) {
        _suNg_times_suNg(tmp,*_FIELD_AT(ws_gtf[1],ix),*_4FIELD_AT(in,ix,mu));
        _suNg_times_suNg_dagger(*_4FIELD_AT(out,ix,mu),tmp,*_FIELD_AT(ws_gtf[1],iup(ix,mu)));
      }
    } /* SITE_FOR */
  } /* PIECE FOR */

  start_gf_sendrecv(out);
  complete_gf_sendrecv(out);
}



void WL_broadcast_polyakov(suNg* poly, suNg_field* gf) {
  int x,y,z;
  if(COORD[0]==NP_T-1) {
    for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
      int i=ipt(T-1,x,y,z);
      poly[_WL_3VOL_INDEX(x,y,z)]=*_4FIELD_AT(gf,i,0);
    }
  }

#ifdef WITH_MPI
  int mpiret;
  if(COORD[0]==NP_T-1) {
    MPI_Request comm_req[NP_T-1];
    int destCID=CID;
    for(int t=0; t<NP_T-1; t++) {
      destCID=proc_up(destCID,0);
      mpiret=MPI_Isend(poly, /* buffer */
          (X*Y*Z)*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
          MPI_DOUBLE, /* basic datatype */
          destCID, /* cid of destination */
          t, /* tag of communication */
          cart_comm, /* use the cartesian communicator */
          &(comm_req[t]) /* handle to communication request */
          );
#ifndef NDEBUG
      if (mpiret != MPI_SUCCESS) {
        char mesg[MPI_MAX_ERROR_STRING];
        int mesglen;
        MPI_Error_string(mpiret,mesg,&mesglen);
        lprintf("MPI",0,"ERROR: %s\n",mesg);
        error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot start send polyakov");
      }
#endif /* NDEBUG */
    }
    
    MPI_Status status[NP_T-1];
    mpiret=MPI_Waitall(NP_T-1, comm_req, status);
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen, k;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      for (k=0; k<NP_T-1; ++k) {
        if (status[k].MPI_ERROR != MPI_SUCCESS) {
          MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
          lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
              k, 
              status[k].MPI_SOURCE, 
              status[k].MPI_TAG, 
              mesg);
        }
      }
      error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot complete communications");
    }
#endif /* NDEBUG */

  } else {

    int sCOORD[4], sCID;
    sCOORD[0]=NP_T-1;sCOORD[1]=COORD[1];sCOORD[2]=COORD[2];sCOORD[3]=COORD[3];
    mpiret=MPI_Cart_rank(cart_comm, sCOORD, &sCID);
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot retrieve source CID");
    }
#endif /* NDEBUG */
    MPI_Status status;
    mpiret=MPI_Recv(poly, /* buffer */
        (X*Y*Z)*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
        MPI_DOUBLE, /* basic datatype */
        sCID, /* cid of destination */
        COORD[0], /* tag of communication */
        cart_comm, /* use the cartesian communicator */
        &status);
#ifndef NDEBUG
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      if (status.MPI_ERROR != MPI_SUCCESS) {
        MPI_Error_string(status.MPI_ERROR,mesg,&mesglen);
        lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
            0, 
            status.MPI_SOURCE, 
            status.MPI_TAG, 
            mesg);
      }
      error(1,1,"WL_broadcast_polyakov [wilsonloops.c]","Cannot receive polyakov");
    }
#endif /* NDEBUG */
  }
#endif /* WITH_MPI */
}



static void WL_parallel_transport(suNg* ret, const suNg_field* gf, int x, const int* path, const int length, const int perm[3], const int sign[3]) {
  suNg tmp;
  int inv[3];
  
  inv[perm[0]]=0;
  inv[perm[1]]=1;
  inv[perm[2]]=2;  
  
  _suNg_unit(*ret);
  for(int i=0; i<length; i++) {
    if(sign[path[i]]>0) {
      _suNg_times_suNg(tmp,*ret,*_4FIELD_AT(gf,x,inv[path[i]]+1));
      *ret=tmp;
      x=iup(x,inv[path[i]]+1);
    } else {
      x=idn(x,inv[path[i]]+1);
      _suNg_times_suNg_dagger(tmp,*ret,*_4FIELD_AT(gf,x,inv[path[i]]+1));
      *ret=tmp;
    }
  }
}



void WL_correlators(double** ret, const suNg_field* gf, const suNg* poly, const int nsteps, const int* path, const int length, const int perm[3], int sign[3]) {
  suNg tmp[2];
  int c[3];
  
  c[0]=c[1]=c[2]=0;
  for(int i=0; i<length; i++) c[path[i]]++;
  
  sign[0]=(sign[0]>0)?1:-1;
  sign[1]=(sign[1]>0)?1:-1;
  sign[2]=(sign[2]>0)?1:-1;
  
  for(int s=0; s<nsteps; s++) for(int t=0;t<GLB_T;t++) ret[s][t]=0.;
  
/*  print=0;*/
  
  for(int s=0; s<nsteps; s++) {
    if(s==0) {
      for(int t=0;t<T;t++) for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
        WL_parallel_transport(_FIELD_AT(ws_gtf[0],ipt(t,x,y,z)), gf, ipt(t,x,y,z), path, length, perm, sign);
        buf_gtf[0][_WL_4VOL_INDEX(t,x,y,z)]=*_FIELD_AT(ws_gtf[0],ipt(t,x,y,z));
      }
    } else {
      for(int t=0;t<T;t++) for(int x0=0;x0<X;x0++) for(int y0=0;y0<Y;y0++) for(int z0=0;z0<Z;z0++) {
        int x1=(x0+sign[perm[0]]*c[perm[0]]*s+GLB_X)%GLB_X;
        int y1=(y0+sign[perm[1]]*c[perm[1]]*s+GLB_Y)%GLB_Y;
        int z1=(z0+sign[perm[2]]*c[perm[2]]*s+GLB_Z)%GLB_Z;
        _suNg_times_suNg(tmp[0],buf_gtf[0][_WL_4VOL_INDEX(t,x0,y0,z0)],*_FIELD_AT(ws_gtf[0],ipt(t,x1,y1,z1)));
        buf_gtf[0][_WL_4VOL_INDEX(t,x0,y0,z0)]=tmp[0];
      }
    }
    
#ifdef WITH_MPI
    MPI_Request comm_req[2];
    MPI_Status status[2];
    int mpiret;
    int destCID=CID;
    int sendCID=CID;
#endif /* WITH_MPI */
    memcpy(buf_gtf[1],buf_gtf[0],(T*X*Y*Z)*sizeof(suNg));
    
    for(int DT=1; DT<=NP_T; DT++) {
#ifdef WITH_MPI
      /* start communication for DT */
      if(DT<NP_T) {
        destCID=proc_dn(destCID,0);
        mpiret=MPI_Isend(buf_gtf[0], /* buffer */
            (T*X*Y*Z)*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
            MPI_DOUBLE, /* basic datatype */
            destCID, /* cid of destination */
            DT, /* tag of communication */
            cart_comm, /* use the cartesian communicator */
            &(comm_req[0]) /* handle to communication request */
            );
#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          error(1,1,"WL_correlators [wilsonloops.c]","Cannot start send buffer");
        }
#endif /* NDEBUG */
  
        sendCID=proc_up(sendCID,0);
        mpiret=MPI_Irecv(buf_gtf[2], /* buffer */
            (T*X*Y*Z)*sizeof(suNg)/sizeof(double), /* lenght in units of doubles */
            MPI_DOUBLE, /* basic datatype */
            sendCID, /* cid of origin */
            DT, /* tag of communication */
            cart_comm, /* use the cartesian communicator */
            &(comm_req[1]) /* handle to communication request */
            );
#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          error(1,1,"WL_correlators [wilsonloops.c]","Cannot start receive buffer");
        }
#endif /* NDEBUG */
      }
#endif /* WITH_MPI */

      /* computation of the Wilson loops for DT-1 */
      double dtmp;
      for(int t0=0;t0<T;t0++) for(int t1=0;t1<T;t1++) {
        if(((COORD[0]+DT-1)*T+t1)%GLB_T>=zerocoord[0]+t0) {
          for(int x0=0;x0<X;x0++) for(int y0=0;y0<Y;y0++) for(int z0=0;z0<Z;z0++) {
            _suNg_times_suNg_dagger(tmp[0],
              buf_gtf[0][_WL_4VOL_INDEX(t0,x0,y0,z0)],
              buf_gtf[1][_WL_4VOL_INDEX(t1,x0,y0,z0)]);
            _suNg_trace_re(dtmp,tmp[0]);
            ret[s][((DT-1)*T+t1-t0+GLB_T)%GLB_T]+=dtmp/GLB_VOLUME;
          }
        } else {
          for(int x0=0;x0<X;x0++) for(int y0=0;y0<Y;y0++) for(int z0=0;z0<Z;z0++) {
            int x1=(x0+sign[perm[0]]*c[perm[0]]*(s+1)+GLB_X)%GLB_X;
            int y1=(y0+sign[perm[1]]*c[perm[1]]*(s+1)+GLB_Y)%GLB_Y;
            int z1=(z0+sign[perm[2]]*c[perm[2]]*(s+1)+GLB_Z)%GLB_Z;
            _suNg_times_suNg(tmp[0],buf_gtf[0][_WL_4VOL_INDEX(t0,x0,y0,z0)],poly[_WL_3VOL_INDEX(x1,y1,z1)]);
            _suNg_times_suNg_dagger(tmp[1],tmp[0],buf_gtf[1][_WL_4VOL_INDEX(t1,x0,y0,z0)]);
            _suNg_times_suNg_dagger(tmp[0],tmp[1],poly[_WL_3VOL_INDEX(x0,y0,z0)]);
            _suNg_trace_re(dtmp,tmp[0]);
            ret[s][((DT-1)*T+t1-t0+GLB_T)%GLB_T]+=dtmp/GLB_VOLUME;
          }
        }
      }
    
#ifdef WITH_MPI
      /* wait for communication for DT */
      if(DT<NP_T) {
        mpiret=MPI_Waitall(2, comm_req, status);

#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen, k;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          for (k=0; k<2; ++k) {
            if (status[k].MPI_ERROR != MPI_SUCCESS) {
              MPI_Error_string(status[k].MPI_ERROR,mesg,&mesglen);
              lprintf("MPI",0,"Req [%d] Source [%d] Tag [%] ERROR: %s\n",
                  k, 
                  status[k].MPI_SOURCE, 
                  status[k].MPI_TAG, 
                  mesg);
            }
          }
          error(1,1,"WL_correlators [wilsonloops.c]","Cannot complete communications");
        }
#endif /* NDEBUG */

        memcpy(buf_gtf[1],buf_gtf[2],(T*X*Y*Z)*sizeof(suNg));
      }
#endif /* WITH_MPI */
    }
    
    global_sum(ret[s],GLB_T);
  }
}



void WL_wilsonloops(double HYP_weight[3]) {
  error(WL_npaths==0,1,"WL_wilsonloops [wilsonloops.c]","No path has been loaded");
  
  HYP_smearing(HYP,u_gauge,HYP_weight);

  WL_Hamiltonian_gauge(HYP,HYP);

  WL_broadcast_polyakov(Polyakov,HYP);

/*
  WL_Hamiltonian_gauge(u_gauge,u_gauge);
  
  WL_broadcast_polyakov(Polyakov,u_gauge);

  HYP_smearing(HYP,u_gauge,HYP_weight);
*/
  
  double** WL;
  WL=amalloc(sizeof(double*)*WL_max_nsteps,ALIGN);
  WL[0]=amalloc(sizeof(double)*WL_max_nsteps*GLB_T,ALIGN);
  for(int s=0; s<WL_max_nsteps; s++)
    WL[s]=WL[0]+s*GLB_T;
  
  double** tmp;
  tmp=amalloc(sizeof(double*)*WL_max_nsteps,ALIGN);
  tmp[0]=amalloc(sizeof(double)*WL_max_nsteps*GLB_T,ALIGN);
  for(int s=0; s<WL_max_nsteps; s++)
    tmp[s]=tmp[0]+s*GLB_T;
  
  for(int n=0; n<WL_npaths; n++) {
    for(int s=0;s<WL_path[n].nsteps;s++) for(int t=0;t<GLB_T;t++) WL[s][t]=0.;
    
    int counter=0;
    for(int p=0;p<WL_path[n].nperms;p++) {
      for(int w=0;w<4;w++) {
        int sign[3];
        sign[0]=(w%2==0)?1:-1;
        sign[1]=((w/2)%2==0)?1:-1;
        sign[2]=1;
        if(WL_path[n].c[0]==0 && sign[0]==-1) continue;
        if(WL_path[n].c[1]==0 && sign[1]==-1) continue;
        if(WL_path[n].c[2]==0 && sign[2]==-1) continue;
      
        WL_correlators(tmp,HYP,Polyakov,WL_path[n].nsteps,WL_path[n].path,WL_path[n].length,WL_path[n].perm[p],sign);
        for(int s=0;s<WL_path[n].nsteps;s++) for(int t=0;t<GLB_T;t++) WL[s][t]+=tmp[s][t];
        counter++;
      }
    }

    for(int s=0;s<WL_path[n].nsteps;s++) {
      for(int t=0;t<GLB_T;t++) WL[s][t]/=counter;

      for(int t=0;t<GLB_T;t++) {
        lprintf("WILSON LOOPS",0,"(T,dx,dy,dz,R,WL) = %d %d %d %d %.8e %.8e\n",
                t,(s+1)*WL_path[n].c[0],(s+1)*WL_path[n].c[1],(s+1)*WL_path[n].c[2],
                (s+1)*sqrt(WL_path[n].c[0]*WL_path[n].c[0]+WL_path[n].c[1]*WL_path[n].c[1]+WL_path[n].c[2]*WL_path[n].c[2]),
                WL[s][t]);
      }
    }
  }

  afree(tmp[0]);
  afree(tmp);
  afree(WL[0]);
  afree(WL);
}

