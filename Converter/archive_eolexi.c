/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include "communications.h"
#include "moreio.h"
#include "observables.h"


static int index_eolexi(int x0,int x1,int x2,int x3)
{
   int y0,y1,y2,y3;
   int ix;
   
   y0=x0%GLB_T;
   y1=x1%GLB_X;
   y2=x2%GLB_Y;
   y3=x3%GLB_Z;

   ix = (y3+GLB_Z*(y2+GLB_Y*(y1+GLB_X*y0)))/2;
   if((y0+y1+y2+y3)&1){
      ix+=(GLB_T*GLB_X*GLB_Y*GLB_Z)/2;
   }

   return ix;
      
}


void write_gauge_field_eolexi_BE(char filename[]) 
{
  FILE *fp=NULL;
  int pid=0;
  int g[4], c[4];
  int mu;
  int ix_eolexi, ix_mpieo;

#ifdef WITH_MPI
  /* MPI variables */
  int p[4];
  MPI_Group wg, cg;
  MPI_Status st;
  int cid;
  int mpiret;
#endif
  
  suNg *eolexi_field=NULL;
  
  if(PID==0) {
    eolexi_field=malloc(sizeof(suNg)*4*GLB_T*GLB_X*GLB_Y*GLB_Z);
  }

#ifdef WITH_MPI
  MPI_Comm_group(GLB_COMM,&wg);
  MPI_Comm_group(cart_comm,&cg);
#endif


  for(g[0]=0; g[0]<GLB_T; g[0]++)
  for(g[1]=0; g[1]<GLB_X; g[1]++)
  for(g[2]=0; g[2]<GLB_Y; g[2]++)
  for(g[3]=0; g[3]<GLB_Z; g[3]++) {
  
#ifdef WITH_MPI
    glb_to_proc(g, p); /* get the processor coordinate */
    MPI_Cart_rank(cart_comm, p, &cid);
    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif

    c[0]=g[0]%T; c[1]=g[1]%X; c[2]=g[2]%Y; c[3]=g[3]%Z;
    ix_eolexi=index_eolexi(g[0],g[1],g[2],g[3]);
    ix_mpieo=ipt(c[0],c[1],c[2],c[3]);

    if(pid==0) {
      for(mu=0; mu<4;mu++)
        eolexi_field[4*ix_eolexi+mu]=*pu_gauge(ix_mpieo,mu);
    }

#ifdef WITH_MPI
    MPI_Barrier(GLB_COMM); 

    if(pid!=0) {

      /* send buffer */
      if(PID==pid) {

        mpiret=MPI_Send(pu_gauge(ix_mpieo,0), 2*4*NG*NG, MPI_DOUBLE, pid, 999, GLB_COMM);
#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          error(1,1,"write_gauge_field_eolexi " __FILE__,"Cannot send u_gauge buffer");
        }
#endif
      }
      /* receive buffer */
      if (PID==0) {
        mpiret=MPI_Recv(eolexi_field+4*ix_eolexi, 2*4*NG*NG, MPI_DOUBLE, 0, 999, GLB_COMM, &st);
#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          if (st.MPI_ERROR != MPI_SUCCESS) {
            MPI_Error_string(st.MPI_ERROR,mesg,&mesglen);
            lprintf("MPI",0,"Source [%d] Tag [%] ERROR: %s\n",
                st.MPI_SOURCE,
                st.MPI_TAG,
                mesg);
          }
          error(1,1,"read_gauge_field_eolexi " __FILE__,"Cannot receive u_gauge buffer");
        }
#endif
      }
    }
#endif /* WITH_MPI */

  }

  if(PID==0) {
    error((fp=fopen(filename,"wb"))==NULL,1,"write_gauge_field_eolexi",
        "Failed to open file for writing");
    error(fwrite_BE_double((double*)eolexi_field,4*GLB_T*GLB_X*GLB_Y*GLB_Z*sizeof(suNg)/sizeof(double),fp)!=4*GLB_T*GLB_X*GLB_Y*GLB_Z*(int)sizeof(suNg)/(int)sizeof(double),
        1,"write_gauge_field_eolexi",
        "Failed to write gauge field to file");
    fclose(fp);
    free(eolexi_field);
  }

  lprintf("IO",0,"Configuration [%s] saved\n",filename);
}



void read_gauge_field_eolexi_BE(char filename[]) 
{
  FILE *fp=NULL;
  int g[4], c[4];
  int mu;
  int ix_eolexi, ix_mpieo;
  double plaq;
 
  suNg *eolexi_field = malloc(sizeof(suNg)*4*GLB_T*GLB_X*GLB_Y*GLB_Z);

  if(PID==0) {
    error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field_eolexi",
        "Failed to open file for reading");
    error(fread_BE_double((double*)eolexi_field,4*GLB_T*GLB_X*GLB_Y*GLB_Z*sizeof(suNg)/sizeof(double),fp)!=4*GLB_T*GLB_X*GLB_Y*GLB_Z*(int)sizeof(suNg)/(int)sizeof(double),
        1,"read_gauge_field_eolexi",
        "Failed to read gauge field from file");
  }

#ifdef WITH_MPI
  bcast((double*)eolexi_field,(sizeof(suNg)/sizeof(double))*4*GLB_T*GLB_X*GLB_Y*GLB_Z);
#endif

  for(c[0]=0; c[0]<T; c[0]++)
  for(c[1]=0; c[1]<X; c[1]++)
  for(c[2]=0; c[2]<Y; c[2]++)
  for(c[3]=0; c[3]<Z; c[3]++) {
  
    g[0]=zerocoord[0]+c[0];
    g[1]=zerocoord[1]+c[1];
    g[2]=zerocoord[2]+c[2];
    g[3]=zerocoord[3]+c[3];

    ix_eolexi=index_eolexi(g[0],g[1],g[2],g[3]);
    ix_mpieo=ipt(c[0],c[1],c[2],c[3]);
 
    for(mu=0; mu<4;mu++)
      *pu_gauge(ix_mpieo,mu) = eolexi_field[4*ix_eolexi+mu];
    
  }
  
  if (PID==0) fclose(fp); 

  /* start sendrecv of global gauge field */
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  plaq=avr_plaquette();
  lprintf("IO",0,"Configuration [%s] read.  Plaquette=%e\n",filename,plaq);

  if(PID==0) {
    int i,j;
    int ix=ipt(0,0,0,0);
    for(mu=0;mu<4;mu++) {
      lprintf("IO",20,"x=(0,0,0,0) mu=%d pu_gauge =\n",mu);
      for(i=0; i<NG; i++) {
        lprintf("IO",20,"[ ");
        for(j=0; j<NG; j++){
#ifdef GAUGE_SON
          lprintf("IO",20,"(%.2f) ",pu_gauge(ix,mu)->c[i*NG+j]);
#else
          lprintf("IO",20,"(%.2f , %.2f) ",pu_gauge(ix,mu)->c[i*NG+j].re,pu_gauge(ix,mu)->c[i*NG+j].im);
#endif //GAUGE_SON
	}
        lprintf("IO",20,"]\n");
      }
    }
    free(eolexi_field);
  }

}


void write_gauge_field_eolexi_LE(char filename[]) 
{
  FILE *fp=NULL;
  int pid=0;
  int g[4], c[4];
  int mu;
  int ix_eolexi, ix_mpieo;

#ifdef WITH_MPI
  /* MPI variables */
  int p[4];
  MPI_Group wg, cg;
  MPI_Status st;
  int cid;
  int mpiret;
#endif
  
  suNg *eolexi_field=NULL;
  
  if(PID==0) {
    eolexi_field=malloc(sizeof(suNg)*4*GLB_T*GLB_X*GLB_Y*GLB_Z);
  }

#ifdef WITH_MPI
  MPI_Comm_group(GLB_COMM,&wg);
  MPI_Comm_group(cart_comm,&cg);
#endif


  for(g[0]=0; g[0]<GLB_T; g[0]++)
  for(g[1]=0; g[1]<GLB_X; g[1]++)
  for(g[2]=0; g[2]<GLB_Y; g[2]++)
  for(g[3]=0; g[3]<GLB_Z; g[3]++) {
  
#ifdef WITH_MPI
    glb_to_proc(g, p); /* get the processor coordinate */
    MPI_Cart_rank(cart_comm, p, &cid);
    MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif

    c[0]=g[0]%T; c[1]=g[1]%X; c[2]=g[2]%Y; c[3]=g[3]%Z;
    ix_eolexi=index_eolexi(g[0],g[1],g[2],g[3]);
    ix_mpieo=ipt(c[0],c[1],c[2],c[3]);

    if(pid==0) {
      for(mu=0; mu<4;mu++)
        eolexi_field[4*ix_eolexi+mu]=*pu_gauge(ix_mpieo,mu);
    }

#ifdef WITH_MPI
    MPI_Barrier(GLB_COMM); 

    if(pid!=0) {

      /* send buffer */
      if(PID==pid) {

        mpiret=MPI_Send(pu_gauge(ix_mpieo,0), 2*4*NG*NG, MPI_DOUBLE, pid, 999, GLB_COMM);
#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          error(1,1,"write_gauge_field_eolexi " __FILE__,"Cannot send u_gauge buffer");
        }
#endif
      }
      /* receive buffer */
      if (PID==0) {
        mpiret=MPI_Recv(eolexi_field+4*ix_eolexi, 2*4*NG*NG, MPI_DOUBLE, 0, 999, GLB_COMM, &st);
#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          if (st.MPI_ERROR != MPI_SUCCESS) {
            MPI_Error_string(st.MPI_ERROR,mesg,&mesglen);
            lprintf("MPI",0,"Source [%d] Tag [%] ERROR: %s\n",
                st.MPI_SOURCE,
                st.MPI_TAG,
                mesg);
          }
          error(1,1,"read_gauge_field_eolexi " __FILE__,"Cannot receive u_gauge buffer");
        }
#endif
      }
    }
#endif /* WITH_MPI */

  }

  if(PID==0) {
    error((fp=fopen(filename,"wb"))==NULL,1,"write_gauge_field_eolexi",
        "Failed to open file for writing");
    error(fwrite_LE_double((double*)eolexi_field,4*GLB_T*GLB_X*GLB_Y*GLB_Z*sizeof(suNg)/sizeof(double),fp)!=4*GLB_T*GLB_X*GLB_Y*GLB_Z*(int)sizeof(suNg)/(int)sizeof(double),
        1,"write_gauge_field_eolexi",
        "Failed to write gauge field to file");
    fclose(fp);
    free(eolexi_field);
  }

  lprintf("IO",0,"Configuration [%s] saved\n",filename);
}



void read_gauge_field_eolexi_LE(char filename[]) 
{
  FILE *fp=NULL;
  int g[4], c[4];
  int mu;
  int ix_eolexi, ix_mpieo;
  double plaq;
 
  suNg *eolexi_field = malloc(sizeof(suNg)*4*GLB_T*GLB_X*GLB_Y*GLB_Z);

  if(PID==0) {
    error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field_eolexi",
        "Failed to open file for reading");
    error(fread_LE_double((double*)eolexi_field,4*GLB_T*GLB_X*GLB_Y*GLB_Z*sizeof(suNg)/sizeof(double),fp)!=4*GLB_T*GLB_X*GLB_Y*GLB_Z*(int)sizeof(suNg)/(int)sizeof(double),
        1,"read_gauge_field_eolexi",
        "Failed to read gauge field from file");
  }

#ifdef WITH_MPI
  bcast((double*)eolexi_field,(sizeof(suNg)/sizeof(double))*4*GLB_T*GLB_X*GLB_Y*GLB_Z);
#endif

  for(c[0]=0; c[0]<T; c[0]++)
  for(c[1]=0; c[1]<X; c[1]++)
  for(c[2]=0; c[2]<Y; c[2]++)
  for(c[3]=0; c[3]<Z; c[3]++) {
  
    g[0]=zerocoord[0]+c[0];
    g[1]=zerocoord[1]+c[1];
    g[2]=zerocoord[2]+c[2];
    g[3]=zerocoord[3]+c[3];

    ix_eolexi=index_eolexi(g[0],g[1],g[2],g[3]);
    ix_mpieo=ipt(c[0],c[1],c[2],c[3]);
 
    for(mu=0; mu<4;mu++)
      *pu_gauge(ix_mpieo,mu) = eolexi_field[4*ix_eolexi+mu];
    
  }
  
  if (PID==0) fclose(fp); 

  /* start sendrecv of global gauge field */
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  plaq=avr_plaquette();
  lprintf("IO",0,"Configuration [%s] read.  Plaquette=%e\n",filename,plaq);

  if(PID==0) {
    int i,j;
    int ix=ipt(0,0,0,0);
    for(mu=0;mu<4;mu++) {
      lprintf("IO",20,"x=(0,0,0,0) mu=%d pu_gauge =\n",mu);
      for(i=0; i<NG; i++) {
        lprintf("IO",20,"[ ");
        for(j=0; j<NG; j++){
#ifdef GAUGE_SON
          lprintf("IO",20,"(%.2f) ",pu_gauge(ix,mu)->c[i*NG+j]);
#else
          lprintf("IO",20,"(%.2f , %.2f) ",pu_gauge(ix,mu)->c[i*NG+j].re,pu_gauge(ix,mu)->c[i*NG+j].im);
#endif //GAUGE_SON
	}
        lprintf("IO",20,"]\n");
      }
    }
    free(eolexi_field);
  }

}



