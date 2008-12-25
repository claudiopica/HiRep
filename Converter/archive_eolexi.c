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


void write_gauge_field_eolexi(char filename[]) 
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
  MPI_Comm_group(MPI_COMM_WORLD,&wg);
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

    if(pid == 0) {
      for(mu=0; mu<4;mu++)
        eolexi_field[4*ix_eolexi+mu]=*pu_gauge(ipt(c[0],c[1],c[2],c[3]),mu);
    }

#ifdef WITH_MPI
    MPI_Barrier(MPI_COMM_WORLD); 

    if(pid != 0) {

      /* send buffer */
      if (PID==pid) {

        mpiret=MPI_Send(pu_gauge(ipt(c[0],c[1],c[2],c[3]),0), 2*4*NG*NG, MPI_DOUBLE, pid, 999, MPI_COMM_WORLD);
#ifndef NDEBUG
        if (mpiret != MPI_SUCCESS) {
          char mesg[MPI_MAX_ERROR_STRING];
          int mesglen;
          MPI_Error_string(mpiret,mesg,&mesglen);
          lprintf("MPI",0,"ERROR: %s\n",mesg);
          error(1,1,"read_gauge_field_eolexi " __FILE__,"Cannot send u_gauge buffer");
        }
#endif
      }
      /* receive buffer */
      if (PID==0) {
        mpiret=MPI_Recv(eolexi_field+4*ix_eolexi, 2*4*NG*NG, MPI_DOUBLE, 0, 999, MPI_COMM_WORLD, &st);
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
    error(fwrite(eolexi_field,sizeof(suNg),4*GLB_T*GLB_X*GLB_Y*GLB_Z,fp)!=4*GLB_T*GLB_X*GLB_Y*GLB_Z,
        1,"write_gauge_field_eolexi",
        "Failed to write gauge field to file");
    fclose(fp);
  }
}



void read_gauge_field_eolexi(char filename[]) 
{
  FILE *fp=NULL;
  int g[4], c[4];
  int mu;
  int ix_eolexi, ix_mpieo;
 
  suNg *eolexi_field = malloc(sizeof(suNg)*4*GLB_T*GLB_X*GLB_Y*GLB_Z);

  if(PID==0) {
    error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field_eolexi",
        "Failed to open file for reading");
    error(fread(eolexi_field,sizeof(suNg),4*GLB_T*GLB_X*GLB_Y*GLB_Z,fp)!=4*GLB_T*GLB_X*GLB_Y*GLB_Z,
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
  
    g[0]=COORD[0]*T+c[0];
    g[1]=COORD[1]*X+c[1];
    g[2]=COORD[2]*Y+c[2];
    g[3]=COORD[3]*Z+c[3];

    ix_eolexi=index_eolexi(g[0],g[1],g[2],g[3]);
    ix_mpieo=ipt(c[0],c[1],c[2],c[3]);
 
    for(mu=0; mu<4;mu++)
      *pu_gauge(ix_mpieo,mu) = eolexi_field[4*ix_eolexi+mu];
    
  }
  
  if (PID==0) fclose(fp); 

  /* start sendrecv of global gauge field */
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  free(eolexi_field);

  lprintf("IO",0,"Configuration [%s] read\n",filename);

}


