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
#include "utils.h"

/*
u(row,col,x,y,z,t,dir)
*/

void read_gauge_field_henty(char filename[]) 
{
  FILE *fp=NULL;
  int g[5], p[4], c[4];
  int pid=0;
  int row, col;
  int maxrow, maxcol;
  int counter;
  struct timeval start, end, etime;
  suNg tmpmat;


#ifdef WITH_MPI
  /* MPI variables */
  MPI_Group wg, cg;
  MPI_Status st;
  int cid;
  int mpiret;
#endif

  gettimeofday(&start,0);

  if(PID==0) {
    error((fp=fopen(filename,"r"))==NULL,1,"read_gauge_field_for_henty",
        "Failed to open file for reading");
  }

#ifdef WITH_MPI
  MPI_Comm_group(MPI_COMM_WORLD,&wg);
  MPI_Comm_group(cart_comm,&cg);
#endif

  maxrow=maxcol=0;
  p[0]=p[1]=p[2]=p[3]=0;
  counter=0;
  while(1) {
    if (PID==0) {
      float re, im;
      /* u( row , col , x , y , z , t , dir) = (re,im) */
      int hm=fscanf(fp," u( %d , %d , %d , %d , %d , %d , %d ) =  (%f,%f)\n",
        &row,&col,&g[1],&g[2],&g[3],&g[0],&g[4],&re,&im);
      if(hm != 9) break;
      counter++;
      tmpmat.c[row*NG+col].re=(double)re;
      tmpmat.c[row*NG+col].im=(double)im;
      if(row>=maxrow) maxrow=row;
      else error(maxrow!=NG-1,1,"read_gauge_field_for_henty",
        "Bad number of rows");
      if(col>=maxcol) maxcol=col;
      else error(maxcol!=NG-1,1,"read_gauge_field_for_henty",
        "Bad number of columns");
    }
    if(row==NG-1 && col==NG-1) {
#ifdef WITH_MPI
      glb_to_proc(g, p); /* get the processor coordinate */
      MPI_Cart_rank(cart_comm, p, &cid);
      MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
      if(pid == 0) {
        c[0]=g[0]%T; c[1]=g[1]%X; c[2]=g[2]%Y; c[3]=g[3]%Z;
        *pu_gauge(ipt(c[0],c[1],c[2],c[3]),g[4])=tmpmat;
      }
      
#ifdef WITH_MPI
      MPI_Barrier(MPI_COMM_WORLD); 

      if(pid != 0) {

        /* send buffer */
        if (PID==0) {
          mpiret=MPI_Send(g, 5, MPI_INT, pid, 998, MPI_COMM_WORLD);
#ifndef NDEBUG
          if (mpiret != MPI_SUCCESS) {
            char mesg[MPI_MAX_ERROR_STRING];
            int mesglen;
            MPI_Error_string(mpiret,mesg,&mesglen);
            lprintf("MPI",0,"ERROR: %s\n",mesg);
            error(1,1,"read_gauge_field_for_henty " __FILE__,"Cannot send (coord,mu) buffer");
          }
#endif
          mpiret=MPI_Send(&tmpmat, 2*NG*NG, MPI_DOUBLE, pid, 999, MPI_COMM_WORLD);
#ifndef NDEBUG
          if (mpiret != MPI_SUCCESS) {
            char mesg[MPI_MAX_ERROR_STRING];
            int mesglen;
            MPI_Error_string(mpiret,mesg,&mesglen);
            lprintf("MPI",0,"ERROR: %s\n",mesg);
            error(1,1,"read_gauge_field_for_henty " __FILE__,"Cannot send u_gauge buffer");
          }
#endif
        }
        /* receive buffer */
        if (PID==pid) {
          mpiret=MPI_Recv(g, 5, MPI_INT, 0, 998, MPI_COMM_WORLD, &st);
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
            error(1,1,"read_gauge_field_for_henty " __FILE__,"Cannot receive (coord,mu) buffer");
          }
#endif
          mpiret=MPI_Recv(&tmpmat, 2*NG*NG, MPI_DOUBLE, 0, 999, MPI_COMM_WORLD, &st);
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
            error(1,1,"read_gauge_field_for_henty " __FILE__,"Cannot receive u_gauge buffer");
          }
#endif
          c[0]=g[0]%T; c[1]=g[1]%X; c[2]=g[2]%Y; c[3]=g[3]%Z;
          *pu_gauge(ipt(c[0],c[1],c[2],c[3]),g[4])=tmpmat;
        }
      }
#endif /* WITH_MPI */
    }
  }

  lprintf("IO",0,"Read %d lines\n",counter);
  error(counter!=NG*NG*4*T*X*Y*Z,1,"read_gauge_field_for_henty " __FILE__,"Bad number of lines in file");

  if (PID==0) fclose(fp); 

  /* start sendrecv of global gauge field */
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);

}

