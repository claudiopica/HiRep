/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File archive.c
 *
 * Write and read routines for archiving configurations
 *
 *******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include "observables.h"
#include "communications.h"
#include "utils.h"
#include "ranlux.h"


void read_gauge_field_nocheck(char filename[]) 
{
  FILE *fp=NULL;
  int g[4], p[4];
  double *buff=NULL;
  int pid=0;
  int zsize, rz;
  double plaq, testplaq;
  struct timeval start, end, etime;


#ifdef WITH_MPI
  /* MPI variables */
  MPI_Group wg, cg;
  MPI_Status st;
  int cid;
  int mpiret;
#endif

  gettimeofday(&start,0);

  if(PID==0) {
    int d[5]={0}; /* contains NG,GLB_T,GLB_X,GLB_Y,GLB_Z */
    error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field",
        "Failed to open file for reading");
    /* read NG and global size */
    error(fread_BE_int(d,(size_t)(5),fp)!=(5),
        1,"read_gauge_field",
        "Failed to read gauge field geometry");
    /* Check Gauge group and Lattice dimesions */
    if (NG!=d[0]) {
      lprintf("ERROR",0,"Read value of NG [%d] do not match this code [NG=%d].\nPlease recompile.\n",d[0],NG);
      error(1,1,"read_gauge_field " __FILE__,"Gauge group mismatch");
    }
    if (GLB_T!=d[1] ||GLB_X!=d[2] ||GLB_Y!=d[3] ||GLB_Z!=d[4]) {
      lprintf("ERROR",0,"Read value of global lattice size (%d,%d,%d,%d) do not match input file (%d,%d,%d,%d).\n",
          d[1],d[2],d[3],d[4],GLB_T,GLB_X,GLB_Y,GLB_Z);
      error(1,1,"read_gauge_field " __FILE__,"Gauge group mismatch");
    }
    error(fread_BE_double(&plaq,(size_t)(1),fp)!=(1),
        1,"read_gauge_field",
        "Failed to read gauge field plaquette");
  }

#ifdef WITH_MPI
  MPI_Comm_group(GLB_COMM,&wg);
  MPI_Comm_group(cart_comm,&cg);
#endif

  zsize=GLB_Z/NP_Z; rz=GLB_Z-zsize*NP_Z;
  buff=malloc(sizeof(suNg)*4*(GLB_Z/NP_Z+((rz>0)?1:0)));

  g[3]=0;
  for (g[0]=0;g[0]<GLB_T;++g[0]) { /* loop over T, X and Y direction */
    for (g[1]=0;g[1]<GLB_X;++g[1]) {
      for (g[2]=0;g[2]<GLB_Y;++g[2]) {
#ifdef WITH_MPI
        glb_to_proc(g, p); /* get the processor coordinate */
#endif
        for (p[3]=0;p[3]<NP_Z;++p[3]) { /* loop over processors in Z direction */
          int bsize=sizeof(suNg)/sizeof(double)*4*(GLB_Z/NP_Z+((p[3]<rz)?1:0)); /* buffer size in doubles */
#ifdef WITH_MPI
          MPI_Cart_rank(cart_comm, p, &cid);
          MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
          /* read buffer from file */
          if (PID==0) {
            error(fread_BE_double(buff,(size_t)(bsize),fp)!=(bsize),
                1,"read_gauge_field",
                "Failed to read gauge field from file");
          }

#ifdef WITH_MPI
          /* MPI_Barrier(GLB_COMM); */
          if (pid!=0) { /* do send/receive only if the data is not on PID 0 */ 
            /* send buffer */
            if (PID==0) {
              mpiret=MPI_Send(buff, bsize, MPI_DOUBLE, pid, 999, GLB_COMM);
#ifndef NDEBUG
              if (mpiret != MPI_SUCCESS) {
                char mesg[MPI_MAX_ERROR_STRING];
                int mesglen;
                MPI_Error_string(mpiret,mesg,&mesglen);
                lprintf("MPI",0,"ERROR: %s\n",mesg);
                error(1,1,"read_gauge_field " __FILE__,"Cannot send buffer");
              }
#endif
            }
            /* receive buffer */
            if (PID==pid) {
#ifndef NDEBUG
              error(Z!=(GLB_Z/NP_Z+((p[3]<rz)?1:0)),1,"read_gauge_field", "Local lattice size mismatch!");
#endif
              mpiret=MPI_Recv(buff, bsize, MPI_DOUBLE, 0, 999, GLB_COMM, &st);
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
                error(1,1,"read_gauge_field " __FILE__,"Cannot receive buffer");
              }
#endif
            }
          }
#endif

          if (pid==PID) { /* copy buffer into memory */
            int lsite[4];
            suNg *cm;

            /* convert global to local coordinate */
            origin_coord(lsite);
            lsite[0]=g[0]-lsite[0];
            lsite[1]=g[1]-lsite[1];
            lsite[2]=g[2]-lsite[2];

            /* copy buffer in place */
            cm=(suNg*)buff;
            for (lsite[3]=0; lsite[3]<Z; ++lsite[3]) { /* loop on local Z */
              int ix=ipt(lsite[0],lsite[1],lsite[2],lsite[3]);
              suNg *pm=pu_gauge(ix,0);
              *(pm++)=*(cm++); /* copy 4 directions */
              *(pm++)=*(cm++);
              *(pm++)=*(cm++);
              *(pm)=*(cm++);
            }
          }

        } /* end loop over processors in Z direction */
      }
    }
  } /* end loop over T, X and Y direction */

  if (PID==0) fclose(fp); 
  free(buff);

  /* start sendrecv of global gauge field */
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  /* check average plaquette */
  testplaq=avr_plaquette();

  if (PID==0) {
    if (fabs(testplaq-plaq)>1.e-12) {
      lprintf("ERROR",0,"Stored plaquette value [%e] do not match the configuration! [diff=%e]\n",plaq,fabs(testplaq-plaq));
    }
  }

  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec] Plaquette=%e\n",filename,etime.tv_sec,etime.tv_usec,testplaq);

}
