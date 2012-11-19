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
#include "moreio.h"
#include "linear_algebra.h"
#include "spinor_field.h"



void write_spinor_field(char filename[],spinor_field* sp) 
{
  FILE *fp=NULL;
  int g[4], p[4];
  double *buff=NULL;
  int pid=0;
  int zsize, rz;
  double norm2;
  struct timeval start, end, etime;


#ifdef WITH_MPI
  /* MPI variables */
  MPI_Group wg, cg;
  MPI_Status st;
  int cid;
  int mpiret;
#endif

  norm2=spinor_field_sqnorm_f(sp); /* to use as a checksum in the header */
  if(PID==0) {
    int d[6]={NG,NF,GLB_T,GLB_X,GLB_Y,GLB_Z};
    error((fp=fopen(filename,"wb"))==NULL,1,"write_spinor_field",
        "Failed to open file for writing");
    /* write NG, NF and global size */
    error(fwrite(d,(size_t) sizeof(int),(size_t)(6),fp)!=(6),
        1,"write_spinor_field",
        "Failed to write spinor field geometry");
    /* write sq. norm */
    error(fwrite(&norm2,(size_t) sizeof(double),(size_t)(1),fp)!=(1),
        1,"write_spinor_field",
        "Failed to write spinor field plaquette");
  }

#ifdef WITH_MPI
  MPI_Comm_group(GLB_COMM,&wg);
  MPI_Comm_group(cart_comm,&cg);
#endif

  gettimeofday(&start,0);

  zsize=GLB_Z/NP_Z; rz=GLB_Z-zsize*NP_Z;
  buff=malloc(sizeof(suNf_spinor)*(GLB_Z/NP_Z+((rz>0)?1:0)));

  g[3]=0;
  for (g[0]=0;g[0]<GLB_T;++g[0]) { /* loop over T, X and Y direction */
    for (g[1]=0;g[1]<GLB_X;++g[1]) {
      for (g[2]=0;g[2]<GLB_Y;++g[2]) {
#ifdef WITH_MPI
        glb_to_proc(g, p); /* get the processor coordinate */
#endif
        for (p[3]=0;p[3]<NP_Z;++p[3]) { /* loop over processors in Z direction */
          unsigned int bsize=sizeof(suNf_spinor)/sizeof(double)*(GLB_Z/NP_Z+((p[3]<rz)?1:0)); /* buffer size in doubles */
#ifdef WITH_MPI
          MPI_Cart_rank(cart_comm, p, &cid);
          MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
          if (pid==PID) { /* fill spinor buffer */
            int lsite[4];
            suNf_spinor *cm;

            /* convert global to local coordinate */
            origin_coord(lsite);
            lsite[0]=g[0]-lsite[0];
            lsite[1]=g[1]-lsite[1];
            lsite[2]=g[2]-lsite[2];

            /* fill buffer */
            cm=(suNf_spinor*)buff;
            for (lsite[3]=0; lsite[3]<Z; ++lsite[3]) { /* loop on local Z */
              int ix=ipt(lsite[0],lsite[1],lsite[2],lsite[3]);
              suNf_spinor *pm=_FIELD_AT(sp,ix);
              *(cm++)=*(pm++);
            }
          }
#ifdef WITH_MPI
          MPI_Barrier(GLB_COMM); 
          if (pid!=0) { /* do send/receive only if the data is not on PID 0 */ 
            /* send buffer */
            if (pid==PID) {
#ifndef NDEBUG
              error(Z!=(GLB_Z/NP_Z+((p[3]<rz)?1:0)),1,"write_spinor_field", "Local lattice size mismatch!");
#endif
              mpiret=MPI_Send(buff, bsize, MPI_DOUBLE, 0, 999, GLB_COMM);
#ifndef NDEBUG
              if (mpiret != MPI_SUCCESS) {
                char mesg[MPI_MAX_ERROR_STRING];
                int mesglen;
                MPI_Error_string(mpiret,mesg,&mesglen);
                lprintf("MPI",0,"ERROR: %s\n",mesg);
                error(1,1,"write_spinor_field " __FILE__,"Cannot send buffer");
              }
#endif
            }
            /* receive buffer */
            if (PID==0) {
              mpiret=MPI_Recv(buff, bsize, MPI_DOUBLE, pid, 999, GLB_COMM, &st);
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
                error(1,1,"write_spinor_field " __FILE__,"Cannot receive buffer");
              }
#endif
            }
          }
#endif

          /* write buffer to file */
          if (PID==0) {
            error(fwrite(buff,(size_t) sizeof(double),(size_t)(bsize),fp)!=(bsize),
                1,"write_spinor_field",
                "Failed to write spinor field to file");
          }

        } /* end loop over processors in Z direction */
      }
    }
  } /* end loop over T, X and Y direction */

  if (PID==0) fclose(fp); 
  free(buff);

  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Pseudofermion [%s] saved [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);

}


void read_spinor_field(char filename[], spinor_field *sp) 
{
  FILE *fp=NULL;
  int g[4], p[4];
  double *buff=NULL;
  int pid=0;
  int zsize, rz;
  double norm2, testnorm2;
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
    int d[6]={0}; /* contains NG,NF,GLB_T,GLB_X,GLB_Y,GLB_Z */
    error((fp=fopen(filename,"rb"))==NULL,1,"read_spinor_field",
        "Failed to open file for reading");
    /* read NG, NF and global size */
    error(fread(d,(size_t) sizeof(int),(size_t)(6),fp)!=(6),
        1,"read_spinor_field",
        "Failed to read spinor field geometry");
    /* Check Gauge group and Lattice dimesions */
    if (NG!=d[0]) {
      lprintf("ERROR",0,"Read value of NG [%d] do not match this code [NG=%d].\nPlease recompile.\n",d[0],NG);
      error(1,1,"read_spinor_field " __FILE__,"Gauge group mismatch");
    }
    if (NF!=d[1]) {
      lprintf("ERROR",0,"Read value of NF [%d] do not match this code [NF=%d].\nPlease recompile.\n",d[1],NF);
      error(1,1,"read_spinor_field " __FILE__,"Representation mismatch");
    }
    if (GLB_T!=d[2] ||GLB_X!=d[3] ||GLB_Y!=d[4] ||GLB_Z!=d[5]) {
      lprintf("ERROR",0,"Read value of global lattice size (%d,%d,%d,%d) do not match input file (%d,%d,%d,%d).\n",
          d[2],d[3],d[4],d[5],GLB_T,GLB_X,GLB_Y,GLB_Z);
      error(1,1,"read_spinor_field " __FILE__,"Lattice mismatch");
    }
    error(fread(&norm2,(size_t) sizeof(double),(size_t)(1),fp)!=(1),
        1,"read_spinor_field",
        "Failed to read spinor field sq. norm");
  }

#ifdef WITH_MPI
  MPI_Comm_group(GLB_COMM,&wg);
  MPI_Comm_group(cart_comm,&cg);
#endif

  zsize=GLB_Z/NP_Z; rz=GLB_Z-zsize*NP_Z;
  buff=malloc(sizeof(suNf_spinor)*(GLB_Z/NP_Z+((rz>0)?1:0)));

  g[3]=0;
  for (g[0]=0;g[0]<GLB_T;++g[0]) { /* loop over T, X and Y direction */
    for (g[1]=0;g[1]<GLB_X;++g[1]) {
      for (g[2]=0;g[2]<GLB_Y;++g[2]) {
#ifdef WITH_MPI
        glb_to_proc(g, p); /* get the processor coordinate */
#endif
        for (p[3]=0;p[3]<NP_Z;++p[3]) { /* loop over processors in Z direction */
          unsigned int bsize=sizeof(suNf_spinor)/sizeof(double)*(GLB_Z/NP_Z+((p[3]<rz)?1:0)); /* buffer size in doubles */
#ifdef WITH_MPI
          MPI_Cart_rank(cart_comm, p, &cid);
          MPI_Group_translate_ranks(cg, 1, &cid, wg, &pid);
#endif
          /* read buffer from file */
          if (PID==0) {
            error(fread(buff,(size_t) sizeof(double),(size_t)(bsize),fp)!=(bsize),
                1,"read_spinor_field",
                "Failed to read spinor field from file");
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
                error(1,1,"read_spinor_field " __FILE__,"Cannot send buffer");
              }
#endif
            }
            /* receive buffer */
            if (PID==pid) {
#ifndef NDEBUG
              error(Z!=(GLB_Z/NP_Z+((p[3]<rz)?1:0)),1,"read_spinor_field", "Local lattice size mismatch!");
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
                error(1,1,"read_spinor_field " __FILE__,"Cannot receive buffer");
              }
#endif
            }
          }
#endif

          if (pid==PID) { /* copy buffer into memory */
            int lsite[4];
            suNf_spinor *cm;

            /* convert global to local coordinate */
            origin_coord(lsite);
            lsite[0]=g[0]-lsite[0];
            lsite[1]=g[1]-lsite[1];
            lsite[2]=g[2]-lsite[2];

            /* copy buffer in place */
            cm=(suNf_spinor*)buff;
            for (lsite[3]=0; lsite[3]<Z; ++lsite[3]) { /* loop on local Z */
              int ix=ipt(lsite[0],lsite[1],lsite[2],lsite[3]);
              suNf_spinor *pm=_FIELD_AT(sp,ix);
              *(pm++)=*(cm++);
            }
          }

        } /* end loop over processors in Z direction */
      }
    }
  } /* end loop over T, X and Y direction */

  if (PID==0) fclose(fp); 
  free(buff);

  /* check sq. norm */
  testnorm2=spinor_field_sqnorm_f(sp);
  if (PID==0) {
    if (fabs(testnorm2-norm2)>1.e-14) {
      lprintf("ERROR",0,"Stored sq. norm value [%e] do not match the pseudofermion! [diff=%e]\n",norm2,fabs(testnorm2-norm2));
      error(1,1,"read_spinor_field " __FILE__,"Sq. norm value mismatch");
    }
  }

  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Pseudofermion [%s] read [%ld sec %ld usec] Sq. norm=%e\n",filename,etime.tv_sec,etime.tv_usec,testnorm2);

}


