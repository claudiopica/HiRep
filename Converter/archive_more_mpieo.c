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
#include "utils.h"
#include "ranlux.h"


void write_gauge_field_mpieo_BE(char filename[]) 
{
  write_gauge_field(filename);
}



void read_gauge_field_mpieo_BE(char filename[]) 
{
  read_gauge_field(filename);
}


void write_gauge_field_mpieo_LE(char filename[]) 
{
  FILE *fp=NULL;
  int g[4], p[4];
  double *buff=NULL;
  int pid=0;
  int zsize, rz;
  double plaq;
  struct timeval start, end, etime;


#ifdef WITH_MPI
  /* MPI variables */
  MPI_Group wg, cg;
  MPI_Status st;
  int cid;
  int mpiret;
#endif

  plaq=avr_plaquette(); /* to use as a checksum in the header */
  if(PID==0) {
    int d[5]={NG,GLB_T,GLB_X,GLB_Y,GLB_Z};
    error((fp=fopen(filename,"wb"))==NULL,1,"write_gauge_field",
        "Failed to open file for writing");
    /* write NG and global size */
    error(fwrite_LE_int(d,(size_t)(5),fp)!=(5),
        1,"write_gauge_field",
        "Failed to write gauge field geometry");
    /* write average plaquette */
    error(fwrite_LE_double(&plaq,(size_t)(1),fp)!=(1),
        1,"write_gauge_field",
        "Failed to write gauge field plaquette");
  } else {
    null_error(); /* fopen */
    null_error(); /* NG, global size */
    null_error(); /* average plaquette */
  }

#ifdef WITH_MPI
  MPI_Comm_group(GLB_COMM,&wg);
  MPI_Comm_group(cart_comm,&cg);
#endif

  gettimeofday(&start,0);

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
          if (pid==PID) { /* fill link buffer */
            int lsite[4];
            suNg *cm;

            /* convert global to local coordinate */
            origin_coord(lsite);
            lsite[0]=g[0]-lsite[0];
            lsite[1]=g[1]-lsite[1];
            lsite[2]=g[2]-lsite[2];

            /* fill buffer */
            cm=(suNg*)buff;
            for (lsite[3]=0; lsite[3]<Z; ++lsite[3]) { /* loop on local Z */
              int ix=ipt(lsite[0],lsite[1],lsite[2],lsite[3]);
              suNg *pm=pu_gauge(ix,0);
              *(cm++)=*(pm++); /* copy 4 directions */
              *(cm++)=*(pm++);
              *(cm++)=*(pm++);
              *(cm++)=*(pm);
            }
          }
#ifdef WITH_MPI
          MPI_Barrier(GLB_COMM); 
          if (pid!=0) { /* do send/receive only if the data is not on PID 0 */ 
            /* send buffer */
            if (pid==PID) {
#ifndef NDEBUG
              error(Z!=(GLB_Z/NP_Z+((p[3]<rz)?1:0)),1,"write_gauge_field", "Local lattice size mismatch!");
#endif
              mpiret=MPI_Send(buff, bsize, MPI_DOUBLE, 0, 999, GLB_COMM);
#ifndef NDEBUG
              if (mpiret != MPI_SUCCESS) {
                char mesg[MPI_MAX_ERROR_STRING];
                int mesglen;
                MPI_Error_string(mpiret,mesg,&mesglen);
                lprintf("MPI",0,"ERROR: %s\n",mesg);
              }
	      error(mpiret != MPI_SUCCESS, 1, "write_gauge_field " __FILE__,
		    "Cannot send buffer");
	    } else {
	      null_error(); /* Z!=(GLB_Z/NP_Z+((p[3]<rz)?1:0)) */
	      null_error(); /* mpiret != MPI_SUCCESS */
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
              }
	      error(mpiret != MPI_SUCCESS, 1, "write_gauge_field " __FILE__,
		    "Cannot receive buffer");
	    } else {
	      null_error();
#endif
            }
          }
#endif

          /* write buffer to file */
          if (PID==0) {
            error(fwrite_LE_double(buff,(size_t)(bsize),fp)!=(bsize),
                1,"write_gauge_field",
                "Failed to write gauge field to file");
	  } else {
	    null_error();
          }

        } /* end loop over processors in Z direction */
      }
    }
  } /* end loop over T, X and Y direction */

  if (PID==0) fclose(fp); 
  free(buff);

  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] saved [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);

}

void read_gauge_field_mpieo_LE(char filename[]) 
{
  FILE *fp=NULL;
  int g[4], p[4];
  double *buff=NULL;
  int pid=0;
  int local_err;
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
    error(fread_LE_int(d,(size_t)(5),fp)!=(5),
        1,"read_gauge_field",
        "Failed to read gauge field geometry");

    /* Check Gauge group and Lattice dimesions */
    local_err = 0;
    if (NG!=d[0]) {
      lprintf("ERROR",0,"Read value of NG [%d] do not match this code [NG=%d].\nPlease recompile.\n",d[0],NG);
      local_err = 1;
    }
    error(local_err, 1, "read_gauge_field " __FILE__, "Gauge group mismatch");

    local_err = 0;
    if (GLB_T!=d[1] ||GLB_X!=d[2] ||GLB_Y!=d[3] ||GLB_Z!=d[4]) {
      lprintf("ERROR",0,"Read value of global lattice size (%d,%d,%d,%d) do not match input file (%d,%d,%d,%d).\n",
          d[1],d[2],d[3],d[4],GLB_T,GLB_X,GLB_Y,GLB_Z);
      local_err = 1;
    }
    error(local_err, 1, "read_gauge_field " __FILE__, "Gauge group mismatch");

    error(fread_LE_double(&plaq,(size_t)(1),fp)!=(1),
        1,"read_gauge_field",
        "Failed to read gauge field plaquette");
  } else {
    null_error(); /* fopen */
    null_error(); /* NG, global size */
    null_error(); /* gauge group */
    null_error(); /* lattice dimensions */
    null_error(); /* read plaquette */
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
            error(fread_LE_double(buff,(size_t)(bsize),fp)!=(bsize),
                1,"read_gauge_field",
                "Failed to read gauge field from file");
	  } else {
	    null_error();
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
              }
	      error(mpiret != MPI_SUCCESS, 1, "read_gauge_field " __FILE__,
		    "Cannot send buffer");
	    } else {
	      null_error();
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
              }
	      error(mpiret != MPI_SUCCESS, 1, "read_gauge_field " __FILE__,
		    "Cannot receive buffer");
	    } else {
	      null_error(); /* local lattice size mismatch */
	      null_error(); /* cannot receive buffer */
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
    local_error = 0;
    if (fabs(testplaq-plaq)>1.e-14) {
      lprintf("ERROR",0,"Stored plaquette value [%e] do not match the configuration! [diff=%e]\n",plaq,fabs(testplaq-plaq));
      local_error = 1;
    }
    error(local_error, 1, "read_gauge_field " __FILE__,
	  "Plaquette value mismatch");
  } else {
    null_error();
  }

  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec] Plaquette=%e\n",filename,etime.tv_sec,etime.tv_usec,testplaq);

}
