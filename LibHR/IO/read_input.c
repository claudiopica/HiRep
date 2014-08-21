/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include <memory.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static void mpi_broadcast_parameters(input_record_t crec[]) {
#ifdef WITH_MPI
  int stringlen=0;
  for (; crec->name!=NULL; ++crec) {
    if(crec->descr==NULL) continue;
    switch(crec->type) {
      case INT_T:
        MPI_Bcast(crec->ptr,1,MPI_INT,0,GLB_COMM);
        break;
      case UNSIGNED_T:
        MPI_Bcast(crec->ptr,1,MPI_UNSIGNED,0,GLB_COMM);
        break;
      case DOUBLE_T:
        MPI_Bcast(crec->ptr,1,MPI_DOUBLE,0,GLB_COMM);
        break;
      case STRING_T:
        stringlen=strlen(crec->ptr)+1;
        MPI_Bcast(&stringlen,1,MPI_INT,0,GLB_COMM);
        MPI_Bcast(crec->ptr,stringlen,MPI_CHAR,0,GLB_COMM);
        break;
      default:
        error(1,1,"read_input " __FILE__, "Unknown type in input descriptor\n");
    }
  }
#endif
  
}

void read_input(input_record_t irec[], char *filename) {
  FILE *inputfile;
  fpos_t pos;
  char buf[256];
  int npar=0, *found=NULL;
  input_record_t *crec=NULL;

  if (irec==NULL || irec->name==NULL) return; /* no input parameter specified */

  /* when using mpi only PID==0 reads the input file */
  if (PID!=0) {
    mpi_broadcast_parameters(irec);
    return;
  }

  error((inputfile=fopen(filename,"r"))==NULL,1,"read_input " __FILE__ ,
      "Failed to open input file\n");

  /* count the input parameters */
  for (crec=irec, npar=0; crec->name!=NULL; ++crec) { ++npar; }
  found=calloc(npar,sizeof(*found)); /* set found to zero */

  do {
    int count=0;
    int nowarn;
    (void) nowarn;
    /*skip white space*/
    nowarn=fscanf(inputfile," ");

    fgetpos(inputfile,&pos);

    /* skip comments */
    nowarn=fscanf(inputfile,"//%n",&count);
    if(count==2) { nowarn=fscanf(inputfile,"%*[^\n]"); goto NEXTLOOP; } if(feof(inputfile)) break; fsetpos(inputfile,&pos);

    /* read variables as described in irec */
    for (count=0; count<npar; ++count) {
      if(irec[count].descr==NULL) continue;
      if(fscanf(inputfile, irec[count].descr, irec[count].ptr)==1) {
        found[count]=1;
        goto NEXTLOOP;
      } 
      if(feof(inputfile)) goto ENDLOOP; 
      fsetpos(inputfile,&pos);
    }

    nowarn=fscanf(inputfile,"%s",buf);
    lprintf("READINPUT",10000,"Ignoring unknown token: [%s]\n",buf);

NEXTLOOP:

    if(ferror(inputfile)) {
      lprintf("READINPUT",0,"Cannot read input file.\n");
      perror(0);
    }

  } while (!feof(inputfile));

ENDLOOP:

  fclose(inputfile);

  while(npar>0) {
    --npar;
    if (found[npar]==0 && irec[npar].descr!=NULL) 
      lprintf("READINPUT",0,
          "Warning: input parameter [%s] not found in [%s]!\n",
          irec[npar].name, filename );
  }
  free(found);

  mpi_broadcast_parameters(irec);

}





