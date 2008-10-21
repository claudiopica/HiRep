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
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"


void read_gauge_field(char filename[])
{
  FILE *fp;

  error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field",
        "Failed to open file for reading gauge fields\n");

  fread(u_gauge->ptr,(size_t) sizeof(suNg),(size_t)(VOLUME*4),fp);
  error(ferror(fp)!=0,1,"read_gauge_field",
        "Failed to read gauge field from file");

  fclose(fp);
}

void write_gauge_field(char filename[])
{
  FILE *fp;

  error((fp=fopen(filename,"wb"))==NULL,1,"write_gauge_field",
        "Failed to open file for writing");

  error(fwrite(u_gauge->ptr,(size_t) sizeof(suNg),(size_t)(VOLUME*4),fp) \
        !=(VOLUME*4),1,"write_gauge_field",
        "Failed to write gauge field to file");

  fclose(fp);
}


