/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File error.c
* 
* Error handling functions
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "error.h"
#include "logger.h"
#include "geometry.h"
#include "setup.h"

void error(int test, int no, char *name, char *text)
{
	int glb_test = test;

	/* If this rank is causing a failure, print an error.
	 * Place this first so there is a message in case a bug causes a
	 * deadlock preventing exit. */
	if(test != 0)
	{
		lprintf("ERROR", 0, "%s:\n%s\n", name, text);
		lprintf("ERROR", 0, "Exiting program...\n");
	}

	/* if exit code is negative, then exit immediately, dirtily */
        if (no < 0)
        {
          if (test != 0)
          {
	    exit(0);
	  }
	} else { /* otherwise check other ranks and try to exit cleanly */
#ifdef WITH_MPI
	  /* Co-ordinate with other ranks to decide if any wants to exit */
	  MPI_Allreduce(&test, &glb_test, 1, MPI_INTEGER, MPI_SUM, GLB_COMM);
#endif /* WITH_MPI */

	  if(glb_test != 0)
	  {
            finalize_process();
            exit(no);
          }
        }
}
