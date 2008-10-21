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
#include "error.h"
#include "logger.h"

void error(int test,int no,char *name,char *text)
{
   if (test!=0)
   {
      lprintf("ERROR",0,"%s:\n%s\n",name,text);
      lprintf("ERROR",0,"Exiting program...\n");
      exit(no);
   }
}

