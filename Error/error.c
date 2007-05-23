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

void error(int test,int no,char *name,char *text)
{
   if (test!=0)
   {
      printf("\n");
      printf("Error in %s\n%s\n",name,text);
      printf("Program aborted\n\n");
      exit(no);
   }
}

