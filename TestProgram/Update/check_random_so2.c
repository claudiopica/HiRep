/******************************************************************************* 
* Check that the molecular dynamics evolution is reversible
* 
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "update.h"
#include "global.h"
#include "logger.h"

int main(int argc,char *argv[])
{


for (int i=0;i<100000;i++){
   lprintf("THETA",0," %f \n", random_so2(3.56531));
}


   lprintf("THETA",0," \n");
   exit(0);
}
