/*******************************************************************************
*
* Testing geometry
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "geometry.h"
#include "global.h"


int main(int argc,char *argv[])
{
	geometry_mpi_eo();
	printf("The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
	printf("\n");
	
	test_geometry_mpi_eo();

	exit(0);
}
