/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File print_wdmatrix.c
*
* print pbm image of the Wilson-Dirac matrix
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "geometry.h"
#include "global.h"
#include "error.h"

/* write a .pbm file containing an image of the wilson-dirac matrix:
 * every black dot represents a non-zero entry of the matrix in (x,mu) space
 * (flavor in not considered here)
 */
void print_wdmatrix(char *filename)
{
  int x0,x1,x2;
	FILE *fp;

  error((fp=fopen(filename,"w"))==NULL,1,"print_wdmatrix",
        "Failed to open file for writing WD matrix image\n");

	fprintf(fp,"P1\n%ld %ld\n",VOLUME,VOLUME);
  for (x0=0;x0<VOLUME;x0++){ /* loop over rows */
    int nei[8];
		for (x1=0;x1<4;++x1){
			nei[2*x1]=iup(x0,x1);
			nei[2*x1+1]=idn(x0,x1);
		}
		for(x2=0;x2<VOLUME;++x2){ /* loop over columns */
			char c='0';
			if(x2==x0) c='1';
			else
				for(x1=0;x1<8;++x1){
					if (x2==nei[x1]) c='1';
				}
			fprintf(fp,"%c ",c);
		}
		fprintf(fp,"\n");
  }  

	fclose (fp);
}
