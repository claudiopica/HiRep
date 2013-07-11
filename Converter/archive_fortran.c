/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "suN.h"
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include "communications.h"
#include "moreio.h"
#include "utils.h"
#include "observables.h"


void read_gauge_field_fortran(char filename[]) 
{
  FILE *fp=NULL;
  int g[4];
  int  mu,i;
  struct timeval start, end, etime;
  float test[2*NG*NG];
  float info[16];


  gettimeofday(&start,0);
  
  error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field_fortran",
	"Failed to open file for reading");
  
  error(fread(info,sizeof(float),16,fp)!=16,
        1,"read_gauge_field_fortran",
        "Failed to read header from file");
  
  for(mu=0;mu<4;mu++)	     
    for(g[0]=0;g[0]<GLB_T;g[0]++)
      for(g[3]=0;g[3]<GLB_Z;g[3]++)
	for(g[2]=0;g[2]<GLB_Y;g[2]++)
	  for(g[1]=0;g[1]<GLB_X;g[1]++) {
	    
	    
	    error(fread(test,sizeof(float),2*NG*NG,fp)!= 2*NG*NG,
		  1,"read_gauge_field_asci",
		  "Failed to read header from file");
    
	    int j ;
	    for(i=0;i<NG;i++)
	      for(j=0;j<NG;j++){
		pu_gauge(ipt(g[0],g[1],g[2],g[3]),(mu+1)%4)->c[i+j*NG].re=test[2*(j+i*NG)];
		pu_gauge(ipt(g[0],g[1],g[2],g[3]),(mu+1)%4)->c[i+j*NG].im=test[2*(j+i*NG)+1];
	      }
	    
	  }   
   
   fclose(fp); 
   full_plaquette();
   gettimeofday(&end,0);

   timeval_subtract(&etime,&end,&start);
   lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);
   lprintf("IO",0,"Plaquette eval(%f) read(%f)\n",avr_plaquette(),info[14]);   
}


