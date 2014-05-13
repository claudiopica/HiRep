/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include "communications.h"
#include "moreio.h"
#include "utils.h"


void read_spinor_field_ascii(char filename[],spinor_field * sf) 
{
  FILE *fp=NULL;
  int g[4];
  int mu;
  int counter=0,counter0=0,pointcounter=0;
  struct timeval start, end, etime;
  float re, im;
  int Vdone[GLB_T][GLB_X][GLB_Y][GLB_Z][12];
  
  
  gettimeofday(&start,0);
  
  error((fp=fopen(filename,"r"))==NULL,1,"read_gauge_field_ascii",
	"Failed to open file for reading");
  
  
  for(g[0]=0;g[0]<GLB_T;g[0]++)
    for(g[1]=0;g[1]<GLB_X;g[1]++)
      for(g[2]=0;g[2]<GLB_Y;g[2]++)
	for(g[3]=0;g[3]<GLB_Z;g[3]++)
	  for(mu=0;mu<12;mu++)
	    Vdone[g[0]][g[1]][g[2]][g[3]][mu]=0;
  
  int linecounter=0;
  while(1) {
    /* u( row , col , x , y , z , t , dir) = (re,im) */
    int hm=fscanf(fp," %d %d %d %d\n",
		  &g[0],&g[1],&g[2],&g[3]);
    if(hm != 4) break;
    pointcounter++;
    linecounter++;
    for(mu=0;mu<12;mu++){
      hm=fscanf(fp," %f %f\n",&re,&im);
      if(hm != 2){
	lprintf("IO",0,"Found problem in file %s at line %d\n",filename,linecounter);
	error(0,1,"read_spinor_field_ascii",
			"Bad number of element in the gauge field\n");
      }

      if((g[0]+1)/T==COORD[0] && (g[1]-1)/X==COORD[1] && (g[2]-1)/Y==COORD[2] && (g[3]-1)/Z==COORD[3] ){
	
	_FIELD_AT(sf,ipt((g[0]+1)%T,(g[1]-1)%X,(g[2]-1)%Y,(g[3]-1)%Z))->c[mu/3].c[mu%3].re=re;
	_FIELD_AT(sf,ipt((g[0]+1)%T,(g[1]-1)%X,(g[2]-1)%Y,(g[3]-1)%Z))->c[mu/3].c[mu%3].im=im;
	Vdone[g[0]+1][g[1]-1][g[2]-1][g[3]-1][mu]=1; 
      }
      
      counter++;
    }
  }
  
  for(g[0]=0;g[0]<GLB_T;g[0]++)
    for(g[1]=0;g[1]<GLB_X;g[1]++)
      for(g[2]=0;g[2]<GLB_Y;g[2]++)
	for(g[3]=0;g[3]<GLB_Z;g[3]++)
	  for(mu=0;mu<12;mu++)
	    if(Vdone[g[0]][g[1]][g[2]][g[3]][mu]!=0)
	      counter0+=1;
  

  global_sum_int(&counter0,1);
  counter0=12*GLB_Y*GLB_X*GLB_Z*GLB_T-counter0;

  lprintf("IO",0,"Read %d lines\n",pointcounter+counter);
  lprintf("IO",0,"counter0= %d/%d\n",counter0,12*GLB_Y*GLB_X*GLB_Z);
  error(counter0!=12*GLB_Y*GLB_X*GLB_Z,1,"read_spinor_field_ascii " __FILE__,"Bad number of lines in file");
  
  fclose(fp); 
  
  start_sf_sendrecv(sf);
  
  complete_sf_sendrecv(sf);

  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);
  
}

