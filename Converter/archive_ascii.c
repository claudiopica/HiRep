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


void read_gauge_field_ascii(char filename[]) 
{
  FILE *fp=NULL;
  int g[4];
  int alpha, gamma, mu;
  int counter=0,counter0=0,pointcounter=0;
  struct timeval start, end, etime;
  suNg tmpmat;
  float re, im;
  int Vdone[GLB_T][GLB_X][GLB_Y][GLB_Z][4];
      

  gettimeofday(&start,0);
  
  error((fp=fopen(filename,"r"))==NULL,1,"read_gauge_field_ascii",
	"Failed to open file for reading");
  
  
  for(g[0]=0;g[0]<GLB_T;g[0]++)
    for(g[1]=0;g[1]<GLB_X;g[1]++)
      for(g[2]=0;g[2]<GLB_Y;g[2]++)
	for(g[3]=0;g[3]<GLB_Z;g[3]++)
	  for(mu=0;mu<4;mu++)
	    Vdone[g[0]][g[1]][g[2]][g[3]][mu]=0;
  
  while(1) {
    /* u( row , col , x , y , z , t , dir) = (re,im) */
    int hm=fscanf(fp," %d %d %d %d\n",
		  &g[0],&g[1],&g[2],&g[3]);
    if(hm != 4) break;
    pointcounter++;
    for(mu=0;mu<4;mu++){
      for(gamma=0;gamma<NG;gamma++)
	for(alpha=0;alpha<NG;alpha++){
	  hm=fscanf(fp," %f %f\n",&re,&im);
	  if(hm != 2) error(0,1,"read_gauge_field_ascii",
			    "Bad number of element in the gauge field\n");
	  tmpmat.c[gamma*NG+alpha].re=re;
	  tmpmat.c[gamma*NG+alpha].im=im;
	  counter++;
	}
      *pu_gauge(ipt(g[0]+1,g[1]-1,g[2]-1,g[3]-1),mu)=tmpmat;
      Vdone[g[0]+1][g[1]-1][g[2]-1][g[3]-1][mu]=1; 
    }
  }
  
  for(g[0]=0;g[0]<GLB_T;g[0]++)
    for(g[1]=0;g[1]<GLB_X;g[1]++)
      for(g[2]=0;g[2]<GLB_Y;g[2]++)
	for(g[3]=0;g[3]<GLB_Z;g[3]++)
	  for(mu=0;mu<4;mu++){
	    if(Vdone[g[0]][g[1]][g[2]][g[3]][mu]==0){
	      counter0+=NG*NG;
	      _suNg_zero(*pu_gauge(ipt(g[0],g[1],g[2],g[3]),mu));
	      }
	  }


  lprintf("IO",0,"Read %d lines\n",counter+pointcounter);
  error(counter+counter0!=NG*NG*4*GLB_T*GLB_X*GLB_Y*GLB_Z,1,"read_gauge_field_ascii " __FILE__,"Bad number of lines in file");
  
  fclose(fp); 
  
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);
  
}

