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


void read_gauge_field_milc(char filename[]) 
{
  FILE *fp=NULL;
  int g[4];
  int  mu,i;
  struct timeval start, end, etime;
  float test[2*NG*NG];
  int discard[24];

  gettimeofday(&start,0);
  
  error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field_milc",
	"Failed to open file for reading");
  
  error(fread_LE_int(discard,24,fp)!=24,
        1,"read_gauge_field_asci",
        "Failed to read header from file");
  
   for(g[0]=1;g[0]<GLB_T-1;g[0]++)
     for(g[3]=0;g[3]<GLB_Z;g[3]++)
       for(g[2]=0;g[2]<GLB_Y;g[2]++)
	 for(g[1]=0;g[1]<GLB_X;g[1]++)
	   for(mu=0;mu<4;mu++){	     

	     error(fread(test,sizeof(float),2*NG*NG,fp)!= 2*NG*NG,
		   1,"read_gauge_field_asci",
		   "Failed to read header from file");

	     for(i=0;i<NG*NG;i++) {
	       pu_gauge(ipt(g[0],g[1],g[2],g[3]),(mu+1)%4)->c[i].re=test[2*i];
	       pu_gauge(ipt(g[0],g[1],g[2],g[3]),(mu+1)%4)->c[i].im=test[2*i+1];
	     }
	   }
   
   for(g[3]=0;g[3]<GLB_Z;g[3]++)
     for(g[2]=0;g[2]<GLB_Y;g[2]++)
       for(g[1]=0;g[1]<GLB_X;g[1]++)
	 for(mu=1;mu<4;mu++)
	   *pu_gauge(ipt(GLB_T-1,g[1],g[2],g[3]),mu)=*pu_gauge(ipt(1,g[1],g[2],g[3]),mu);


   for(g[1]=0;g[1]<GLB_X;g[1]++)
     for(g[2]=0;g[2]<GLB_Y;g[2]++)
       for(g[3]=0;g[3]<GLB_Z;g[3]++) {
	 for(mu=0;mu<4;mu++)
	   _suNg_zero(*pu_gauge(ipt(0,g[1],g[2],g[3]),mu));
	 
	 _suNg_zero(*pu_gauge(ipt(GLB_T-1,g[1],g[2],g[3]),0));
       }  
  fclose(fp); 
  full_plaquette();
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);
  
}


void read_gauge_field_milc_no3row(char filename[]) 
{
  FILE *fp=NULL;
  int g[4];
  int  mu,i;
  struct timeval start, end, etime;
  float test[2*(NG)*(NG)];

  gettimeofday(&start,0);
  
  error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field_milc",
	"Failed to open file for reading");
  
  /* error(fread_LE_int(discard,24,fp)!=24, */
  /*       1,"read_gauge_field_asci", */
  /*       "Failed to read header from file"); */
  char line[512],string1[512];
  long cur_pos=0;
  int found =0;
  double head_plaq=0.0;

  while(1){
    error(fgets(line, sizeof(line), fp)==NULL,1,"read_gauge_field_milc_no3row",
	  "Failed to open file for reading");
    if(found==0){
      if(strstr(line,"BEGIN_HEADER")!=NULL) found=1;
      else break;
    }
    if(strstr(line,"PLAQUETTE")!=NULL){
      sscanf(line,"%s %s %lf", string1,string1,&head_plaq);
    }
    cur_pos = ftell(fp);
    lprintf("read_gauge_field_milc_no3row",10,"Header Report: %s",line);
    if(strstr(line,"END_HEADER")!=NULL) {  lprintf("read_gauge_field_milc_no3row",10,"\n",line);; break;}
  }
  fseek(fp,cur_pos,SEEK_SET);
    
   for(g[0]=0;g[0]<GLB_T;g[0]++)
     for(g[3]=0;g[3]<GLB_Z;g[3]++)
       for(g[2]=0;g[2]<GLB_Y;g[2]++)
	 for(g[1]=0;g[1]<GLB_X;g[1]++)
	   for(mu=0;mu<4;mu++){	     

	     error(fread_BE_float(test,(size_t)(12),fp)!= 12,
		   1,"read_gauge_field_asci",
		   "Failed to read header from file");
	     
	     /* for(i=0;i<2*NG;i++)  printf("re=%e im=%f\n",test[2*i],test[2*i+1]); */

	     
	     test[12]= -test[4]*test[8] + test[5]*test[9] + test[2]*test[10] -  test[3]*test[11];
	     test[13]=  test[5]*test[8] + test[4]*test[9] - test[3]*test[10] -  test[2]*test[11];
	     test[14]=  test[4]*test[6] - test[5]*test[7] - test[0]*test[10] +  test[1]*test[11];
	     test[15]= -test[5]*test[6] - test[4]*test[7] + test[1]*test[10] +  test[0]*test[11];
	     test[16]= -test[2]*test[6] + test[3]*test[7] + test[0]*test[8] - test[1]*test[9];
	     test[17]=  test[3]*test[6] + test[2]*test[7] - test[1]*test[8] - test[0]*test[9];

	     double dre=-test[4]*test[8]*test[12] + test[5]*test[9]*test[12] + 
	       test[2]*test[10]*test[12] - test[3]*test[11]*test[12] + 
	       test[5]*test[8]*test[13] + test[4]*test[9]*test[13] - 
	       test[3]*test[10]*test[13] - test[2]*test[11]*test[13] + 
	       test[4]*test[6]*test[14] - test[5]*test[7]*test[14] - 
	       test[0]*test[10]*test[14] + test[1]*test[11]*test[14] - 
	       test[5]*test[6]*test[15] - test[4]*test[7]*test[15] + 
	       test[1]*test[10]*test[15] + test[0]*test[11]*test[15] - 
	       test[2]*test[6]*test[16] + test[3]*test[7]*test[16] + 
	       test[0]*test[8]*test[16] - test[1]*test[9]*test[16] + 
	       test[3]*test[6]*test[17] + test[2]*test[7]*test[17] - 
	       test[1]*test[8]*test[17] - test[0]*test[9]*test[17];

	     double dim=-test[5]*test[8]*test[12] - test[4]*test[9]*test[12] + 
	       test[3]*test[10]*test[12] + test[2]*test[11]*test[12] - 
	       test[4]*test[8]*test[13] + test[5]*test[9]*test[13] + 
	       test[2]*test[10]*test[13] - test[3]*test[11]*test[13] + 
	       test[5]*test[6]*test[14] + test[4]*test[7]*test[14] - 
	       test[1]*test[10]*test[14] - test[0]*test[11]*test[14] + 
	       test[4]*test[6]*test[15] - test[5]*test[7]*test[15] - 
	       test[0]*test[10]*test[15] + test[1]*test[11]*test[15] - 
	       test[3]*test[6]*test[16] - test[2]*test[7]*test[16] + 
	       test[1]*test[8]*test[16] + test[0]*test[9]*test[16] - 
	       test[2]*test[6]*test[17] + test[3]*test[7]*test[17] + 
	       test[0]*test[8]*test[17] - test[1]*test[9]*test[17];


	     error(fabs(dre -1.0) > 1.e-6 || fabs(dim) > 1.e-6, 1,"read_gauge_field_asci",
		   "Failed to reconstruct an unitary matrix");

    
	     for(i=0;i<NG*NG;i++) {
	       pu_gauge(ipt(g[0],g[1],g[2],g[3]),(mu+1)%4)->c[i].re=test[2*i];
	       pu_gauge(ipt(g[0],g[1],g[2],g[3]),(mu+1)%4)->c[i].im=test[2*i+1];
	     }

	   }
  
  fclose(fp); 
  if(head_plaq!=0.0)
    error(fabs(head_plaq-avr_plaquette()) > 1.e-9 , 1,"read_gauge_field_asci",
	"Plaquette mismatch\n");




  full_plaquette();
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);
  
}

