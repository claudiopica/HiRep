/*************************************************************************** \
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
#include "linear_algebra.h"


void read_gauge_field_openQCD(char filename[]) 
{
  FILE *fp=NULL;
  int g[4];
  int  mu,i,j;
  struct timeval start, end, etime;
  double test[2*NG*NG];
  int size[4];
  double readplaq;

  gettimeofday(&start,0);
  
  error((fp=fopen(filename,"rb"))==NULL,1,"read_gauge_field_openQCD",
	"Failed to open file for reading");
  
  error(fread_LE_int(size,4,fp)!=4,
        1,"read_gauge_field_openQCD",
        "Failed to read lattice size from the header of the conf file");

  error(fread(&readplaq,sizeof(double),1,fp)!=1,1,"read_gauge_field_openQCD",
        "Failed to read the plaquette value from the header  of the conf file");

  int id;
  error(size[0]!=GLB_T,1,"read_gauge_field_openQCD","Wrong lattice size");
  error(size[1]!=GLB_X,1,"read_gauge_field_openQCD","Wrong lattice size");
  error(size[2]!=GLB_Y,1,"read_gauge_field_openQCD","Wrong lattice size");
  error(size[3]!=GLB_Z,1,"read_gauge_field_openQCD","Wrong lattice size");
  


  for(g[0]=0;g[0]<GLB_T;g[0]++)
    for(g[1]=0;g[1]<GLB_X;g[1]++)
      for(g[2]=0;g[2]<GLB_Y;g[2]++)
        for(g[3]=0;g[3]<GLB_Z;g[3]++)
          if((g[0]+g[1]+g[2]+g[3])%2==1)
            for(mu=0;mu<4;mu++){	     
              
              id = ipt(g[0],g[1],g[2],g[3]);
              
              error(fread(test,sizeof(double),2*NG*NG,fp)!= 2*NG*NG,
                    1,"read_gauge_field_openQCD",
                    "Failed to read header from file");               
              for(j=0;j<NG;j++) 
                for(i=0;i<NG;i++) {
                  int k = j+i*NG;
                  pu_gauge(id,mu)->c[k].re=test[2*k];
                  pu_gauge(id,mu)->c[k].im=test[2*k+1];
                }
              
              
              
              id = idn(id,mu);
              
              error(fread(test,sizeof(double),2*NG*NG,fp)!= 2*NG*NG,
                    1,"read_gauge_field_openQCD",
                    "Failed to read header from file");
              for(j=0;j<NG;j++) 
                for(i=0;i<NG;i++) {
                  int k = j+i*NG;
                   pu_gauge(id,mu)->c[k].re=test[2*k];
                  pu_gauge(id,mu)->c[k].im=test[2*k+1];
                }
              
              
            }
 
  
  fclose(fp); 

  if((NG*avr_plaquette()-readplaq)*(NG*avr_plaquette()-readplaq)>1.e-14)
    error(1,1,"read_gauge_field_openQCD",
                    "Wrong plaquette checksum");
  else
    lprintf("IO",0,"Plaquette checksum matches\n");
  
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] read [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);
  
}

void write_gauge_field_openQCD(char filename[]) 
{
  FILE *fp=NULL;
  int g[4];
  int  mu,i,j;
  struct timeval start, end, etime;
  double test[2*NG*NG];
  int size[4]={GLB_T,GLB_X,GLB_Y,GLB_Z};
  double writeplaq=NG*avr_plaquette();

  gettimeofday(&start,0);
  
  error((fp=fopen(filename,"wb"))==NULL,1,"write_gauge_field_openQCD",
	"Failed to open file for writing");
  
  error(fwrite_LE_int(size,4,fp)!=4,
        1,"write_gauge_field_openQCD",
        "Failed to write lattice size into the header of the conf file");

  error(fwrite_LE_double(&writeplaq,1,fp)!=1,1,"write_gauge_field_openQCD",
        "Failed to write the plaquette value into the header of the conf file");

  int id;

  for(g[0]=0;g[0]<GLB_T;g[0]++)
    for(g[1]=0;g[1]<GLB_X;g[1]++)
      for(g[2]=0;g[2]<GLB_Y;g[2]++)
        for(g[3]=0;g[3]<GLB_Z;g[3]++)
          if((g[0]+g[1]+g[2]+g[3])%2==1)
            for(mu=0;mu<4;mu++){	     
              
              id = ipt(g[0],g[1],g[2],g[3]);
              
             for(j=0;j<NG;j++) 
               for(i=0;i<NG;i++) {
                 int k = j+i*NG;
                 test[2*k]=pu_gauge(id,mu)->c[k].re;
                 test[2*k+1]=pu_gauge(id,mu)->c[k].im;
               }
             
             error(fwrite_LE_double(test,2*NG*NG,fp)!= 2*NG*NG,
                   1,"read_gauge_field_openQCD",
                   "Failed to read header from file");               
             
             
             id = idn(id,mu);
              
             for(j=0;j<NG;j++) 
               for(i=0;i<NG;i++) {
                 int k = j+i*NG;
                 test[2*k]=pu_gauge(id,mu)->c[k].re;
                 test[2*k+1]=pu_gauge(id,mu)->c[k].im;
               }
             error(fwrite_LE_double(test,2*NG*NG,fp)!= 2*NG*NG,
                   1,"read_gauge_field_openQCD",
                   "Failed to read header from file");               
             
             
            }
  
  
  fclose(fp); 
  
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("IO",0,"Configuration [%s] wrote [%ld sec %ld usec]\n",filename,etime.tv_sec,etime.tv_usec);
  
}
