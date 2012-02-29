#include <stdlib.h>
#include "suN.h"
#include "global.h"
#include "logger.h"

#ifdef TWISTED_BC

#ifdef UPDATE_EO
#error TWISTED_BC must be used without UPDATE_EO
#endif

void init_twbc(){

  int g[4],i,j,k;
  int coord[3];
  int counter,counterplus;

  lprintf("MAIN",0,"TWISTED boundary conditions\n");

  error((twbc_staples!=NULL),1,"init_twbc [init_twbc.c]",
	"Double initialization of the structure twbc_staples");
  
  twbc_staples=calloc(sizeof(int*),glattice.gsize*4);
/*   for(i=0;i<4*glattice.gsize;i++) twbc_staples[i]=malloc(sizeof(int)*6); */

  error((twbc_plaq!=NULL),1,"init_twbc [init_twbc.c]",
	"Double initialization of the structure twbc_plaq");
  
  twbc_plaq=malloc(sizeof(int)*glattice.gsize*16);
  for(i=0;i<16*glattice.gsize;i++) twbc_plaq[i]=1;
  
  for(g[0]=0;g[0]<T;g[0]++)
    for(g[1]=0;g[1]<X;g[1]++)
      for(g[2]=0;g[2]<Y;g[2]++)
	for(g[3]=0;g[3]<Z;g[3]++){
	  i=ipt(g[0],g[1],g[2],g[3]);
	  counter=0;
	  counterplus=0;

	  if(COORD[1]*X+g[1]==0 ) {coord[counter]=-1; counter++; }
	  else if (COORD[1]*X+g[1]==GLB_X-1) {coord[counter]=1; counter++; counterplus++; }
	  if(COORD[2]*Y+g[2]==0 ) {coord[counter]=-2; counter++; }
	  else if (COORD[2]*Y+g[2]==GLB_Y-1) {coord[counter]=2; counter++; counterplus++; }
	  if(COORD[3]*Z+g[3]==0 ) {coord[counter]=-3; counter++; }
	  else if (COORD[3]*Z+g[3]==GLB_Z-1) {coord[counter]=3; counter++; counterplus++; }

	  if(counter>1 && counterplus>0){
	    for(j=0;j<counter;j++) {
	      if(coord[j]>0){
		error(twbc_staples[i*4+coord[j]]!=NULL,1,"init_twbc","Already allocated\n");
 		twbc_staples[i*4+coord[j]]=malloc(sizeof(int)*6);
 		for(k=0;k<6;k++) twbc_staples[i*4+coord[j]][k]=1;
		for(k=0;k<counter;k++){
		  if(k==j) continue;
		  if(coord[k]>0) twbc_staples[i*4+coord[j]][(coord[k]-coord[j]+4)%4-1]=-1;
		  if(coord[k]<0) twbc_staples[i*4+coord[j]][(-coord[k]-coord[j]+8)%4+2]=-1;
		}
	      }
	    }
	  }


	  if(counterplus>1)
	    for(j=0;j<counter;j++)
	      if(coord[j]>0)
		for(k=0;k<counter;k++)
		  if(coord[k]>0){
		    if(k==j) continue;
		    twbc_plaq[i*16+coord[j]*4+coord[k]]=-1;
		  }
	  
	  
	}
  

}


void free_twbc(){
  int i;
  
  for(i=0;i<glattice.gsize*4;i++)
    if(twbc_staples[i]!=NULL) free(twbc_staples[i]);
  free(twbc_staples);
  
  free(twbc_plaq);
}
#endif
