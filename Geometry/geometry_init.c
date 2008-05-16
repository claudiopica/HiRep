/*******************************************************************************
*
* File geometry_init.c
*
* Inizialization of geometry structures
*
*******************************************************************************/

#include "geometry.h"
#include "global.h"
#include "error.h"
#include <stdlib.h>

static int init=1;
static int *alloc_mem=NULL;

static void free_memory() {
	if(alloc_mem!=NULL) {
		free(alloc_mem);
		alloc_mem=NULL;
		iup=idn=NULL;
		ipt=NULL;
		ipt_4d=NULL;
	}
}




void geometry_set()
{
  if (init) {
    T=5;
    X=7;
    Y=7;
    Z=8;
    T_BORDER=1;
    X_BORDER=1;  
    Y_BORDER=1;
    Z_BORDER=1;
    VOL3=X*Y*Z;
    VOLUME=VOL3*T;
    GLOBAL_T=np_t*T;
    GLOBAL_X=np_x*X;
    GLOBAL_Y=np_y*Y;
    GLOBAL_Z=np_z*Z;
  }	
}

void geometry_mem_alloc() {
	if (init) {
		int *cur;
		size_t req_mem=0;
		unsigned int VOL_SIZE=glattice.gsize;
		
		req_mem+=2*4*VOL_SIZE; /* for iup and idn */
		req_mem+=(X+2*X_BORDER)*(Y+2*Y_BORDER)*(Z+2*Z_BORDER)*(T+2*T_BORDER);     /* for ipt */
		
		alloc_mem=malloc(req_mem*sizeof(int));
		error((alloc_mem==NULL),1,"geometry_init [geometry_init.c]",
		      "Cannot allocate memory");
		
		cur=alloc_mem;
#define ALLOC(ptr,size) ptr=cur; cur+=(size) 
		
		/* iup and idn */
		ALLOC(iup,4*VOL_SIZE);
		ALLOC(idn,4*VOL_SIZE);
		/* ipt */

		ALLOC(ipt,(X+2*X_BORDER)*(Y+2*Y_BORDER)*(Z+2*Z_BORDER)*(T+2*T_BORDER));
		/* ipt_4d */
		

		atexit(&free_memory);

		init=0;
#undef ALLOC

	}
}







