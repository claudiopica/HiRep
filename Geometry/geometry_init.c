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

void geometry_init() {
	if (init) {
		int *cur;
		size_t req_mem=0;

		T=2;
		X=4;
		Y=2;
		Z=6;
		VOL3=X*Y*Z;
		VOLUME=VOL3*T;

		req_mem+=2*4*VOLUME; /* for iup and idn */
		req_mem+=VOLUME;     /* for ipt */
		req_mem+=VOLUME;     /* for ipt_4d */

		alloc_mem=malloc(req_mem*sizeof(int));
		error((alloc_mem==NULL),1,"geometry_init [geometry_init.c]",
         "Cannot allocate memory");

		cur=alloc_mem;
#define ALLOC(ptr,size) ptr=cur; cur+=(size) 

		/* iup and idn */
		ALLOC(iup,4*VOLUME);
		ALLOC(idn,4*VOLUME);
		/* ipt */
		ALLOC(ipt,VOLUME);
		/* ipt_4d */
		ALLOC(ipt_4d,VOLUME);


		atexit(&free_memory);

		init=0;
	}
}






void geometry_set()
{
  if (init) {
    T=4;
    X=4;
    Y=4;
    Z=4;
    VOL3=X*Y*Z;
    VOLUME=VOL3*T;
  }	
}

void geometry_mem_alloc(geometry_descriptor * gd) {
	if (init) {
		int *cur;
		size_t req_mem=0;
		unsigned int VOL_SIZE=gd->gsize;
		
		req_mem+=2*4*VOL_SIZE; /* for iup and idn */
		req_mem+=VOLUME;     /* for ipt */
		
		alloc_mem=malloc(req_mem*sizeof(int));
		error((alloc_mem==NULL),1,"geometry_init [geometry_init.c]",
		      "Cannot allocate memory");
		
		cur=alloc_mem;
		
		/* iup and idn */
		ALLOC(iup,4*VOL_SIZE);
		ALLOC(idn,4*VOL_SIZE);
		/* ipt */
		ALLOC(ipt,VOLUME);
		/* ipt_4d */
		

		atexit(&free_memory);

		init=0;
	}
}



#undef ALLOC




