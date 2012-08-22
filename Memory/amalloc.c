/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File amalloc.c
* 
* Align memory (de)allocation functions
*
*******************************************************************************/

#include <stdlib.h>
#include "memory.h"


#ifdef AMALLOC_MEASURE

#include <string.h>
#include "logger.h"

typedef struct {
  void* addr;
  size_t size;
} mem_map_t;

static mem_map_t *allocated = NULL;
static size_t nallocated = 0;
static double memory_Mb = 0.;
static double max_memory_Mb = 0.;

void insert(void *addr, size_t size) {
  mem_map_t *tmp;
  tmp = malloc(sizeof(mem_map_t)*(nallocated+1));
  
  if(nallocated>0) memcpy(tmp,allocated,sizeof(mem_map_t)*nallocated);
  tmp[nallocated].addr=addr;
  tmp[nallocated].size=size;
  nallocated++;
  
  if(nallocated>0) free(allocated);
  allocated=tmp;
  
  memory_Mb += size/1048576.;
  
  if(memory_Mb>max_memory_Mb) {
    max_memory_Mb=memory_Mb;
    lprintf("AMALLOC",0,"Total memory allocated with amalloc: %.6f Mb\n", max_memory_Mb);
  }

  /*lprintf("AMALLOC",0,"Total memory allocated with amalloc: %.6f Mb\n", memory_Mb);*/
}

void remove(void *addr) {
  error(nallocated==0,1,"amalloc.c","[remove] nallocated==0!");
  int found=(1==0);
  size_t at;
  mem_map_t *tmp;
  
  for(at=0;at<nallocated;at++)
  if(allocated[at].addr==addr) {
    found=(1==1);
    break;
  }
  
  error(!found,1,"amalloc.c","[remove] !found");

  memory_Mb -= allocated[at].size/1048576.;
  
  tmp=NULL;
  if(nallocated!=1) tmp = malloc(sizeof(mem_map_t)*(nallocated-1));
  if(at!=0) memcpy(tmp,allocated,sizeof(mem_map_t)*at);
  if(nallocated-at-1!=0) memcpy(tmp+at,allocated+at+1,sizeof(mem_map_t)*(nallocated-at-1));
  nallocated--;
  free(allocated);
  allocated=tmp;

  /*lprintf("AMALLOC",0,"Total memory allocated with amalloc: %.6f Mb\n", memory_Mb);*/
}

#endif



struct addr_t
{
   char *addr;
   char *true_addr;
   struct addr_t *next;
};

static struct addr_t *first=NULL;

void *amalloc(size_t size,int p)
{
   int shift;
   char *true_addr,*addr;
   unsigned long mask;
   struct addr_t *new;

   if ((size<=0)||(p<0))
      return(NULL);

   shift=1<<p;
   mask=(unsigned long)(shift-1);

   true_addr=malloc(size+shift);
   new=malloc(sizeof(*first));
   
   if ((true_addr==NULL)||(new==NULL))
   {
      free(true_addr);
      free(new);
      return(NULL);
   }

   addr=(char*)(((unsigned long)(true_addr+shift))&(~mask));
   (*new).addr=addr;
   (*new).true_addr=true_addr;
   (*new).next=first;
   first=new;

#ifdef AMALLOC_MEASURE
   insert((void*)addr,size);
#endif

   return((void*)(addr));
}


void afree(void *addr)
{
   struct addr_t *p,*q;

#ifdef AMALLOC_MEASURE
   remove(addr);
#endif

   q=NULL;
   
   for (p=first;p!=NULL;p=(*p).next)
   {
      if ((*p).addr==addr)
      {
         if (q!=NULL)
            (*q).next=(*p).next;
         else
            first=(*p).next;
         
         free((*p).true_addr);
         free(p);
         return;
      }
      
      q=p;
   }
}

