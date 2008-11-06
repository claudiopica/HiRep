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
                   
   return((void*)(addr));
}


void afree(void *addr)
{
   struct addr_t *p,*q;

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

