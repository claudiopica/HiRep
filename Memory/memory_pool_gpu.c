/***************************************************************************\
* Copyright (c) 2012 Ari Hietanen                                           *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifdef WITH_GPU
#include "gpu.h"

typedef struct _memory_pool_elem memory_pool_elem;

struct _memory_pool_elem{
  void* mem;
  size_t size;
  memory_pool_elem* next;
};

static memory_pool_elem* memory_pool =NULL;
static memory_pool_elem* memory_allocated=NULL;

cudaError_t alloc_pool_gpu(void** gpu_ptr, size_t size){
  memory_pool_elem*  mem_pool_it1 = memory_pool;
  memory_pool_elem*  mem_pool_it2;
  cudaError_t err;
  while (mem_pool_it1 != NULL){
    if (mem_pool_it1->size == size){
      if (mem_pool_it1 == memory_pool){
	memory_pool = memory_pool->next;
      }
      else{
	mem_pool_it2->next = mem_pool_it1->next;
      }
      /*Move to allocated memory*/
      mem_pool_it1->next = memory_allocated;
      memory_allocated = mem_pool_it1;
      *gpu_ptr =  mem_pool_it1->mem;
      return cudaSuccess;
    }
    mem_pool_it2 = mem_pool_it1; 
    mem_pool_it1 = mem_pool_it1->next; 
  }
  err=cudaMalloc((void **) gpu_ptr,size);
  mem_pool_it1 = (memory_pool_elem*) malloc(sizeof(memory_pool_elem));
  mem_pool_it1->next = memory_allocated;
  mem_pool_it1->size = size ;
  mem_pool_it1->mem = *gpu_ptr;
  memory_allocated = mem_pool_it1;
  return err;
}

void free_pool_gpu(void* gpu_ptr){
  memory_pool_elem*  mem_pool_it1 = memory_allocated;
  memory_pool_elem*  mem_pool_it2;
  while (mem_pool_it1 != NULL){
    if (mem_pool_it1->mem == gpu_ptr){
      if (mem_pool_it1 == memory_allocated){
	memory_allocated = memory_allocated->next;
      }
      else{
	mem_pool_it2->next = mem_pool_it1->next;
      }
      //Move to pool
      mem_pool_it1->next = memory_pool;
      memory_pool = mem_pool_it1;
      return;
    }
    mem_pool_it2 = mem_pool_it1; 
    mem_pool_it1 = mem_pool_it1->next;     
  }
  error(1,1,"gpu_free_pool [memory_pool_gpu.c]","Memory pool corrupted");
}

void erase_pool_gpu(){
  memory_pool_elem*  mem_pool_it1 = memory_pool;
  memory_pool_elem*  mem_pool_it2;  
  while(mem_pool_it1 != NULL){
    mem_pool_it2=mem_pool_it1;
    mem_pool_it1 = mem_pool_it1->next;
    cudaFree(mem_pool_it2->mem);
    free(mem_pool_it2);
  }
}

#endif
