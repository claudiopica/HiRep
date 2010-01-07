#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include "observables.h"
#include "communications.h"
#include "memory.h"
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <math.h>



void swapendian16(void* ptr) {
  uint16_t x;
  memcpy(&x,ptr,16/CHAR_BIT);
  x = (x>>8) | 
      (x<<8);
  memcpy(ptr,&x,16/CHAR_BIT);
}

void swapendian32(void* ptr) {
  uint32_t x;
  memcpy(&x,ptr,32/CHAR_BIT);
  x = (x>>24) | 
      ((x<<8) & 0x00FF0000) |
      ((x>>8) & 0x0000FF00) |
      (x<<24);
  memcpy(ptr,&x,32/CHAR_BIT);
}

void swapendian64(void* ptr) {
  uint64_t x;
  memcpy(&x,ptr,64/CHAR_BIT);
  x = (x>>56) | 
      ((x<<40) & 0x00FF000000000000) |
      ((x<<24) & 0x0000FF0000000000) |
      ((x<<8)  & 0x000000FF00000000) |
      ((x>>8)  & 0x00000000FF000000) |
      ((x>>24) & 0x0000000000FF0000) |
      ((x>>40) & 0x000000000000FF00) |
      (x<<56);
  memcpy(ptr,&x,64/CHAR_BIT);
}

enum { _BIG_ENDIAN_HRP, _LITTLE_ENDIAN_HRP };

static int which_endian() {
  uint16_t one=1;
  char c=*((char*)(&one));
  if(c==0) lprintf("IO",100,"Big Endian on machine.\n");
  else lprintf("IO",100,"Little Endian on machine.\n");
  return (c==0)?_BIG_ENDIAN_HRP:_LITTLE_ENDIAN_HRP;
}

int fwrite_BE_int(int* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==_BIG_ENDIAN_HRP) {
    ret=fwrite(ptr,sizeof(int),n,fp);
    lprintf("IO",100,"Written %d integers.\n",n);
  } else {
    int* ptr_l=malloc(sizeof(int)*n);
    memcpy(ptr_l,ptr,sizeof(int)*n);
    int i;
    if(sizeof(int)*CHAR_BIT==16)
      for(i=0; i<n; i++) swapendian16(ptr_l+i);
    else if(sizeof(int)*CHAR_BIT==32)
      for(i=0; i<n; i++) swapendian32(ptr_l+i);
    else if(sizeof(int)*CHAR_BIT==64)
      for(i=0; i<n; i++) swapendian64(ptr_l+i);
    ret=fwrite(ptr_l,sizeof(int),n,fp);
    lprintf("IO",100,"Written %d integers with swapped endian.\n",n);
    free(ptr_l);
  }
  return ret;
}

int fwrite_LE_int(int* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==_LITTLE_ENDIAN_HRP) {
    ret=fwrite(ptr,sizeof(int),n,fp);
    lprintf("IO",100,"Written %d integers.\n",n);
  } else {
    int* ptr_l=malloc(sizeof(int)*n);
    memcpy(ptr_l,ptr,sizeof(int)*n);
    int i;
    if(sizeof(int)*CHAR_BIT==16)
      for(i=0; i<n; i++) swapendian16(ptr_l+i);
    else if(sizeof(int)*CHAR_BIT==32)
      for(i=0; i<n; i++) swapendian32(ptr_l+i);
    else if(sizeof(int)*CHAR_BIT==64)
      for(i=0; i<n; i++) swapendian64(ptr_l+i);
    ret=fwrite(ptr_l,sizeof(int),n,fp);
    lprintf("IO",100,"Written %d integers with swapped endian.\n",n);
    free(ptr_l);
  }
  return ret;
}

int fwrite_BE_double(double* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==_BIG_ENDIAN_HRP) {
    ret=fwrite(ptr,sizeof(double),n,fp);
    lprintf("IO",100,"Written %d doubles.\n",n);
  } else {
    double* ptr_l=malloc(sizeof(double)*n);
    memcpy(ptr_l,ptr,sizeof(double)*n);
    int i;
    if(sizeof(double)*CHAR_BIT==16)
      for(i=0; i<n; i++) swapendian16(ptr_l+i);
    else if(sizeof(double)*CHAR_BIT==32)
      for(i=0; i<n; i++) swapendian32(ptr_l+i);
    else if(sizeof(double)*CHAR_BIT==64)
      for(i=0; i<n; i++) swapendian64(ptr_l+i);
    ret=fwrite(ptr_l,sizeof(double),n,fp);
    lprintf("IO",100,"Written %d doubles with swapped endian.\n",n);
    free(ptr_l);
  }
  return ret;
}

int fwrite_LE_double(double* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==_LITTLE_ENDIAN_HRP) {
    ret=fwrite(ptr,sizeof(double),n,fp);
    lprintf("IO",100,"Written %d doubles.\n",n);
  } else {
    double* ptr_l=malloc(sizeof(double)*n);
    memcpy(ptr_l,ptr,sizeof(double)*n);
    int i;
    if(sizeof(double)*CHAR_BIT==16)
      for(i=0; i<n; i++) swapendian16(ptr_l+i);
    else if(sizeof(double)*CHAR_BIT==32)
      for(i=0; i<n; i++) swapendian32(ptr_l+i);
    else if(sizeof(double)*CHAR_BIT==64)
      for(i=0; i<n; i++) swapendian64(ptr_l+i);
    ret=fwrite(ptr_l,sizeof(double),n,fp);
    lprintf("IO",100,"Written %d doubles with swapped endian.\n",n);
    free(ptr_l);
  }
  return ret;
}

int fread_BE_int(int* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==_BIG_ENDIAN_HRP) {
    ret=fread(ptr,sizeof(int),n,fp);
    lprintf("IO",100,"Read %d integers.\n",n);
  } else {
    int i;
    ret=fread(ptr,sizeof(int),n,fp);
    if(sizeof(int)*CHAR_BIT==16)
      for(i=0; i<n; i++) swapendian16(ptr+i);
    else if(sizeof(int)*CHAR_BIT==32)
      for(i=0; i<n; i++) swapendian32(ptr+i);
    else if(sizeof(int)*CHAR_BIT==64)
      for(i=0; i<n; i++) swapendian64(ptr+i);
    lprintf("IO",100,"Read %d integers with swapped endian.\n",n);
  }
  return ret;
}

int fread_LE_int(int* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==_LITTLE_ENDIAN_HRP) {
    ret=fread(ptr,sizeof(int),n,fp);
    lprintf("IO",100,"Read %d integers.\n",n);
  } else {
    int i;
    ret=fread(ptr,sizeof(int),n,fp);
    if(sizeof(int)*CHAR_BIT==16)
      for(i=0; i<n; i++) swapendian16(ptr+i);
    else if(sizeof(int)*CHAR_BIT==32)
      for(i=0; i<n; i++) swapendian32(ptr+i);
    else if(sizeof(int)*CHAR_BIT==64)
      for(i=0; i<n; i++) swapendian64(ptr+i);
    lprintf("IO",100,"Read %d integers with swapped endian.\n",n);
  }
  return ret;
}

int fread_BE_double(double* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==_BIG_ENDIAN_HRP) {
    ret=fread(ptr,sizeof(double),n,fp);
    lprintf("IO",100,"Read %d doubles.\n",n);
  } else {
    int i;
    ret=fread(ptr,sizeof(double),n,fp);
    if(sizeof(double)*CHAR_BIT==16)
      for(i=0; i<n; i++) swapendian16(ptr+i);
    else if(sizeof(double)*CHAR_BIT==32)
      for(i=0; i<n; i++) swapendian32(ptr+i);
    else if(sizeof(double)*CHAR_BIT==64)
      for(i=0; i<n; i++) swapendian64(ptr+i);
    lprintf("IO",100,"Read %d doubles with swapped endian.\n",n);
  }
  return ret;
}

int fread_LE_double(double* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==_LITTLE_ENDIAN_HRP) {
    ret=fread(ptr,sizeof(double),n,fp);
    lprintf("IO",100,"Read %d doubles.\n",n);
  } else {
    int i;
    ret=fread(ptr,sizeof(double),n,fp);
    if(sizeof(double)*CHAR_BIT==16)
      for(i=0; i<n; i++) swapendian16(ptr+i);
    else if(sizeof(double)*CHAR_BIT==32)
      for(i=0; i<n; i++) swapendian32(ptr+i);
    else if(sizeof(double)*CHAR_BIT==64)
      for(i=0; i<n; i++) swapendian64(ptr+i);
    lprintf("IO",100,"Read %d doubles with swapped endian.\n",n);
  }
  return ret;
}

