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

enum { BIG_ENDIAN, LITTLE_ENDIAN };
static int which_endian() {
  uint16_t one=1;
  char c=*((char*)(&one));
  if(c==0) lprintf("IO",100,"Big Endian on machine.\n");
  else lprintf("IO",100,"Little Endian on machine.\n");
  return (c==0)?BIG_ENDIAN:LITTLE_ENDIAN;
}

int fwrite_BE_int(int* ptr, size_t n, FILE* fp) {
  int ret;
  if(which_endian()==BIG_ENDIAN) {
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
  if(which_endian()==LITTLE_ENDIAN) {
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
  if(which_endian()==BIG_ENDIAN) {
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
  if(which_endian()==LITTLE_ENDIAN) {
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
  if(which_endian()==BIG_ENDIAN) {
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
  if(which_endian()==LITTLE_ENDIAN) {
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
  if(which_endian()==BIG_ENDIAN) {
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
  if(which_endian()==LITTLE_ENDIAN) {
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


void gaugefield_swapendian() {
  _DECLARE_INT_ITERATOR(ix);
  int a;
  double plaq;

  if(sizeof(double)*CHAR_BIT==32) {
    _MASTER_FOR(&glattice,ix) {
      for(a=0; a<2*NG*NG; a++) {
        swapendian32((double*)pu_gauge(ix,0)+a);
        swapendian32((double*)pu_gauge(ix,1)+a);
        swapendian32((double*)pu_gauge(ix,2)+a);
        swapendian32((double*)pu_gauge(ix,3)+a);
      }
    }
  } else if(sizeof(double)*CHAR_BIT==64) {
    _MASTER_FOR(&glattice,ix) {
      for(a=0; a<2*NG*NG; a++) {
        swapendian64((double*)pu_gauge(ix,0)+a);
        swapendian64((double*)pu_gauge(ix,1)+a);
        swapendian64((double*)pu_gauge(ix,2)+a);
        swapendian64((double*)pu_gauge(ix,3)+a);
      }
    }
  } else {
    error(1,1,"gaugefield_swapendian","sizeof(double)*CHAR_BIT != 32, 64 !\n");
  }

  /* start sendrecv of global gauge field */
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  plaq=avr_plaquette();
  lprintf("IO",0,"Swap endian.  New plaquette=%e\n",plaq);
}

void gaugefield_swapendian_auto() {
  double plaq;

  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  plaq=avr_plaquette();
  if(!isnan(plaq)) {
    lprintf("IO",0,"Endianness not swapped.\n");
    return;
  }
  
  gaugefield_swapendian();
  lprintf("IO",0,"Endianness swapped.\n");
  plaq=avr_plaquette();
  error(isnan(plaq),1,"gaugefield_swapendian_auto",
        "Endianness not determined.\n");
}


