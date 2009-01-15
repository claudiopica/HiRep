#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include "observables.h"
#include "communications.h"
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <math.h>


void swapendian32(double* fl) {
  uint32_t x;
  memcpy(&x,fl,sizeof(double));
  x = (x>>24) | 
      ((x<<8) & 0x00FF0000) |
      ((x>>8) & 0x0000FF00) |
      (x<<24);
  memcpy(fl,&x,sizeof(double));
}

void swapendian64(double* fl) {
  uint64_t x;
  memcpy(&x,fl,sizeof(double));
  x = (x>>56) | 
      ((x<<40) & 0x00FF000000000000) |
      ((x<<24) & 0x0000FF0000000000) |
      ((x<<8)  & 0x000000FF00000000) |
      ((x>>8)  & 0x00000000FF000000) |
      ((x>>24) & 0x0000000000FF0000) |
      ((x>>40) & 0x000000000000FF00) |
      (x<<56);
  memcpy(fl,&x,sizeof(double));
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


