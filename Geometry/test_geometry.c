/*******************************************************************************
*
* File geometry.c
*
* Definition of the lattice geometry
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "geometry.h"
#include "global.h"
#include "safe_mod.h"
#include "logger.h"

void test_geometry(void)
{
  int x0,x1,x2,x3,ix,mu;
  int check[VOLUME] ={0};
  int ckdir[4];

  for (x0=0;x0<VOLUME;x0++){
    check[x0]=0;
  }  

  /* TEST safe_mod first */
  if (safe_mod(-1,T)!=(T-1)) {
    lprintf("TESTING",50,"Error in safe_mod (test0)!\n");
    return;
  }
  if (safe_mod(-1,L)!=(L-1)) {
    lprintf("TESTING",50,"Error in safe_mod (test1)!\n");
    return;
  }
  if (safe_mod(T,T)!=0) {
    lprintf("TESTING",50,"Error in safe_mod (test2)!\n");
    return;
  }
  if (safe_mod(L,L)!=0) {
    lprintf("TESTING",50,"Error in safe_mod (test3)!\n");
    return;
  }


	lprintf("TESTING",50,"Checking geometry...");
  for (x0=0;x0<T;x0++){
    for (x1=0;x1<L;x1++){
      for (x2=0;x2<L;x2++){
	for (x3=0;x3<L;x3++){
	  ix = ipt[x0][x1][x2][x3];
	  check[ix]++;
	  ckdir[0] = ipt[safe_mod(x0+1,T)][x1][x2][x3];
	  ckdir[1] = ipt[x0][safe_mod(x1+1,L)][x2][x3];
	  ckdir[2] = ipt[x0][x1][safe_mod(x2+1,L)][x3];
	  ckdir[3] = ipt[x0][x1][x2][safe_mod(x3+1,L)];
	  for (mu=0;mu<4;++mu){
	    if (ckdir[mu]!=iup[ix][mu]) {
	      lprintf("TESTING",50," FAILED. [site %d=(%d,%d,%d,%d) dir %d up]\n",ix,x0,x1,x2,x3,mu);
	      return;
	    }
	  }
	  ckdir[0] = ipt[safe_mod(x0-1,T)][x1][x2][x3];
	  ckdir[1] = ipt[x0][safe_mod(x1-1,L)][x2][x3];
	  ckdir[2] = ipt[x0][x1][safe_mod(x2-1,L)][x3];
	  ckdir[3] = ipt[x0][x1][x2][safe_mod(x3-1,L)];
	  for (mu=0;mu<4;++mu){
	    if (ckdir[mu]!=idn[ix][mu]) {
	      lprintf("TESTING",50," FAILED. [site %d=(%d,%d,%d,%d) dir %d dn]\n",ix,x0,x1,x2,x3,mu);
	      return;
	    }
	  }
	}
      }
    }
  }

  for (x0=0;x0<VOLUME;x0++){
    if(check[x0]!=1) {
      lprintf("TESTING",50," FAILED. [site %d counted %d times]\n",x0, check[x0]);
      return;
    }
  }

  lprintf("TESTING",50," PASSED.\n");
}

static void find_coord(int ix, int *x0r, int *x1r, int *x2r, int *x3r)
{
  int found = 0;
  int x0, x1, x2, x3;

  for (x0=0;x0<T;x0++){
    for (x1=0;x1<L;x1++){
      for (x2=0;x2<L;x2++){
	for (x3=0;x3<L;x3++){
	  if (ix == ipt[x0][x1][x2][x3]){
	    if (found) {
	      lprintf("TESTING",50,"Geometry: errore in find_coord sul sito %d\n", ix);
	    } else {
	      ++found;
	      *x0r = x0; 
	      *x1r = x1; 
	      *x2r = x2; 
	      *x3r = x3; 
	    }
	  }
	}
      }
    }
  }  
}

void print_geometry(void)
{
  int x0,x1,x2,x3,ix,n;
  int x0r,x1r,x2r,x3r;

  for (x0=0;x0<T;x0++){
    for (x1=0;x1<L;x1++){
      for (x2=0;x2<L;x2++){
	for (x3=0;x3<L;x3++){
	  ix = ipt[x0][x1][x2][x3];
	  printf("sito %d=(%d,%d,%d,%d) ",ix, x0, x1, x2, x3);
	  n =iup[ix][0];
	  find_coord(n, &x0r, &x1r, &x2r, &x3r);
	  if ((safe_mod(x0+1,T)!=x0r)||(x1!=x1r)||(x2!=x2r)||(x3!=x3r))
	    printf("0up (%d,%d,%d,%d) ",x0r, x1r, x2r, x3r);
	  n =idn[ix][0];
	  find_coord(n, &x0r, &x1r, &x2r, &x3r);
	  if ((safe_mod(x0-1,T)!=x0r)||(x1!=x1r)||(x2!=x2r)||(x3!=x3r))
	    printf("0dn (%d,%d,%d,%d) ",x0r, x1r, x2r, x3r);
	  n =iup[ix][1];
	  find_coord(n, &x0r, &x1r, &x2r, &x3r);
	  if ((x0!=x0r)||(safe_mod(x1+1,L)!=x1r)||(x2!=x2r)||(x3!=x3r))
	    printf("1up (%d,%d,%d,%d) ",x0r, x1r, x2r, x3r);
	  n =idn[ix][1];
	  find_coord(n, &x0r, &x1r, &x2r, &x3r);
	  if ((x0!=x0r)||(safe_mod(x1-1,L)!=x1r)||(x2!=x2r)||(x3!=x3r))
	    printf("1dn (%d,%d,%d,%d) ",x0r, x1r, x2r, x3r);
	  n =iup[ix][2];
	  find_coord(n, &x0r, &x1r, &x2r, &x3r);
	  if ((x0!=x0r)||(x1!=x1r)||(safe_mod(x2+1,L)!=x2r)||(x3!=x3r))
	    printf("2up (%d,%d,%d,%d) ",x0r, x1r, x2r, x3r);
	  n =idn[ix][2];
	  find_coord(n, &x0r, &x1r, &x2r, &x3r);
	  if ((x0!=x0r)||(x1!=x1r)||(safe_mod(x2-1,L)!=x2r)||(x3!=x3r))
	    printf("2dn (%d,%d,%d,%d) ",x0r, x1r, x2r, x3r);
	  n =iup[ix][3];
	  find_coord(n, &x0r, &x1r, &x2r, &x3r);
	  if ((x0!=x0r)||(x1!=x1r)||(x2!=x2r)||(safe_mod(x3+1,L)!=x3r))
	    printf("3up (%d,%d,%d,%d) ",x0r, x1r, x2r, x3r);
	  n =idn[ix][3];
	  find_coord(n, &x0r, &x1r, &x2r, &x3r);
	  if ((x0!=x0r)||(x1!=x1r)||(x2!=x2r)||(safe_mod(x3-1,L)!=x3r))
	    printf("3dn (%d,%d,%d,%d)",x0r, x1r, x2r, x3r);
	  printf("\n");
	}
      }
    }
  }

}
