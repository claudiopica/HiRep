#include "bootstrap.h"
#include "ranlxs.h"


void b_sample(int len, int elen, int *bs){
  const int rl=128;
  static float r[rl];
  for (int i=0,k=elen;k>0;k-=rl) {
    int n=rl; if(k<rl) n=k;
    ranlxs(r,rl);
    for (;n>0;) {
      --n;
      bs[i++]=(int)(r[n]*float(len));
    }
  }
}

