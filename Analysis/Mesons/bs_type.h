#ifndef __BS_TYPE
#define __BS_TYPE
#include <map>
#include <string>
#include "datasample.h"

class Corr_t {
public:
  int length;
  datasample* d;
  Corr_t(int lt) { length = lt/2+1; d = new datasample[length]; }
  ~Corr_t() { delete[] d; }
  void purge() { for(int i=0;i<length;i++) d[i].clear(); }
};


#endif
