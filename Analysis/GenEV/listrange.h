#ifndef _OPLIST
#define _OPLIST
#include <list>
#include "type.h"

using namespace std;

class oplist : public list<op_bl> {
 public:
  oplist(par * apar);
  
  bool add_op();
  void populate(char * str);
  op_bl * find_op(int op_index);
  void report();
  double * cut(double * mat,int ntnbin);
  int size();

 private:
  int Bl_max;
  int Bl_min;
  int nblock;
  int Op_max;
  int Op_min;
  int numop;

  vector<int> oprange;
};

#endif
