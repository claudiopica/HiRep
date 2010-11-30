#include "progbar.h"
#include <iostream>

using std::flush;

const char CSI[]="[";
const char eraseline[]="[2K";
const char spin[]="(O)|"; //"/-\|", "(O)O", ".oOo", "(O)|" 

void DrawBar(double p, std::ostream &out) {
  const double w=79.;
  static unsigned int ss=0;
  int f=static_cast<int>(w);
  //out<<CSI<<"?25l";
  out<<CSI<<"0G";
  out<<eraseline<<'[';
  int n=static_cast<int>(p*(w-2.)/100.);
  for(int i=0; i<n; i++) out<<'#';
  out<<CSI<<f<<'G';
  out<<']'<<spin[(ss++)&3];
  int c=f/2-2;
  out<<CSI<<c<<'G';
  out<<CSI<<"7m";
  if(p>=100.) {
    out<<"DONE\n";
  } else {
    out<<static_cast<int>(p)<<'%';
  }
  out<<CSI<<"27m";
  out<<flush;
}


