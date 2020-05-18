#ifndef _TYPE
#define _TYPE

#include <vector>
#include <string>
#include <map>

struct datafile
{
  std::string filename;
  int nmeas;
};

class cpoints
{
public:
  cpoints(int i1, int i2)
  {
    p1 = i1;
    p2 = i2;
  };

  int p1;
  int p2;
};

struct dtcor
{
  int index;
  std::vector<cpoints> points;
};

struct par
{
  par() : iflag(0), vflag(0), Dflag(0), rflag(0), Gflag(0), fflag(0), Jflag(0), xflag(0), dflag(0), eflag(0), Tflag(0), numop(0), ndt(0), nt(0), binwidth(1), numbinjk(0), lt(100), n_states(6){};

  char list_dir[256];
  char cor_def_filename[256];
  char vev_name[100];
  char opstring[512];
  int iflag, vflag, Dflag, rflag, Gflag, fflag, Jflag, xflag, dflag, eflag, Tflag;
  int numop;
  int ndt;
  int nt;
  int binwidth, numbinjk;
  int p_diag, p_inv;
  int lt;
  int n_states;
  std::map<int, dtcor> corrdef;
  std::map<int, int> t2idx;
  std::vector<int> activeop;
  std::vector<datafile> vevfiles;
};

#endif
