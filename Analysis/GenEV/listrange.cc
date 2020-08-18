#include <list>
#include <string>
#include <iostream>
#include "listrange.h"
#include <sstream>
#include <algorithm> // std::sort
#include <cstdio>

oplist::oplist(par *apar)
{

  numop = apar->numop;
  Op_max = apar->Op_max;
  Op_min = apar->Op_min;
  oprange = apar->oprange;

  this->populate(apar->opstring);
};

op_bl *oplist::find_op(int op_index)
{
  list<op_bl>::iterator it1;
  int counter = 0;
  op_bl *ret = NULL;
  for (it1 = this->begin(); it1 != this->end(); it1++)
  {
    if (it1->op_index == op_index)
    {
      ret = &(*it1);
      counter++;
    }
  }
  if (counter <= 1)
    return ret;
  else
  {
    cerr << "[Error][find_op] Duplicate operator in the op_list." << endl;
    exit(1);
  }
};

bool oplist::add_op()
{
  op_bl supp;
  vector<int>::iterator it;

  if (this->back().is_set)
  {
    supp.good_bl.clear();
    supp.bad_bl.clear();
    supp.is_set = false;

    it = find(this->oprange.begin(), this->oprange.end(), this->back().op_index);
    if (*it < this->Op_max)
      supp.op_index = *(++it);
    else
      return false;

    for (int j = 1; j <= nblock; j++)
    {
      if (j >= Bl_min && j <= Bl_max)
        supp.good_bl.push_back(j);
      else
        supp.bad_bl.push_back(j);
    }
    this->push_back(supp);
  }
  else
  {

    int bl_change = this->back().good_bl.back();
    this->back().good_bl.pop_back();
    this->back().bad_bl.push_back(bl_change);
    std::sort(this->back().bad_bl.begin(), this->back().bad_bl.end());
    if (this->back().good_bl.empty())
      this->back().is_set = true;
  }

  return true;
};

static inline void help_report(int numop, int nblock)
{
  cerr << "Operators string not correct." << endl;
  cerr << "Should be written as n1-n2,...,ni-nj, with n1 <= n2 <= n3 <= " << numop << " ." << endl;
  cerr << "or n1:b11/b12,n2:b21/b31,n3-n4,with n1 <= n2 <= n3 <= " << numop << " and bi1,bi2 <= " << nblock << " ." << endl;
};

void oplist::report()
{
  string myrep;
  char tmp[256];

  oplist::iterator it1;
  vector<int>::iterator it2;

  myrep.clear();

  for (it1 = this->begin(); it1 != this->end(); it1++)
  {
    cout << "op : " << (it1->op_index) << " bl: ";
    for (it2 = it1->good_bl.begin(); it2 != it1->good_bl.end(); it2++)
      cout << *it2 << " ";
    cout << endl;
  }

  int oprange[2] = {-1, -1};
  int blrange[2] = {0, 0};
  for (it1 = this->begin(); it1 != this->end(); it1++)
  {
    if (it1->good_bl.empty())
      continue;
    int op = it1->op_index;
    if (oprange[0] == -1)
      oprange[0] = oprange[1] = op;
    bool fullop = true;
    std::sort(it1->good_bl.begin(), it1->good_bl.end());

    if ((int)(it1->good_bl.size()) != Bl_max - Bl_min + 1)
      fullop = false;

    //cout << "fullop " << op << " " << fullop << endl;

    if (fullop)
    {
      if (oprange[1] != op)
      {
        if (oprange[1] + 1 == op && blrange[1] == 0)
          oprange[1] = op;
        else
        {
          if (blrange[1] == 0)
          {
            sprintf(tmp, "%d-%d,", oprange[0], oprange[1]);
            myrep.append(tmp);
          }
          blrange[0] = blrange[1] = 0;
          oprange[0] = oprange[1] = op;
        }
      }
    }
    else
    {
      //      cout << "ci entro for op " << op << " " << myrep << endl;
      if (oprange[1] != op && blrange[1] == 0)
      {
        sprintf(tmp, "%d-%d,", oprange[0], oprange[1]);
        //else sprintf(tmp,"%d:%d/%d,",oprange[0],blrange[0],blrange[1]);
        myrep.append(tmp);
      }

      oprange[0] = oprange[1] = op;
      blrange[0] = blrange[1] = *(it1->good_bl.begin());

      it2 = it1->good_bl.begin();
      it2++;

      for (; it2 != it1->good_bl.end(); it2++)
      {
        if (blrange[1] + 1 == *it2)
          blrange[1] = *it2;
        else
        {
          sprintf(tmp, "%d:%d/%d,", oprange[0], blrange[0], blrange[1]);
          myrep.append(tmp);
          blrange[0] = blrange[1] = *it2;
        }
      }
      sprintf(tmp, "%d:%d/%d,", oprange[0], blrange[0], blrange[1]);
      myrep.append(tmp);
      //   cout << "ci esco " << myrep << endl;
    }
  }
  if (blrange[1] == 0)
  {
    sprintf(tmp, "%d-%d,", oprange[0], oprange[1]);
    myrep.append(tmp);
  }

  if (!myrep.empty())
    myrep.erase(myrep.end() - 1);

  cout << "[INFO][OPLIST->REPORT] Op string: " << myrep << endl;
};

void oplist::populate(char *input)
{
  string str2;
  size_t found = 0, minus, slash, colon, n1_start = 0, n1_end = 0, n2_start = 0, n2_end = 0;
  int op_start, op_end, bl_start, bl_end;
  bool control = true;
  op_bl supp, *tsupp;

  string str(input);
  this->clear();

  if (str.empty())
  {
    for (int i = 1; i < this->Op_min; i++)
    {
      supp.good_bl.clear();
      supp.bad_bl.clear();
      supp.op_index = i;
      supp.is_set = true;
      for (int j = 1; j <= nblock; j++)
        supp.bad_bl.push_back(j);

      this->push_back(supp);
    }
    supp.good_bl.clear();
    supp.bad_bl.clear();
    supp.op_index = Op_min;
    supp.is_set = false;
    for (int j = 1; j <= nblock; j++)
    {
      if (j >= Bl_min && j <= Bl_max)
        supp.good_bl.push_back(j);
      else
        supp.bad_bl.push_back(j);
    }
    this->push_back(supp);
    return;
  }

  while (control)
  {
    found = str.find(",", n1_start);
    n2_end = found;
    if (found == string::npos)
    {
      control = false;
      n2_end = str.length();
    }
    str2 = str.substr(n1_start, n2_end - n1_start);
    colon = str2.find(":", 0);

    if (colon == string::npos)
    {
      //cerr<<"--------------------------Non ci sono i due punti! qui->"<<str2<<endl;
      minus = str.find("-", n1_start);
      //cerr<<minus<<endl;
      if (minus == string::npos)
      {
        help_report(numop, nblock);
        exit(1);
      }
      //	  cerr<<"--------------------------Ho trovato un intervallo di operatori! qui->"<<str2<<endl;
      n2_start = minus + 1;
      n1_end = minus;
      str2 = str.substr(n1_start, n1_end - n1_start);
      //	  cerr<<"--------------------------Questo qui e\' il primo ->"<<str2<<endl;
      istringstream stream1(str2);
      stream1 >> op_start;
      str2 = str.substr(n2_start, n2_end - n2_start);
      istringstream stream2(str2);
      stream2 >> op_end;

      if (op_start > Op_max || op_start < Op_min || op_end > Op_max || op_end < Op_min || op_end < op_start)
      {
        help_report(numop, nblock);
        exit(1);
      }

      if (!this->empty())
        if (this->back().op_index >= op_start)
        {
          help_report(numop, nblock);
          exit(1);
        }

      for (int i = op_start; i <= op_end; i++)
      {
        supp.good_bl.clear();
        supp.bad_bl.clear();
        supp.op_index = i;
        supp.is_set = true;
        for (int j = 1; j <= nblock; j++)
        {
          if (j >= Bl_min && j <= Bl_max)
            supp.good_bl.push_back(j);
          else
            supp.bad_bl.push_back(j);
        }
        this->push_back(supp);
      }
    }
    else
    {
      colon = colon + n1_start;
      //	  cerr<<colon<<endl;
      //	  cerr<<"--------------------------Ci sono i due punti! qui->"<<str2<<endl;
      n2_start = colon + 1;
      n1_end = colon;
      str2 = str.substr(n1_start, n1_end - n1_start);
      // 	  cerr<<"--------------------------Questo qui e\' l\'operatore ->"<<str2<<endl;
      istringstream stream1(str2);
      stream1 >> op_start;

      //	  supp.op_index=op_start;
      //  supp.is_set=true;

      str2 = str.substr(n2_start, n2_end - n2_start);

      slash = str2.find("/", 0);

      if (slash == string::npos)
      {
        help_report(numop, nblock);
        exit(1);
      }

      n1_start = n2_start;
      n1_end = n1_start + slash;
      n2_start = n1_end + 1;

      str2 = str.substr(n1_start, n1_end - n1_start);
      istringstream stream2(str2);
      stream2 >> bl_start;

      str2 = str.substr(n2_start, n2_end - n2_start);
      istringstream stream3(str2);
      stream3 >> bl_end;

      if (op_start > Op_max || op_start < Op_min)
      {
        help_report(numop, nblock);
        exit(1);
      }

      if (!this->empty())
        if (this->back().op_index > op_start)
        {
          help_report(numop, nblock);
          exit(1);
        }

      tsupp = find_op(op_start);
      if (tsupp == NULL)
      {
        supp.good_bl.clear();
        supp.bad_bl.clear();
        supp.op_index = op_start;
        supp.is_set = true;
        for (int j = 1; j <= nblock; j++)
          supp.bad_bl.push_back(j);
        this->push_back(supp);
        tsupp = find_op(op_start);
      }

      vector<int>::iterator it;
      for (int j = bl_start; j <= bl_end; j++)
      {
        it = find(tsupp->good_bl.begin(), tsupp->good_bl.end(), j);
        if (it != tsupp->good_bl.end())
        {
          cerr << "[Error][oplist::populate] Blocking level already present in the positive list" << endl;
          exit(1);
        }
        tsupp->good_bl.push_back(j);

        it = find(tsupp->bad_bl.begin(), tsupp->bad_bl.end(), j);
        if (it == tsupp->bad_bl.end())
        {
          cerr << "[Error][oplist::populate] Blocking level not present in the negative list" << endl;
          exit(1);
        }
        tsupp->bad_bl.erase(it);
      }
      std::sort(tsupp->good_bl.begin(), tsupp->good_bl.end());
      std::sort(tsupp->bad_bl.begin(), tsupp->bad_bl.end());
    }
    n1_start = found + 1;
  }

  for (int i = this->begin()->op_index - 1; i > 0; i--)
  {
    supp.good_bl.clear();
    supp.bad_bl.clear();
    supp.op_index = i;
    supp.is_set = true;
    for (int j = 1; j <= nblock; j++)
      supp.bad_bl.push_back(j);

    this->push_front(supp);
  }

  /*Setting the start value*/
  this->back().is_set = false;

  // oplist::iterator it1;
  // vector<int>::iterator it2;
  // cout<< "good op"<< endl;
  // for (it1 = this->begin(); it1 != this->end(); it1++)
  //   for (it2 = it1->good_bl.begin(); it2 != it1->good_bl.end(); it2++)
  //     cout<<it1->op_index<<":"<<*it2<<endl;

  // cout<< "bad op"<< endl;
  // for (it1 = this->begin(); it1 != this->end(); it1++)
  //   for (it2 = it1->bad_bl.begin(); it2 != it1->bad_bl.end(); it2++)
  //     cout<<it1->op_index<<":"<<*it2<<endl;
}

int oplist::size()
{
  oplist::iterator it1;
  vector<int>::iterator it2;
  int size = 0;
  for (it1 = this->begin(); it1 != this->end(); it1++)
    for (it2 = it1->good_bl.begin(); it2 != it1->good_bl.end(); it2++)
      size++;
  return size;
}

double *oplist::cut(double *mat, int ntnbin)
{
  oplist::iterator it1;
  vector<int>::iterator it2;
  oplist::iterator it3;
  vector<int>::iterator it4;

  int size = this->size();

  int id = 0;
  double *ret = new double[size * size * ntnbin];
  for (int iout = 0; iout < ntnbin; iout++)
    for (it1 = this->begin(); it1 != this->end(); it1++)
      for (it2 = it1->good_bl.begin(); it2 != it1->good_bl.end(); it2++)
        for (it3 = this->begin(); it3 != this->end(); it3++)
          for (it4 = it3->good_bl.begin(); it4 != it3->good_bl.end(); it4++)
          {
            ret[id] = mat[(it1->op_index - 1) + (*it2 - 1) * numop + ((it3->op_index) - 1) * numop * nblock + (*it4 - 1) * numop * nblock * numop + iout * nblock * numop * nblock * numop];
            id++;
          }
  return ret;
}
