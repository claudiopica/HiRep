#ifndef DATASAMPLE_H
#define DATASAMPLE_H

#include <vector>
#include <iostream>

using std::vector;
using std::ostream;

struct estimate {
  double val;
  double err;

  friend ostream &operator<<(ostream &out, estimate e){
    return out<<e.val<<" +- "<<e.err;
  }
};

class datasample: public vector<double> {
public:
  datasample(){}
  ~datasample(){}

  estimate avr() const;
  estimate stderr() const;
  estimate autocorr() const;
  void ci(double low, double med, double hi, double *il, double *im, double *ih) const;
  void hist(double *n, int npt, double *x, double w) const;
  void cumhist(double *n, int npt, double *x) const;

  estimate JKavr(const int blsize) const;
  estimate JKvar(const int blsize) const;

  friend estimate JKcov(const datasample &dt1, const datasample &dt2, const int blsize);

};

#endif //#ifndef DATASAMPLE_H
