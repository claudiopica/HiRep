#ifndef __BS_CTRL
#define __BS_CTRL

#include <iostream>
#include <set>
#include <string>
#include <map>
#include "bs_type.h"


#define MAX_NLEVELS 9

class eval_ctrl {
 protected:
  bool level[MAX_NLEVELS];
  std::string name;
    
 public:
  virtual ~eval_ctrl() {};
  
  virtual void prop_dep()=0;

  virtual void prop_cuts() {};
  
  virtual void allocate_datamemory(int lt) {};
  
  virtual void purge_b2() {};
  
  virtual void print()=0;
  
  virtual bool isprimary() {return false;}
  virtual bool isderived() {return false;}
  virtual bool isratio() {return false;}

  static std::set<eval_ctrl*> all;
  static bool bs_2nd;
  static std::map<std::string, int> left_cut;
  static std::map<std::string, int> right_cut;
  static void normalize_cuts(int lt);
  static bool fill_dep(const char* arg);
  static std::map<std::string, Corr_t*> channel_map;

  virtual void set_level(int lv)=0;
  bool get_level(int i){return level[i];}
  int get_left_cut();
  int get_right_cut();

  virtual std::ostream& print_utility(std::ostream& os)=0;
  friend std::ostream& operator<<(std::ostream& os, eval_ctrl& ev); 
};



enum { P_READ=0, P_COR_B1, P_EFF_B1, P_COR_B2, P_EFF_B2, P_FIT, P_AVE_COR, P_AVE_EFF, P_AVE_FIT };
enum { D_EFF_B1=0, D_EFF_B2, D_FIT, D_AVE_EFF, D_AVE_FIT };
enum { R_ACTIVE=0 };

class primary_ctrl : public eval_ctrl {
 private:
  int chan_counter;
 
  void prop_dep();

  std::ostream& print_utility(std::ostream& os);  
  
 public:
  int nchan;
  char** channel;
  
  Corr_t** cor_b0;
  Corr_t* cor_b1;
  Corr_t* cor_b2;
  Corr_t* eff_b1;
  Corr_t* eff_b2;
  datasample fit;

  double** ws;

  static std::set<primary_ctrl*> all_pr;

  primary_ctrl(const char* n, int nc);
  ~primary_ctrl();
  void add_channel(const char * c1);
  void set_level(int lv);
  void allocate_datamemory(int lt);
  void purge_b2();
  void print();
  bool isprimary() {return true;};
  
  bool eval_cor_eff(const int bsl, const int* bs, const int len, int Lt, int effmass_method);
  void store_cor_eff(int bs);
};


class derived_ctrl : public eval_ctrl {
 private:
  void prop_dep();
  void prop_cuts();

  std::ostream& print_utility(std::ostream& os);
  
 public:
  int ndep[3];
  eval_ctrl** dep[3];
  int* dep_level[3];
  
  Corr_t* eff_b1;
  Corr_t* eff_b2;
  datasample fit;

  static std::set<derived_ctrl*> all_der;

  derived_ctrl(const char* n);
  ~derived_ctrl();
  void set_ndeps(int lv, int nd);
  void set_level(int lv);
  void allocate_datamemory(int lt);
  void purge_b2();
  void print();
  bool isderived() {return true;};
};


class ratio_ctrl : public eval_ctrl {
 private:
  void prop_dep();

  std::ostream& print_utility(std::ostream& os);  
  
 public:
  int ndep;
  eval_ctrl** dep;
  int* dep_level;
  
  datasample fit;

  ratio_ctrl(const char* n);
  ~ratio_ctrl();
  void set_ndeps(int nd);
  void set_level(int lv);
  void allocate_datamemory(int lt);
  void print();
  bool isratio() {return true;};
};

#endif
