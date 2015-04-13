#include <iostream>
#include <string.h>
#include <stdlib.h>
#include "bs_ctrl.h"
#include "effective_mass.h"


using namespace std;


extern int blsize;

set<eval_ctrl*> eval_ctrl::all;
set<primary_ctrl*> primary_ctrl::all_pr;
set<derived_ctrl*> derived_ctrl::all_der;
bool eval_ctrl::bs_2nd=false;
std::map<std::string, int> eval_ctrl::left_cut;
std::map<std::string, int> eval_ctrl::right_cut;
std::map<std::string, Corr_t*> eval_ctrl::channel_map;

  
primary_ctrl id("id",1);
primary_ctrl g0("g0",1);
  
primary_ctrl g5("g5",1);
primary_ctrl g0g5("g0g5",1);

primary_ctrl g1("g1",1);
primary_ctrl g2("g2",1);
primary_ctrl g3("g3",1);
primary_ctrl gk("gk",3);
primary_ctrl g0g1("g0g1",1);
primary_ctrl g0g2("g0g2",1);
primary_ctrl g0g3("g0g3",1);
primary_ctrl g0gk("g0gk",3);

primary_ctrl g5g1("g5g1",1);
primary_ctrl g5g2("g5g2",1);
primary_ctrl g5g3("g5g3",1);
primary_ctrl g5gk("g5gk",3);

primary_ctrl g0g5g1("g0g5g1",1);
primary_ctrl g0g5g2("g0g5g2",1);
primary_ctrl g0g5g3("g0g5g3",1);
primary_ctrl g0g5gk("g0g5gk",3);

primary_ctrl g5_g0g5("g5_g0g5",1);


derived_ctrl mpcac("mpcac");
derived_ctrl gps("gps");
derived_ctrl fps("fps");
derived_ctrl fv("fv");
derived_ctrl fvk("fvk");
derived_ctrl fak("fak");


ratio_ctrl mvmps("mvmps");
ratio_ctrl mvkmps("mvkmps");
ratio_ctrl mpsfps("mpsfps");
ratio_ctrl mvfps("mvfps");
ratio_ctrl mvkfps("mvkfps");
ratio_ctrl gmor("gmor");
ratio_ctrl mpsmpcac("mpsmpcac");
ratio_ctrl mps2mpcac("mps2mpcac");
ratio_ctrl fvfps("fvfps");
ratio_ctrl fvkfps("fvkfps");
ratio_ctrl mamps("mamps");
ratio_ctrl makmps("makmps");
ratio_ctrl mamv("mamv");
ratio_ctrl makmvk("makmvk");


bool eval_ctrl::fill_dep(const char* arg){

  id.add_channel("id");
  g0.add_channel("g0");
  
  g5.add_channel("g5");
  g0g5.add_channel("g0g5");

  g1.add_channel("g1");
  g2.add_channel("g2");
  g3.add_channel("g3");
  gk.add_channel("g1");
  gk.add_channel("g2");
  gk.add_channel("g3");

  g0g1.add_channel("g0g1");
  g0g2.add_channel("g0g2");
  g0g3.add_channel("g0g3");
  g0gk.add_channel("g0g1");
  g0gk.add_channel("g0g2");
  g0gk.add_channel("g0g3");

  g5g1.add_channel("g5g1");
  g5g2.add_channel("g5g2");
  g5g3.add_channel("g5g3");
  g5gk.add_channel("g5g1");
  g5gk.add_channel("g5g2");
  g5gk.add_channel("g5g3");

  g0g5g1.add_channel("g0g5g1");
  g0g5g2.add_channel("g0g5g2");
  g0g5g3.add_channel("g0g5g3");
  g0g5gk.add_channel("g0g5g1");
  g0g5gk.add_channel("g0g5g2");
  g0g5gk.add_channel("g0g5g3");

  g5_g0g5.add_channel("g5_g0g5_re");


  mpcac.set_ndeps(D_EFF_B1,3);
  mpcac.dep[D_EFF_B1][0]=&g5; mpcac.dep_level[D_EFF_B1][0]=P_COR_B1;
  mpcac.dep[D_EFF_B1][1]=&g5; mpcac.dep_level[D_EFF_B1][1]=P_EFF_B1;
  mpcac.dep[D_EFF_B1][2]=&g5_g0g5; mpcac.dep_level[D_EFF_B1][2]=P_COR_B1;
  mpcac.set_ndeps(D_EFF_B2,3);
  mpcac.dep[D_EFF_B2][0]=&g5; mpcac.dep_level[D_EFF_B2][0]=P_COR_B2;
  mpcac.dep[D_EFF_B2][1]=&g5; mpcac.dep_level[D_EFF_B2][1]=P_EFF_B2;
  mpcac.dep[D_EFF_B2][2]=&g5_g0g5; mpcac.dep_level[D_EFF_B2][2]=P_COR_B2;
  mpcac.set_ndeps(D_FIT,1);
  mpcac.dep[D_FIT][0]=&mpcac; mpcac.dep_level[D_FIT][0]=D_EFF_B2;

  gps.set_ndeps(D_EFF_B1,2);
  gps.dep[D_EFF_B1][0]=&g5; gps.dep_level[D_EFF_B1][0]=P_COR_B1;
  gps.dep[D_EFF_B1][1]=&g5; gps.dep_level[D_EFF_B1][1]=P_EFF_B1;
  gps.set_ndeps(D_EFF_B2,2);
  gps.dep[D_EFF_B2][0]=&g5; gps.dep_level[D_EFF_B2][0]=P_COR_B2;
  gps.dep[D_EFF_B2][1]=&g5; gps.dep_level[D_EFF_B2][1]=P_EFF_B2;
  gps.set_ndeps(D_FIT,1);
  gps.dep[D_FIT][0]=&gps; gps.dep_level[D_FIT][0]=D_EFF_B2;

  fps.set_ndeps(D_EFF_B1,3);
  fps.dep[D_EFF_B1][0]=&g5; fps.dep_level[D_EFF_B1][0]=P_COR_B1;
  fps.dep[D_EFF_B1][1]=&g5; fps.dep_level[D_EFF_B1][1]=P_EFF_B1;
  fps.dep[D_EFF_B1][2]=&mpcac; fps.dep_level[D_EFF_B1][2]=D_EFF_B1;
  fps.set_ndeps(D_EFF_B2,3);
  fps.dep[D_EFF_B2][0]=&g5; fps.dep_level[D_EFF_B2][0]=P_COR_B2;
  fps.dep[D_EFF_B2][1]=&g5; fps.dep_level[D_EFF_B2][1]=P_EFF_B2;
  fps.dep[D_EFF_B2][2]=&mpcac; fps.dep_level[D_EFF_B2][2]=D_EFF_B2;
  fps.set_ndeps(D_FIT,1);
  fps.dep[D_FIT][0]=&fps; fps.dep_level[D_FIT][0]=D_EFF_B2;

  fv.set_ndeps(D_EFF_B1,2);
  fv.dep[D_EFF_B1][0]=&g1; fv.dep_level[D_EFF_B1][0]=P_COR_B1;
  fv.dep[D_EFF_B1][1]=&g1; fv.dep_level[D_EFF_B1][1]=P_EFF_B1;
  fv.set_ndeps(D_EFF_B2,2);
  fv.dep[D_EFF_B2][0]=&g1; fv.dep_level[D_EFF_B2][0]=P_COR_B2;
  fv.dep[D_EFF_B2][1]=&g1; fv.dep_level[D_EFF_B2][1]=P_EFF_B2;
  fv.set_ndeps(D_FIT,1);
  fv.dep[D_FIT][0]=&fv; fv.dep_level[D_FIT][0]=D_EFF_B2;

  fvk.set_ndeps(D_EFF_B1,2);
  fvk.dep[D_EFF_B1][0]=&gk; fvk.dep_level[D_EFF_B1][0]=P_COR_B1;
  fvk.dep[D_EFF_B1][1]=&gk; fvk.dep_level[D_EFF_B1][1]=P_EFF_B1;
  fvk.set_ndeps(D_EFF_B2,2);
  fvk.dep[D_EFF_B2][0]=&gk; fvk.dep_level[D_EFF_B2][0]=P_COR_B2;
  fvk.dep[D_EFF_B2][1]=&gk; fvk.dep_level[D_EFF_B2][1]=P_EFF_B2;
  fvk.set_ndeps(D_FIT,1);
  fvk.dep[D_FIT][0]=&fvk; fvk.dep_level[D_FIT][0]=D_EFF_B2;

	fak.set_ndeps(D_EFF_B1,2);
	fak.dep[D_EFF_B1][0]=&g5gk; fak.dep_level[D_EFF_B1][0]=P_COR_B1;
	fak.dep[D_EFF_B1][1]=&g5gk; fak.dep_level[D_EFF_B1][1]=P_EFF_B1;
	fak.set_ndeps(D_EFF_B2,2);
	fak.dep[D_EFF_B2][0]=&g5gk; fak.dep_level[D_EFF_B2][0]=P_COR_B2;
	fak.dep[D_EFF_B2][1]=&g5gk; fak.dep_level[D_EFF_B2][1]=P_EFF_B2;
	fak.set_ndeps(D_FIT,1);
	fak.dep[D_FIT][0]=&fak; fak.dep_level[D_FIT][0]=D_EFF_B2;


  mvmps.set_ndeps(2);
  mvmps.dep[0]=&g5; mvmps.dep_level[0]=P_FIT;
  mvmps.dep[1]=&g1; mvmps.dep_level[1]=P_FIT;

  mvkmps.set_ndeps(2);
  mvkmps.dep[0]=&g5; mvkmps.dep_level[0]=P_FIT;
  mvkmps.dep[1]=&gk; mvkmps.dep_level[1]=P_FIT;

  mpsfps.set_ndeps(2);
  mpsfps.dep[0]=&fps; mpsfps.dep_level[0]=D_FIT;
  mpsfps.dep[1]=&g5; mpsfps.dep_level[1]=P_FIT;

  mvfps.set_ndeps(2);
  mvfps.dep[0]=&fps; mvfps.dep_level[0]=D_FIT;
  mvfps.dep[1]=&g1; mvfps.dep_level[1]=P_FIT;

  mvkfps.set_ndeps(2);
  mvkfps.dep[0]=&fps; mvkfps.dep_level[0]=D_FIT;
  mvkfps.dep[1]=&gk; mvkfps.dep_level[1]=P_FIT;

  mpsmpcac.set_ndeps(2);
  mpsmpcac.dep[0]=&mpcac; mpsmpcac.dep_level[0]=D_FIT;
  mpsmpcac.dep[1]=&g5; mpsmpcac.dep_level[1]=P_FIT;

  mps2mpcac.set_ndeps(2);
  mps2mpcac.dep[0]=&mpcac; mps2mpcac.dep_level[0]=D_FIT;
  mps2mpcac.dep[1]=&g5; mps2mpcac.dep_level[1]=P_FIT;

  fvfps.set_ndeps(2);
  fvfps.dep[0]=&fps; fvfps.dep_level[0]=D_FIT;
  fvfps.dep[1]=&fv; fvfps.dep_level[1]=D_FIT;

  fvkfps.set_ndeps(2);
  fvkfps.dep[0]=&fps; fvkfps.dep_level[0]=D_FIT;
  fvkfps.dep[1]=&fvk; fvkfps.dep_level[1]=D_FIT;

  mamps.set_ndeps(2);
  mamps.dep[0]=&g5; mamps.dep_level[0]=P_FIT;
  mamps.dep[1]=&g5g1; mamps.dep_level[1]=P_FIT;

  makmps.set_ndeps(2);
  makmps.dep[0]=&g5; makmps.dep_level[0]=P_FIT;
  makmps.dep[1]=&g5gk; makmps.dep_level[1]=P_FIT;

  mamv.set_ndeps(2);
  mamv.dep[0]=&g1; mamv.dep_level[0]=P_FIT;
  mamv.dep[1]=&g5g1; mamv.dep_level[1]=P_FIT;

  makmvk.set_ndeps(2);
  makmvk.dep[0]=&gk; makmvk.dep_level[0]=P_FIT;
  makmvk.dep[1]=&g5gk; makmvk.dep_level[1]=P_FIT;

  gmor.set_ndeps(3);
  gmor.dep[0]=&mpcac; gmor.dep_level[0]=D_FIT;
  gmor.dep[1]=&g5; gmor.dep_level[1]=P_FIT;
  gmor.dep[2]=&fps; gmor.dep_level[2]=D_FIT;


  bool found=false;

  /* CORRELATORS */
#define IF_COR(ch) if(strcmp(arg,#ch"_cor")==0 || strcmp(arg,"all_cor")==0){ch.set_level(P_AVE_COR); found=true;}
  
  IF_COR(id)
  IF_COR(g0)
  IF_COR(g5)
  IF_COR(g0g5)
  IF_COR(g1)
  IF_COR(g2)
  IF_COR(g3)
  IF_COR(gk)
  IF_COR(g0g1)
  IF_COR(g0g2)
  IF_COR(g0g3)
  IF_COR(g0gk)
  IF_COR(g5g1)
  IF_COR(g5g2)
  IF_COR(g5g3)
  IF_COR(g5gk)
  IF_COR(g0g5g1)
  IF_COR(g0g5g2)
  IF_COR(g0g5g3)
  IF_COR(g0g5gk)
  
#undef IF_COR
  
    /* EFFECTIVE MASSES */
#define IF_EFF(ch) if(strcmp(arg,#ch"_eff")==0 || strcmp(arg,"all_eff")==0){ch.set_level(P_AVE_EFF); found=true;}
  
  IF_EFF(id)
  IF_EFF(g0)
  IF_EFF(g5)
  IF_EFF(g0g5)
  IF_EFF(g1)
  IF_EFF(g2)
  IF_EFF(g3)
  IF_EFF(gk)
  IF_EFF(g0g1)
  IF_EFF(g0g2)
  IF_EFF(g0g3)
  IF_EFF(g0gk)
  IF_EFF(g5g1)
  IF_EFF(g5g2)
  IF_EFF(g5g3)
  IF_EFF(g5gk)
  IF_EFF(g0g5g1)
  IF_EFF(g0g5g2)
  IF_EFF(g0g5g3)
  IF_EFF(g0g5gk)
  
#undef IF_EFF

      /* FITS */

#define IF_FIT(ch) if(strcmp(arg,#ch)==0 || strcmp(arg,"all_fit")==0){ch.set_level(P_AVE_FIT); found=true;}
  
  IF_FIT(id)
  IF_FIT(g0)
  IF_FIT(g5)
  IF_FIT(g0g5)
  IF_FIT(g1)
  IF_FIT(g2)
  IF_FIT(g3)
  IF_FIT(gk)
  IF_FIT(g0g1)
  IF_FIT(g0g2)
  IF_FIT(g0g3)
  IF_FIT(g0gk)
  IF_FIT(g5g1)
  IF_FIT(g5g2)
  IF_FIT(g5g3)
  IF_FIT(g5gk)
  IF_FIT(g0g5g1)
  IF_FIT(g0g5g2)
  IF_FIT(g0g5g3)
  IF_FIT(g0g5gk)
  
#undef IF_FIT

	/* DERIVED */
#define IF_DER_EFF(ch) if(strcmp(arg,#ch"_eff")==0 || strcmp(arg,"all_der_eff")==0){ch.set_level(D_AVE_EFF); found=true;}

  IF_DER_EFF(mpcac)
  IF_DER_EFF(fps)
  IF_DER_EFF(gps)
  IF_DER_EFF(fv)
  IF_DER_EFF(fvk)
  IF_DER_EFF(fak)
  
#undef IF_DER_EFF
  
#define IF_DER_FIT(ch) if(strcmp(arg,#ch)==0 || strcmp(arg,"all_der_fit")==0){ch.set_level(D_AVE_FIT); found=true;}
    
  IF_DER_FIT(mpcac)
  IF_DER_FIT(fps)
  IF_DER_FIT(gps)
  IF_DER_FIT(fv)
  IF_DER_FIT(fvk)
  IF_DER_FIT(fak)
  
#undef IF_DER_FIT
  
	    /* RATIOS */
#define IF_RATIO(ch) if(strcmp(arg,#ch)==0 || strcmp(arg,"all_ratio")==0){ch.set_level(R_ACTIVE); found=true;}

    IF_RATIO(mvmps);
    IF_RATIO(mvkmps);
    IF_RATIO(mpsfps);
    IF_RATIO(mvfps);
    IF_RATIO(mvkfps);
    IF_RATIO(gmor);
    IF_RATIO(mpsmpcac);
    IF_RATIO(mps2mpcac);
    IF_RATIO(fvfps);
    IF_RATIO(fvkfps);
    IF_RATIO(mamps);
    IF_RATIO(makmps);
    IF_RATIO(mamv);
    IF_RATIO(makmvk);
    
#undef IF_RATIO

	    return found;
}

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

void eval_ctrl::normalize_cuts(int lt) {
  map<string,int>::iterator end=left_cut.end();
  map<string,int>::iterator g1=left_cut.find("g1");
  map<string,int>::iterator g2=left_cut.find("g2");
  map<string,int>::iterator g3=left_cut.find("g3");
  map<string,int>::iterator gk=left_cut.find("gk");
  if(gk==end && g1!=end && g2!=end && g3!=end){
    left_cut["gk"]=max((*g1).second,max((*g2).second,(*g3).second));
    right_cut["gk"]=min(right_cut["g1"],min(right_cut["g2"],right_cut["g3"]));
  }
  else if(gk!=end && g1==end && g2==end && g3==end){
    left_cut["g1"]=left_cut["g2"]=left_cut["g3"]=left_cut["gk"];
    right_cut["g1"]=right_cut["g2"]=right_cut["g3"]=right_cut["gk"];
  }
  
  g1=left_cut.find("g5g1");
  g2=left_cut.find("g5g2");
  g3=left_cut.find("g5g3");
  gk=left_cut.find("g5gk");
  if(gk==end && g1!=end && g2!=end && g3!=end){
    left_cut["g5gk"]=max((*g1).second,max((*g2).second,(*g3).second));
    right_cut["g5gk"]=min(right_cut["g5g1"],min(right_cut["g5g2"],right_cut["g5g3"]));
  }
  else if(gk!=end && g1==end && g2==end && g3==end){
    left_cut["g5g1"]=left_cut["g5g2"]=left_cut["g5g3"]=left_cut["g5gk"];
    right_cut["g5g1"]=right_cut["g5g2"]=right_cut["g5g3"]=right_cut["g5gk"];
  }
  
  g1=left_cut.find("g0g1");
  g2=left_cut.find("g0g2");
  g3=left_cut.find("g0g3");
  gk=left_cut.find("g0gk");
  if(gk==end && g1!=end && g2!=end && g3!=end){
    left_cut["g0gk"]=max((*g1).second,max((*g2).second,(*g3).second));
    right_cut["g0gk"]=min(right_cut["g0g1"],min(right_cut["g0g2"],right_cut["g0g3"]));
  }
  else if(gk!=end && g1==end && g2==end && g3==end){
    left_cut["g0g1"]=left_cut["g0g2"]=left_cut["g0g3"]=left_cut["g0gk"];
    right_cut["g0g1"]=right_cut["g0g2"]=right_cut["g0g3"]=right_cut["g0gk"];
  }
  
  g1=left_cut.find("g0g5g1");
  g2=left_cut.find("g0g5g2");
  g3=left_cut.find("g0g5g3");
  gk=left_cut.find("g0g5gk");
  if(gk==end && g1!=end && g2!=end && g3!=end){
    left_cut["g0g5gk"]=max((*g1).second,max((*g2).second,(*g3).second));
    right_cut["g0g5gk"]=min(right_cut["g0g5g1"],min(right_cut["g0g5g2"],right_cut["g0g5g3"]));
  }
  else if(gk!=end && g1==end && g2==end && g3==end){
    left_cut["g0g5g1"]=left_cut["g0g5g2"]=left_cut["g0g5g3"]=left_cut["g0g5gk"];
    right_cut["g0g5g1"]=right_cut["g0g5g2"]=right_cut["g0g5g3"]=right_cut["g0g5gk"];
  }
 


  for(set<primary_ctrl*>::iterator pctrl_it=primary_ctrl::all_pr.begin();pctrl_it!=primary_ctrl::all_pr.end();pctrl_it++) {
    if((*pctrl_it)->get_level(P_EFF_B1) || (*pctrl_it)->get_level(P_EFF_B2)) {
      if(left_cut.find((*pctrl_it)->name) == left_cut.end())  {
        cerr << "[normalize_cuts] Cuts not defined for " << (*pctrl_it)->name << endl;
        exit(1);
      }
    }
  }

  for(set<derived_ctrl*>::iterator dctrl_it=derived_ctrl::all_der.begin();dctrl_it!=derived_ctrl::all_der.end();dctrl_it++) {
    if((*dctrl_it)->get_level(D_EFF_B1) || (*dctrl_it)->get_level(D_EFF_B2)) {
      if(left_cut.find((*dctrl_it)->name) == left_cut.end())  {
        cerr << "[normalize_cuts Cuts not defined for " << (*dctrl_it)->name << endl;
        exit(1);
      }
    }
  }
  


  for(set<eval_ctrl*>::iterator ctrl_it=eval_ctrl::all.begin();ctrl_it!=eval_ctrl::all.end();ctrl_it++)
    (*ctrl_it)->prop_cuts();

  for(set<eval_ctrl*>::iterator ctrl_it=eval_ctrl::all.begin();ctrl_it!=eval_ctrl::all.end();ctrl_it++) {
    if( ((*ctrl_it)->isprimary() && ((*ctrl_it)->get_level(P_EFF_B1) || (*ctrl_it)->get_level(P_EFF_B2)) )
	|| ((*ctrl_it)->isderived() && ((*ctrl_it)->get_level(D_EFF_B1) || (*ctrl_it)->get_level(D_EFF_B2)) ) )
      {
	if((*ctrl_it)->get_left_cut() < 0) {
	  cerr << "[normalize_cuts] Left cut for " << (*ctrl_it)->name
	       << " less than 0" << endl;
	  exit(1);
	}
	if((*ctrl_it)->get_right_cut() > lt/2) {
	  cerr << "[normalize_cuts] Right cut for " << (*ctrl_it)->name
	       << " greater than " << lt/2 << endl;
	  exit(1);
	}
      }
  }

}

#undef min
#undef max
primary_ctrl::primary_ctrl(const char* n, int nc){
  for(int i=0;i<MAX_NLEVELS;i++)
    level[i]=false;

  if(nc<=0){
    cerr<<"[primary_ctrl]: Wrong initilization of the operator: nc <= 0"<< endl;
    exit(1);
  }
  name=n;
  nchan=nc;

  channel=NULL;
  chan_counter=0;
  
  channel=new char*[nchan];
  for(int i=0;i<nchan;i++)
    channel[i]=new char[32];

  all.insert(this);
  all_pr.insert(this);
  
  cor_b0=NULL;
  cor_b1=NULL;
  
  ws = NULL;

/*  cout << "3 PRIMARY_CTRL constructor " << name << " nchan=" << nchan << endl;*/
}


derived_ctrl::derived_ctrl(const char* n){
  for(int i=0;i<MAX_NLEVELS;i++)
    level[i]=false;

  name=n;
  
  for(int i=D_EFF_B1; i<=D_FIT; i++) {
    ndep[i]=0;
    dep[i]=NULL;
    dep_level[i]=NULL;
  }

  all.insert(this);
  all_der.insert(this);

/*  cout << "3 DERIVED_CTRL constructor " << name << endl;*/
}


ratio_ctrl::ratio_ctrl(const char* n){
  for(int i=0;i<MAX_NLEVELS;i++)
    level[i]=false;

  name=n;

  ndep=0;
  dep=NULL;
  dep_level=NULL;

  all.insert(this);

/*  cout << "3 RATIO_CTRL constructor " << name << endl;*/
}



primary_ctrl::~primary_ctrl(){
  if(nchan>0){
    for(int i=0;i<nchan;i++)
      delete[] channel[i];

    delete[] channel;
  }
  if(cor_b0!=NULL) delete[] cor_b0;
  if(cor_b1!=NULL) delete cor_b1;

  if(ws != NULL) {
    delete[] ws[0];
    delete[] ws[1];
    delete[] ws[2];
    delete[] ws[3];
    delete[] ws;
  }

  all.erase(this);
  all_pr.erase(this);
}


derived_ctrl::~derived_ctrl(){
  for(int i=D_EFF_B1; i<=D_FIT; i++) {
    if(ndep[i]>0){
      delete[] dep_level[i];
      delete[] dep[i];
    }
  }

  all.erase(this);
  all_der.erase(this);
}


ratio_ctrl::~ratio_ctrl(){
  if(ndep>0){
    delete[] dep_level;
    delete[] dep;
  }

  all.erase(this);
}


void primary_ctrl::set_level(int lv){
  if(lv<P_READ || lv>P_AVE_FIT){
    cerr<<"[primary_ctrl:set_level]:Level must be in the range lv>=P_READ && lv<=P_AVE_FIT for primary operators"<< endl;
    exit(1);
  }

  if(!level[lv]) {
    level[P_READ]=true;
    level[lv]=true;

    if(lv==P_AVE_FIT) set_level(P_FIT);
    if(lv==P_AVE_EFF) set_level(P_EFF_B1);
    if(lv==P_AVE_COR) set_level(P_COR_B1);

    if(lv==P_FIT) set_level(P_EFF_B2);
    if(lv==P_EFF_B2) set_level(P_COR_B2);
    if(lv==P_COR_B2) bs_2nd=true;

    if(lv==P_EFF_B1) set_level(P_COR_B1);

    prop_dep();
  }
}


void derived_ctrl::set_level(int lv){
  if(lv<D_EFF_B1 || lv>D_AVE_FIT){
    cerr<<"[derived_ctrl:set_level]:Level must be lv>=D_EFF_B1 && lv<=D_AVE_FIT for derived operators"<< endl;
    exit(1);
  }
  
  if(!level[lv]) {
    level[lv]=true;

    if(lv==D_AVE_FIT) set_level(D_FIT);
    if(lv==D_AVE_EFF) set_level(D_EFF_B1);

    if(lv==D_FIT) set_level(D_EFF_B2);
    if(lv==D_EFF_B2) bs_2nd=true;

    prop_dep();
  }
}

void ratio_ctrl::set_level(int lv=0){
  if(lv!=R_ACTIVE){
    cerr<<"[ratio_ctrl:set_level]:Level must be lv==R_ACTIVE for ratio operators"<< endl;
    exit(1);
  }

  if(!level[lv]) {
    level[lv]=true;
    prop_dep();
  }
}


void primary_ctrl::add_channel(const char * c1){
  if(chan_counter==nchan){
    cerr<<"[primary_ctrl:add_channel]:Number of channels (nchan) too small"<<endl;
    exit(1);
  }
  strcpy(channel[chan_counter],c1);
  chan_counter++;

  cout << "3 PRIMARY_CTRL add_channel"
       << " name=" << name
       << " channel[" << chan_counter << "]=" << c1
       << endl;
}


void primary_ctrl::prop_dep() {
}


void derived_ctrl::prop_dep() {
  for(int k=D_EFF_B1;k<=D_FIT;k++) {
    if(level[k]) {
      for(int i=0; i< ndep[k];i++)
        dep[k][i]->set_level(dep_level[k][i]);
    }
  }
}


void ratio_ctrl::prop_dep() {
  for(int i=0; i< ndep;i++)
    dep[i]->set_level(dep_level[i]);
}


void derived_ctrl::prop_cuts() {
  for(int k=D_EFF_B1;k<=D_FIT;k++) {
    if(level[k]) {
      for(int i=0; i< ndep[k];i++) {
        if( (dep[k][i]->isprimary() && (dep_level[k][i]==P_EFF_B1 || dep_level[k][i]==P_EFF_B2))
	    || (dep[k][i]->isderived() && (dep_level[k][i]==D_EFF_B1 || dep_level[k][i]==D_EFF_B2)) )
	  {
	    if(left_cut[name] < dep[k][i]->get_left_cut())
	      left_cut[name] = dep[k][i]->get_left_cut();
	    if(right_cut[name] > dep[k][i]->get_right_cut())
	      right_cut[name] = dep[k][i]->get_right_cut();
	  }
      }
    }
  }
}



ostream& operator<<(ostream& os, eval_ctrl& ev) {
  os << ev.name;
  ev.print_utility(os);
  return os;
}


ostream& primary_ctrl::print_utility(ostream& os) {
  return os
    << " P_READ " << level[P_READ]
    << " P_COR_B1 " << level[P_COR_B1]
    << " P_EFF_B1 " << level[P_EFF_B1]
    << " P_COR_B2 " << level[P_COR_B2]
    << " P_EFF_B2 " << level[P_EFF_B2]
    << " P_FIT " << level[P_FIT]
    << " P_AVE_COR " << level[P_AVE_COR]
    << " P_AVE_EFF " << level[P_AVE_EFF]
    << " P_AVE_FIT " << level[P_AVE_FIT];
}


ostream& derived_ctrl::print_utility(ostream& os) {
  return os
    << " D_EFF_B1 " << level[D_EFF_B1]
    << " D_EFF_B2 " << level[D_EFF_B2]
    << " D_FIT " << level[D_FIT]
    << " D_AVE_EFF " << level[D_AVE_EFF]
    << " D_AVE_FIT " << level[D_AVE_FIT];
}


ostream& ratio_ctrl::print_utility(ostream& os) {
  return os
    << " R_ACTIVE " << level[R_ACTIVE];
}



void primary_ctrl::allocate_datamemory(int lt) {
  ws = new double*[4];
  ws[0] = new double[lt/2+1];
  ws[1] = new double[lt/2+1];
  ws[2] = new double[lt/2+1];
  ws[3] = new double[lt/2+1];

  if(level[P_READ]) {
    if(cor_b0 == NULL) {
      cor_b0=new Corr_t*[nchan];
      for(int nc=0;nc<nchan;nc++){
        if(channel_map.find(channel[nc])==channel_map.end())
          channel_map[channel[nc]]=new Corr_t(lt);
        cor_b0[nc]=channel_map[channel[nc]];
      }
    }
  }
  
  if(level[P_COR_B1]) {
    if(cor_b1 == NULL) cor_b1=new Corr_t(lt);
  }
  
  if(level[P_EFF_B1]) {
    if(eff_b1 == NULL) eff_b1=new Corr_t(lt);
  }
  
  if(level[P_COR_B2]) {
    if(cor_b2 == NULL) cor_b2=new Corr_t(lt);
  }
  
  if(level[P_EFF_B2]) {
    if(eff_b2 == NULL) eff_b2=new Corr_t(lt);
  }
}


void derived_ctrl::allocate_datamemory(int lt) {
  if(level[D_EFF_B1]) {
    if(eff_b1 == NULL) eff_b1=new Corr_t(lt);
  }
  
  if(level[D_EFF_B2]) {
    if(eff_b2 == NULL) eff_b2=new Corr_t(lt);
  }
}


void ratio_ctrl::allocate_datamemory(int lt) {
}


void derived_ctrl::set_ndeps(int lv, int nd){
  if(lv<D_EFF_B1 || lv>D_FIT){
    cerr<<"[derived_ctrl::set_ndep]: Wrong level: lv<D_EFF_B1 || lv>D_FIT"<< endl;
    exit(1);
  }
  if(nd<=0){
    cerr<<"[derived_ctrl::set_ndep]: Wrong initilization of the operator: nd <= 0"<< endl;
    exit(1);
  }
  if(ndep[lv]!=0){
    cerr<<"[derived_ctrl::set_ndep]: Level " << lv << " already initialized" << endl;
    exit(1);
  }

  ndep[lv]=nd;
  dep[lv]=new eval_ctrl*[nd];
  dep_level[lv]=new int[nd];
  for(int i=0; i<nd; i++)
    dep_level[lv][i]=-1;
}


void ratio_ctrl::set_ndeps(int nd){
  if(nd<=0){
    cerr<<"[ratio_ctrl::set_ndep]: Wrong initilization of the operator: nd <= 0"<< endl;
    exit(1);
  }
  if(ndep!=0){
    cerr<<"[ratio_ctrl::set_ndep]: Already initialized" << endl;
    exit(1);
  }

  ndep=nd;
  dep=new eval_ctrl*[nd];
  dep_level=new int[nd];
  for(int i=0; i<nd; i++)
    dep_level[i]=-1;
}


void primary_ctrl::purge_b2() {
  if(cor_b2 != NULL) cor_b2->purge();
  if(eff_b2 != NULL) eff_b2->purge();
}


void derived_ctrl::purge_b2() {
  if(eff_b2 != NULL) eff_b2->purge();
}



void primary_ctrl::print() {
  if(level[P_AVE_COR]) {
    for(int t=0; t<cor_b1->length; t++)
      cout << "10 COR " << name << " " << t << " " << cor_b1->d[t].avr().val << " " << cor_b1->d[t].stderr().val << " blksize " << blsize << endl;
  }
  
  if(level[P_AVE_EFF]) {
    for(int t=left_cut[name]; t<=right_cut[name]; t++)
      cout << "10 EFF " << name << " " << t << " " << eff_b1->d[t].avr().val << " " << eff_b1->d[t].stderr().val << " blksize " << blsize << endl;
  }
  
  if(level[P_AVE_FIT])
    cout << "10 FIT " << name << " " << fit.avr().val << " " << fit.stderr().val << " blksize " << blsize << " lcut " << left_cut[name] << " rcut " << right_cut[name] << endl;
}


void derived_ctrl::print() {
  if(level[D_AVE_EFF]) {
    for(int t=left_cut[name]; t<=right_cut[name]; t++)
      cout << "10 EFF " << name << " " << t << " " << eff_b1->d[t].avr().val << " " << eff_b1->d[t].stderr().val << " blksize " << blsize << endl;
  }
  
  if(level[D_AVE_FIT])
    cout << "10 FIT " << name << " " << fit.avr().val << " " << fit.stderr().val << " blksize " << blsize << " lcut " << left_cut[name] << " rcut " << right_cut[name] << endl;
}

void ratio_ctrl::print() {
  if(level[R_ACTIVE]) {
    cout << "10 FIT " << name << " " << fit.avr().val << " " << fit.stderr().val << " blksize " << blsize << " lcut " << left_cut[name] << " rcut " << right_cut[name] << endl;
  }
}


int eval_ctrl::get_left_cut() {
  if(left_cut.find(name) == left_cut.end()) return 0;
  return left_cut[name];
}

int eval_ctrl::get_right_cut() {
  if(right_cut.find(name) == right_cut.end()) return 0;
  return right_cut[name];
}



bool primary_ctrl::eval_cor_eff(const int bsl, const int* bs, const int len, int Lt, int effmass_method) {
  double *ws_cor;
  double *ws_eff;
  int P_COR;
  int P_EFF;
  
  if(bsl==1) {
    P_COR=P_COR_B1;
    P_EFF=P_EFF_B1;
    ws_cor=ws[0];
    ws_eff=ws[1];
  } else if(bsl==2) {
    P_COR=P_COR_B2;
    P_EFF=P_EFF_B2;
    ws_cor=ws[2];
    ws_eff=ws[3];
  } else {
    cerr<<"[primary_ctrl::eval_cor_eff]: bsl = " << bsl << " can be only 1 or 2" << endl;
    exit(1);
  }

  if(get_level(P_COR)) {
    for (int i=0;i<cor_b0[0]->length;++i) {
      double ave=0.;
      for (int k=0;k<len;++k)
        for(int nc=0; nc<nchan; nc++)
          ave+=cor_b0[nc]->d[i][bs[k]];
      ave /= static_cast<double>(len*(nchan));
      ws_cor[i]=ave;
    }
  }

  if(get_level(P_EFF)) {
    if(effmass_method==0) {
      for (int i=get_left_cut(); i<=get_right_cut(); i++){
        ws_eff[i]=plain_eff_mass(
				 Lt/2-i,ws_cor[(i-1+Lt)%Lt]/ws_cor[i]);
        if (ws_eff[i]<0.) return false;
      }
    } else if(effmass_method==1) {
      double m1;
      int ret;
      for (int i=get_left_cut(); i<=get_right_cut(); i++){
        ret=shifted_prony_eff_mass_1(
				     ws_cor,i,get_right_cut(),&m1,Lt);
        if (ret==0) return false;
        ws_eff[i]=m1;
      }
    } else if(effmass_method==2) {
      double m1,m2;
      int ret;
      for (int i=get_left_cut(); i<=get_right_cut(); i++){
        ret=shifted_prony_eff_mass_2(
				     ws_cor,i,get_right_cut(),&m1,&m2,Lt);
        if (ret==0) return false;
        ws_eff[i]=m1;
      }
    }
  }

  return true;
}




void primary_ctrl::store_cor_eff(const int bsl) {
  Corr_t *cor;
  Corr_t *eff;
  double *ws_cor;
  double *ws_eff;
  int P_COR;
  int P_EFF;
  
  if(bsl==1) {
    P_COR=P_COR_B1;
    P_EFF=P_EFF_B1;
    cor=cor_b1;
    eff=eff_b1;
    ws_cor=ws[0];
    ws_eff=ws[1];
  } else if(bsl==2) {
    P_COR=P_COR_B2;
    P_EFF=P_EFF_B2;
    cor=cor_b2;
    eff=eff_b2;
    ws_cor=ws[2];
    ws_eff=ws[3];
  } else {
    cerr<<"[primary_ctrl::store_cor_eff]: bsl = " << bsl << " can be only 1 or 2" << endl;
    exit(1);
  }

  if(get_level(P_COR)) {
    for (int i=0; i<cor_b0[0]->length; i++)
      cor->d[i].push_back(ws_cor[i]);
  }

  if(get_level(P_EFF)) {
    for (int i=get_left_cut(); i<=get_right_cut(); i++)
      eff->d[i].push_back(ws_eff[i]);
  }

}
