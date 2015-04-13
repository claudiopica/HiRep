#include "bs_type.h"
#include "bs_reader.h"
#include "bs_ctrl.h"
#include "datasample.h"
#include "effective_mass.h"
#include "bootstrap.h"
#include "progbar.h"
#include "bs_observables.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <stdlib.h>

using namespace std;

int Lt=0;
int Ls=0;
double svol=.0;
int blsize=0;
int nsamples1=1000;
int nsamples2=100;
string channel;
string inputfilename;
string cutfilename;
int effmass_method;


extern primary_ctrl id;
extern primary_ctrl g0;
extern primary_ctrl g5;
extern primary_ctrl g0g5;
extern primary_ctrl g1;
extern primary_ctrl g2;
extern primary_ctrl g3;
extern primary_ctrl gk;
extern primary_ctrl g0g1;
extern primary_ctrl g0g2;
extern primary_ctrl g0g3;
extern primary_ctrl g0gk;
extern primary_ctrl g5g1;
extern primary_ctrl g5g2;
extern primary_ctrl g5g3;
extern primary_ctrl g5gk;
extern primary_ctrl g0g5g1;
extern primary_ctrl g0g5g2;
extern primary_ctrl g0g5g3;
extern primary_ctrl g0g5gk;
extern primary_ctrl g5_g0g5;
extern derived_ctrl mpcac;
extern derived_ctrl gps;
extern derived_ctrl fps;
extern derived_ctrl fv;
extern derived_ctrl fvk;
extern derived_ctrl fak;
extern ratio_ctrl mvmps;
extern ratio_ctrl mvkmps;
extern ratio_ctrl mpsfps;
extern ratio_ctrl mvfps;
extern ratio_ctrl mvkfps;
extern ratio_ctrl gmor;
extern ratio_ctrl mpsmpcac;
extern ratio_ctrl mps2mpcac;
extern ratio_ctrl fvfps;
extern ratio_ctrl fvkfps;
extern ratio_ctrl mamps;
extern ratio_ctrl makmps;
extern ratio_ctrl mamv;
extern ratio_ctrl makmvk;


int main(int argc, char* argv[]) {
  int cmdline=read_cmdline(argc, argv);
  
  if(cmdline == 0 ){
    cout << "1 CMDLINE channel " << channel << "\n";
    cout << "1 CMDLINE inputfile " << inputfilename << "\n";
    cout << "1 CMDLINE cutfile " << cutfilename << "\n";
    cout << "1 CMDLINE Lt " << Lt << "\n";
    cout << "1 CMDLINE Ls " << Ls << "\n";
    cout << "1 CMDLINE blocksize " << blsize << "\n";
    cout << "1 CMDLINE nsamples1 " << nsamples1 << "\n";
    cout << "1 CMDLINE nsamples2 " << nsamples2 << "\n";
  }
  
  if(!eval_ctrl::fill_dep(channel.c_str())){
    cerr << "[XXX]: No recognized channel ("<< channel <<")"<<endl;
    exit(1);
  }
  
  for(set<eval_ctrl*>::iterator ctrl_it=eval_ctrl::all.begin();ctrl_it!=eval_ctrl::all.end();ctrl_it++)
    cout << "2 EVAL_CTRL " << **ctrl_it << "\n";

  for(set<primary_ctrl*>::iterator pctrl_it=primary_ctrl::all_pr.begin();pctrl_it!=primary_ctrl::all_pr.end();pctrl_it++)
    cout << "2 PRIMARY_CTRL " << **pctrl_it <<endl;


  for(set<eval_ctrl*>::iterator ctrl_it=eval_ctrl::all.begin();ctrl_it!=eval_ctrl::all.end();ctrl_it++)
    (*ctrl_it)->allocate_datamemory(Lt);

  if(eval_ctrl::bs_2nd)
    cout << "2 EVAL_CTRL bs_2nd = true\n";
  else
    cout << "2 EVAL_CTRL bs_2nd = false\n";
  cout << "2 EVAL_CTRL nwanted = " << eval_ctrl::channel_map.size() << "\n";
  for(map<string,Corr_t*>::iterator chan_it=eval_ctrl::channel_map.begin(); chan_it!=eval_ctrl::channel_map.end(); chan_it++)
    cout << "2 EVAL_CTRL wanted " << (*chan_it).first<< endl;

  if(cmdline != 0) exit(0);

  map<string,int>::iterator cut_it;
  int csize=channel.size(),comparesize=4;
  if(csize-comparesize<0) comparesize=csize;

   if( channel.compare(csize-comparesize,comparesize,"_cor")!=0) { 
    read_cut(cutfilename.c_str());
    eval_ctrl::normalize_cuts(Lt);
    for(cut_it=eval_ctrl::left_cut.begin();cut_it!=eval_ctrl::left_cut.end();cut_it++)
      cout << "2 READCUT " << (*cut_it).first << " " << (*cut_it).second << " " << eval_ctrl::right_cut[(*cut_it).first] << "\n";
   }  

   cout << "3 READING " << inputfilename;
   read_input(inputfilename.c_str());
   cout << " DONE"<<endl;


 /* for(map<string,Corr_t*>::iterator chan_it=eval_ctrl::channel_map.begin(); chan_it!=eval_ctrl::channel_map.end(); chan_it++) {
    cout << "4 CORR " << (*chan_it).first;
    for(int t=0; t<(*chan_it).second->length; t++)
      cout << " " << (*chan_it).second->d[t][0];
    cout << endl;
  }*/



  ///////////////////////////////////////////////////////////////////////////
  //set effective lenght for bootstrap samples
  ///////////////////////////////////////////////////////////////////////////
  unsigned int len=eval_ctrl::channel_map.begin()->second->d[0].size();

  int elen=(int)(double(len)/double(blsize));
  //cerr<<"elen="<<elen<<"len="<<len<<"\n";

  ///////////////////////////////////////////////////////////////////////////
  //boostrap
  ///////////////////////////////////////////////////////////////////////////
  int count1=0;
  for (int sample1=0;sample1<nsamples1;++sample1) {
    ++count1;
    DrawBar(double(sample1)/double(nsamples1-1)*100.);

    //check the efficiency to avoide infinite loops
    if(count1>(20*nsamples1)/100 && double(sample1)/count1 < 0.6){
      cerr << endl;
      cout << "20 EFFICIENCY_BS1 0 %" << endl;
      exit(1);
    }



    //create bootstrap sample
    int bs[elen];
    b_sample(len, elen, bs);
      
    bool valid_sample1=true;

    for(set<primary_ctrl*>::iterator pctrl_it=primary_ctrl::all_pr.begin();pctrl_it!=primary_ctrl::all_pr.end();pctrl_it++) {
      valid_sample1=(*pctrl_it)->eval_cor_eff(1,bs,elen,Lt,effmass_method);
      if (valid_sample1==false) break;
    }
  
    if(valid_sample1 && eval_ctrl::bs_2nd){
      int trys=0;
      
      for(set<eval_ctrl*>::iterator ctrl_it=eval_ctrl::all.begin();ctrl_it!=eval_ctrl::all.end();ctrl_it++)
        (*ctrl_it)->purge_b2();
      
      for (int sample2=0;sample2<nsamples2 && trys<100;++sample2) {
        int bs_tmp[elen], bs2[elen];
        b_sample(elen, elen, bs_tmp);
        for(int i=0;i<elen;i++) bs2[i]=bs[bs_tmp[i]];

        bool valid_sample2=true;

        for(set<primary_ctrl*>::iterator pctrl_it=primary_ctrl::all_pr.begin();pctrl_it!=primary_ctrl::all_pr.end();pctrl_it++) {
          valid_sample2=(*pctrl_it)->eval_cor_eff(2,bs2,elen,Lt,effmass_method);
          if (valid_sample2==false) break;
        }

        if (valid_sample2) {
          trys=0;

          for(set<primary_ctrl*>::iterator pctrl_it=primary_ctrl::all_pr.begin();pctrl_it!=primary_ctrl::all_pr.end();pctrl_it++) {
            (*pctrl_it)->store_cor_eff(2);
          }
            
          if(mpcac.get_level(D_EFF_B2))
            mpcac_eval(mpcac.eff_b2,g5.cor_b2,g5.eff_b2,g5_g0g5.cor_b2);

          if(gps.get_level(D_EFF_B2))
            gps_eval(gps.eff_b2,g5.cor_b2,g5.eff_b2);

          if(fps.get_level(D_EFF_B2))
            fps_eval(fps.eff_b2,g5.cor_b2,g5.eff_b2,mpcac.eff_b2);

          if(fv.get_level(D_EFF_B2))
            fv_eval(fv.eff_b2,g1.cor_b2,g1.eff_b2);

          if(fvk.get_level(D_EFF_B2))
            fvk_eval(fvk.eff_b2,gk.cor_b2,gk.eff_b2);

			  if(fak.get_level(D_EFF_B2))
				  fak_eval(fak.eff_b2,g5gk.cor_b2,g5gk.eff_b2);

        } else { --sample2; ++trys; }
      }

      if(trys<100) {
        for(set<primary_ctrl*>::iterator pctrl_it=primary_ctrl::all_pr.begin();pctrl_it!=primary_ctrl::all_pr.end();pctrl_it++) {
          if((*pctrl_it)->get_level(P_FIT))
            (*pctrl_it)->fit.push_back(fit((*pctrl_it)->get_left_cut(),(*pctrl_it)->get_right_cut(),(*pctrl_it)->eff_b2));
        }

        for(set<derived_ctrl*>::iterator dctrl_it=derived_ctrl::all_der.begin();dctrl_it!=derived_ctrl::all_der.end();dctrl_it++) {
          if((*dctrl_it)->get_level(D_FIT))
            (*dctrl_it)->fit.push_back(fit((*dctrl_it)->get_left_cut(),(*dctrl_it)->get_right_cut(),(*dctrl_it)->eff_b2));
        }

        if(mvmps.get_level(R_ACTIVE)) {
          mvmps.fit.push_back(g1.fit.back()/g5.fit.back());
        }

        if(mvkmps.get_level(R_ACTIVE)) {
          mvkmps.fit.push_back(gk.fit.back()/g5.fit.back());
        }

        if(mpsfps.get_level(R_ACTIVE)) {
          mpsfps.fit.push_back(g5.fit.back()/fps.fit.back());
        }

        if(mvfps.get_level(R_ACTIVE)) {
          mvfps.fit.push_back(g1.fit.back()/fps.fit.back());
        }

        if(mvkfps.get_level(R_ACTIVE)) {
          mvkfps.fit.push_back(gk.fit.back()/fps.fit.back());
        }

        if(mpsmpcac.get_level(R_ACTIVE)) {
          mpsmpcac.fit.push_back(g5.fit.back()/mpcac.fit.back());
        }

        if(mps2mpcac.get_level(R_ACTIVE)) {
          mps2mpcac.fit.push_back(g5.fit.back()*g5.fit.back()/mpcac.fit.back());
        }

        if(fvfps.get_level(R_ACTIVE)) {
          fvfps.fit.push_back(fv.fit.back()/fps.fit.back());
        }

        if(fvkfps.get_level(R_ACTIVE)) {
          fvkfps.fit.push_back(fvk.fit.back()/fps.fit.back());
        }

        if(mamps.get_level(R_ACTIVE)) {
          mamps.fit.push_back(g5g1.fit.back()/g5.fit.back());
        }

        if(makmps.get_level(R_ACTIVE)) {
          makmps.fit.push_back(g5gk.fit.back()/g5.fit.back());
        }

        if(mamv.get_level(R_ACTIVE)) {
           mamv.fit.push_back(g5g1.fit.back()/g1.fit.back());
        }

        if(makmvk.get_level(R_ACTIVE)) {
           makmvk.fit.push_back(g5gk.fit.back()/gk.fit.back());
        }
 
        if(gmor.get_level(R_ACTIVE)) {
          gmor.fit.push_back(g5.fit.back()*g5.fit.back()*fps.fit.back()*fps.fit.back()/mpcac.fit.back());
        }
        
      } else
        valid_sample1 = false;
    }
    
    if (valid_sample1) {
      for(set<primary_ctrl*>::iterator pctrl_it=primary_ctrl::all_pr.begin();pctrl_it!=primary_ctrl::all_pr.end();pctrl_it++) {
        (*pctrl_it)->store_cor_eff(1);
      }
           
      if(mpcac.get_level(D_EFF_B1))
        mpcac_eval(mpcac.eff_b1,g5.cor_b1,g5.eff_b1,g5_g0g5.cor_b1);

      if(gps.get_level(D_EFF_B1))
        gps_eval(gps.eff_b1,g5.cor_b1,g5.eff_b1);

      if(fps.get_level(D_EFF_B1))
        fps_eval(fps.eff_b1,g5.cor_b1,g5.eff_b1,mpcac.eff_b1);

      if(fv.get_level(D_EFF_B1))
        fv_eval(fv.eff_b1,g1.cor_b1,g1.eff_b1);

      if(fvk.get_level(D_EFF_B1))
        fvk_eval(fvk.eff_b1,gk.cor_b1,gk.eff_b1);

		 if(fak.get_level(D_EFF_B1))
			 fak_eval(fak.eff_b1,g5gk.cor_b1,g5gk.eff_b1);

    } else --sample1;
        
  }
  
  cout << "20 EFFICIENCY_BS1 " << (nsamples1*100)/count1 << " %" << endl;

  ///////////////////////////////////////////////////////////////////////////
  //write output
  ///////////////////////////////////////////////////////////////////////////

  for(set<eval_ctrl*>::iterator ctrl_it=eval_ctrl::all.begin();ctrl_it!=eval_ctrl::all.end();ctrl_it++)
    (*ctrl_it)->print();


  for(map<string,Corr_t*>::iterator chan_it=eval_ctrl::channel_map.begin(); chan_it!=eval_ctrl::channel_map.end(); chan_it++)
    delete (*chan_it).second;

  return 0;
}

