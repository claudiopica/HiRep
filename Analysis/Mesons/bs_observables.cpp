#include <cmath>

#include "bs_observables.h"
#include "effective_mass.h"
#include "bs_ctrl.h"

extern int Lt;
extern int Ls;
extern double svol;

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

void mpcac_eval(Corr_t* mpcac_eff, Corr_t* g5_cor, Corr_t* g5_eff, Corr_t* g5_g0g5_cor) {
  for (int i=mpcac.get_left_cut(); i<=mpcac.get_right_cut(); i++) {
    double ps=g5_eff->d[i].back();
    double m[3];
    m[0]=g5_g0g5_cor->d[i-1].back();
    m[1]=g5_cor->d[i].back();
    if(i<Lt/2) m[2]=g5_g0g5_cor->d[i+1].back();
    else m[2]=-g5_g0g5_cor->d[Lt/2-1].back();
    mpcac_eff->d[i].push_back(-(fabs(ps)>1.e-15?ps/sinh(ps):1.)*(m[0]-m[2])/m[1]/4.);
  }
}


void gps_eval(Corr_t* gps_eff, Corr_t* g5_cor, Corr_t* g5_eff) {
  for (int i=gps.get_left_cut(); i<=gps.get_right_cut(); i++) {
    double mps=g5_eff->d[i].back();
    double cps=g5_cor->d[i].back();
    gps_eff->d[i].push_back(sqrt(mps*cps/hc(i,mps,Lt))*svol);
  }
}

void fps_eval(Corr_t* fps_eff, Corr_t* g5_cor, Corr_t* g5_eff, Corr_t* mpcac_eff) {
  for (int i=fps.get_left_cut(); i<=fps.get_right_cut(); i++) {
    double mps=g5_eff->d[i].back();
    double cps=g5_cor->d[i].back();
    double _mpcac=mpcac_eff->d[i].back();
    fps_eff->d[i].push_back(_mpcac/mps/mps*2.*sqrt(mps*cps/hc(i,mps,Lt))*svol);
  }
}

void fv_eval(Corr_t* fv_eff, Corr_t* g1_cor, Corr_t* g1_eff) {
  for (int i=fv.get_left_cut(); i<=fv.get_right_cut(); i++) {
    double mv=g1_eff->d[i].back();
    double cv=g1_cor->d[i].back();
    fv_eff->d[i].push_back(sqrt(cv/(mv*hc(i,mv,Lt)))*svol);
  }
}

void fvk_eval(Corr_t* fvk_eff, Corr_t* gk_cor, Corr_t* gk_eff) {
  for (int i=fvk.get_left_cut(); i<=fvk.get_right_cut(); i++) {
    double mv=gk_eff->d[i].back();
    double cv=gk_cor->d[i].back();
    fvk_eff->d[i].push_back(sqrt(cv/(mv*hc(i,mv,Lt)))*svol);
  }
}

void fak_eval(Corr_t* fak_eff, Corr_t* g5gk_cor, Corr_t* g5gk_eff) {
	for (int i=fak.get_left_cut(); i<=fak.get_right_cut(); i++) {
		double ma=g5gk_eff->d[i].back();
		double ca=g5gk_cor->d[i].back();
		fak_eff->d[i].push_back(sqrt(ca/(ma*hc(i,ma,Lt)))*svol);
	}
}
