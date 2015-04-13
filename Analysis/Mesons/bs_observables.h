#ifndef __BS_OBSERVABLES
#define __BS_OBSERVABLES
#include "bs_type.h"


void mpcac_eval(Corr_t* mpcac_eff, Corr_t* g5_cor, Corr_t* g5_eff, Corr_t* g5_g0g5_cor);
void gps_eval(Corr_t* gps_eff, Corr_t* g5_cor, Corr_t* g5_eff);
void fps_eval(Corr_t* fps_eff, Corr_t* g5_cor, Corr_t* g5_eff, Corr_t* mpcac_eff);
void fv_eval(Corr_t* fv_eff, Corr_t* g1_cor, Corr_t* g1_eff);
void fvk_eval(Corr_t* fvk_eff, Corr_t* gk_cor, Corr_t* gk_eff);
void fak_eval(Corr_t* fak_eff, Corr_t* g5gk_cor, Corr_t* g5gk_eff);

#endif
