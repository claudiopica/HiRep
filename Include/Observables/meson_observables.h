/***************************************************************************\
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen                              *
*                                                                           *
*                                                                           *
\***************************************************************************/
#ifndef MESON_OBSERVABLES
#define MESON_OBSERVABLES

#ifdef __cplusplus
	extern "C" {
#endif


typedef enum { _g5= 0,_id, _g0, _g1, _g2, _g3,  _g0g1, _g0g2, _g0g3, _g0g5, _g5g1, _g5g2, _g5g3, _g0g5g1, _g0g5g2, _g0g5g3,NGAMMA_IND,_DISC} gamma_ind;
extern char *meson_channel_names[NGAMMA_IND]; //defined in measure_mesons.c

typedef struct meson_observable_s {
  gamma_ind ind1;
  gamma_ind ind2;
  char channel_name[100];
  char channel_type[100];
  double sign;
  int corr_size;
  double* corr_re;
  double* corr_im;
  struct meson_observable_s* next;
} meson_observable;

#ifdef __cplusplus
  }
#endif
#endif
