#ifndef LOCAL_ACTION_H
#define LOCAL_ACTION_H

#include "spinor_field.h"

#ifdef __cplusplus
extern "C" {
#endif

//local_action.c
/* local action */
typedef enum { NEW = 1, DELTA = 2 } local_action_type;

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
extern void (*local_hmc_action)(local_action_type type, scalar_field *loc_action, suNg_av_field *momenta,
                                suNg_scalar_field *momenta_s);
extern void (*pf_local_action)(scalar_field *loc_action, spinor_field *pf);

void local_hmc_action_cpu(local_action_type type, scalar_field *loc_action, suNg_av_field *momenta,
                          suNg_scalar_field *momenta_s);
void local_hmc_action_gpu(local_action_type type, scalar_field *loc_action, suNg_av_field *momenta,
                          suNg_scalar_field *momenta_s);

void pf_local_action_cpu(scalar_field *loc_action, spinor_field *pf);
void pf_local_action_gpu(scalar_field *loc_action, spinor_field *pf);

#ifdef __cplusplus
}
#endif
#endif //LOCAL_ACTION_H
