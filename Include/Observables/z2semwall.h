#ifndef Z2SEMWALL_H
#define Z2SEMWALL_H
#ifdef __cplusplus
extern "C" {
#endif

// z2semwall.c
void z2semwall_qprop_free();
void z2semwall_mesons(int conf, int nhits, int nm, double *m, double acc);

// z2semwall_new.c
void z2semwall_qprop_free_new();
void z2semwall_mesons_new(int conf, int nhits, int nm, double *m, double acc);

#ifdef __cplusplus
}
#endif
#endif //Z2SEMWALL_H
