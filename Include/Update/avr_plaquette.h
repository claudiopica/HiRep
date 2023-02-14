/**
 * @file avr_plaquette.h
 * @brief Plaquette evaluation functions
 */

#ifndef AVR_PLAQUETTE_H
#define AVR_PLAQUETTE_H

#include "hr_complex.h"

#ifdef __cplusplus
	extern "C" {
#endif

//avr_plaquette.h
double plaq(int ix, int mu, int nu);
void cplaq(hr_complex *ret, int ix, int mu, int nu);
double avr_plaquette(void);
void avr_plaquette_time(double *plaqt, double *plaqs);
double local_plaq(int ix);
void full_plaquette(void);
void avr_ts_plaquette(void);

void cplaq_wrk(hr_complex *ret, int ix, int mu, int nu);
hr_complex avr_plaquette_wrk(void);

//TODO: not defined?
// double rect_1x2(int ix, int mu, int nu);
// void crect_1x2(hr_complex *ret, int ix, int mu, int nu);
// double avr_rect_1x2();
// void full_rect_1x2();
// double local_rect_1x2(int ix);

#ifdef __cplusplus
    }
#endif
#endif
