#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h>

#ifdef __cplusplus
	extern "C" {
#endif

/* Timing */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

#ifdef __cplusplus
	}
#endif
#endif //TIMING_H
