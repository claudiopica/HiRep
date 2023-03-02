#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h>
#include <time.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// High resolution timing functions
typedef struct timespec Instant;

/// @brief returns the high resolution clock time
static inline Instant now() {
    Instant clock;
    clock_gettime(CLOCK_MONOTONIC_RAW, &clock);
    return clock;
}
/// @brief returns the interval in microsecond between start and end
static inline double interval_usec(Instant const *end, Instant const *start) {
    return (end->tv_sec - start->tv_sec) * 1.e6 + (end->tv_nsec - start->tv_nsec) * 1.e-3;
}

/// @brief returns the number of nanoseconds between start and end
static inline uint_fast64_t interval_nsec(Instant const *end, Instant const *start) {
    return (end->tv_sec - start->tv_sec) * (uint_fast64_t)(1000000000U) + (end->tv_nsec - start->tv_nsec);
}

typedef struct Timer {
    Instant start;
    Instant lap;
} Timer;

/// @brief Set the timer to initial time
static inline void timer_set(Timer *t) {
    Instant const n = now();
    t->start = n;
    t->lap = n;
}

/// @brief returns the time in microseconds since the timer was started
static inline double timer_read(Timer const *t) {
    Instant const n = now();
    return interval_usec(&n, &(t->start));
}

/// @brief returns the time in microsecond since the last lap. Restart lap time.
static inline double timer_lap(Timer *t) {
    Instant const l = now();
    double laptime = interval_usec(&l, &(t->lap));
    t->lap = l;
    return laptime;
}

/// @brief returns the nominal timer resolution in microseconds
static inline double timer_res() {
    Instant res;
    clock_getres(CLOCK_MONOTONIC_RAW, &res);
    return res.tv_sec * 1.e6 + res.tv_nsec * 1.e-3;
}

/* old low timing function */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y);

#ifdef __cplusplus
}
#endif
#endif //TIMING_H
