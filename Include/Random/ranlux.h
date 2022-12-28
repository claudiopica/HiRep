/// Headerfile for:
/// - ranlxd.c
/// - ranlxs.c

#ifndef RANLUX_H
#define RANLUX_H
#ifdef __cplusplus
	extern "C" {
#endif

/*******************************************************************************
*
* file ranlxd.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

void ranlxd(double r[], int n);
void rlxd_init(int level, int seed);
int rlxd_size(void);
void rlxd_get(int state[]);
void rlxd_reset(int state[]);

/*******************************************************************************
*
* file ranlxs.h
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
*******************************************************************************/

void ranlxs(float r[], int n);
void rlxs_init(int level, int seed);
int rlxs_size(void);
void rlxs_get(int state[]);
void rlxs_reset(int state[]);

#ifdef __cplusplus
    }
#endif
#endif /* RANLUX_H */
