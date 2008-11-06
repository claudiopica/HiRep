/* This is part of the ranlux 3.3 distribution by Martin Luesher */

#ifndef RANLUX_H
#define RANLUX_H

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

extern void ranlxd(double r[],int n);
extern void rlxd_init(int level,int seed);
extern int rlxd_size(void);
extern void rlxd_get(int state[]);
extern void rlxd_reset(int state[]);

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

extern void ranlxs(float r[],int n);
extern void rlxs_init(int level,int seed);
extern int rlxs_size(void);
extern void rlxs_get(int state[]);
extern void rlxs_reset(int state[]);

#endif /* RANLUX_H */
