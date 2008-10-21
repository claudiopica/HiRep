/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File statistics.h
* 
* Functions for statistical analysis of data series
*
*******************************************************************************/

#ifndef EXTRAS_H
#define EXTRAS_H

double average(int n,double a[]);
double sigma0(int n,double a[]);
void auto_corr(int n,double a[],int tmax,double gamma[]);
double sigma(int n,double a[],double *tau,int *flag);

double auto_corr_time(int n,int tmax,double g[],int *flag);
double sigma_bin(int n, int binsize, double a[]);
double sigma_replicas(int n,int r,double a[],double *tau,int *flag);
double sigma_jackknife(int nobs,int n,double a[],double *ave_j,double (*pobs)(double v[]));


#endif

 
